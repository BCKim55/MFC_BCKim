!>
!! @file
!! @brief Contains module m_lso_filter

#:include 'macros.fpp'

!> @brief LSO variable-weight Gaussian filter for conserved variables at save steps.
!!
!! 9-point symmetric FIR stencil
!!   H(xi) = a(1) + 2*[a(2)*cos(xi) + a(3)*cos(2*xi) + a(4)*cos(3*xi) + a(5)*cos(4*xi)],
!! composed over several passes to approximate G(xi) = exp(-sigma^2*xi^2/2). Per-pass
!! weights a(1:5) come from the Python BCD design. Applied direction by direction.
module m_lso_filter

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_constants
    use m_ibm, only: ib_markers
    use m_nvtx

    implicit none

    private

    public :: s_initialize_lso_filter_module, s_copy_and_apply_lso_filter, s_lso_stride_sample, s_finalize_lso_filter_module, &
        & q_filt_vf, q_filt_ds_vf, q_lso_mask_vf, q_lso_mask_ds_vf, s_lso_filter_stage2, s_apply_lso_filter_coarse

    ! Floor on the normalized-convolution denominator filter(m).
    real(wp), parameter :: lso_w_floor = 1.0e-3_wp

    ! Scratch buffer for one directional pass (interior only).
    real(wp), allocatable, dimension(:,:,:) :: lso_tmp
    $:GPU_DECLARE(create='[lso_tmp]')

    ! Filtered copy of the conserved variables (lso_filter_wrt = T).
    type(scalar_field), allocatable :: q_filt_vf(:)

    ! Gas mask (1 in fluid, 0 in solid), the denominator of the normalized convolution filter(m*q)/filter(m).
    type(scalar_field), allocatable :: q_lso_mask_vf(:)

    ! Coarsened version (lso_down_sample_factor > 1).
    type(scalar_field), allocatable :: q_filt_ds_vf(:)

    ! Coarse-grid (stage-2) arrays carry lso_crs_gw ghost layers per active direction; stage 2 runs on the host.
    integer, parameter    :: lso_crs_gw = 4
    integer               :: crs_lo(3), crs_hi(3)
    real(wp), allocatable :: lso2_tmp(:,:,:)
#ifdef MFC_MPI
    real(wp), allocatable :: buff_crs_send(:), buff_crs_recv(:)
#endif

    ! Coarsened filtered gas mask w = filter(m), written so post_process can compose filter2(w*qhat)/filter2(w).
    type(scalar_field), allocatable :: q_lso_mask_ds_vf(:)

contains

    impure subroutine s_initialize_lso_filter_module()

        integer :: i

        @:ALLOCATE(lso_tmp(0:m, 0:n, 0:p))

        if (lso_filter_wrt) then
            @:ALLOCATE(q_filt_vf(1:sys_size))
            do i = 1, sys_size
                @:ALLOCATE(q_filt_vf(i)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
            end do
            do i = 1, sys_size
                @:ACC_SETUP_SFs(q_filt_vf(i))
            end do

            ! Gas-mask denominator field for IB-aware normalized filtering.
            if (ib) then
                @:ALLOCATE(q_lso_mask_vf(1:1))
                @:ALLOCATE(q_lso_mask_vf(1)%sf(idwbuff(1)%beg:idwbuff(1)%end, idwbuff(2)%beg:idwbuff(2)%end, &
                           & idwbuff(3)%beg:idwbuff(3)%end))
                @:ACC_SETUP_SFs(q_lso_mask_vf(1))
            end if

            if (lso_down_sample_factor > 1) then
                crs_lo = 0; crs_hi = 0
                crs_lo(1) = -lso_crs_gw; crs_hi(1) = m_lso_ds + lso_crs_gw
                if (n_lso_ds > 0) then
                    crs_lo(2) = -lso_crs_gw; crs_hi(2) = n_lso_ds + lso_crs_gw
                end if
                if (p_lso_ds > 0) then
                    crs_lo(3) = -lso_crs_gw; crs_hi(3) = p_lso_ds + lso_crs_gw
                end if
                @:ALLOCATE(q_filt_ds_vf(1:sys_size))
                do i = 1, sys_size
                    @:ALLOCATE(q_filt_ds_vf(i)%sf(crs_lo(1):crs_hi(1), crs_lo(2):crs_hi(2), crs_lo(3):crs_hi(3)))
                end do
                if (ib) then
                    @:ALLOCATE(q_lso_mask_ds_vf(1:1))
                    @:ALLOCATE(q_lso_mask_ds_vf(1)%sf(crs_lo(1):crs_hi(1), crs_lo(2):crs_hi(2), crs_lo(3):crs_hi(3)))
                end if
                if (lso2_n_passes_x > 0) then
                    allocate (lso2_tmp(0:m_lso_ds,0:n_lso_ds,0:p_lso_ds))
#ifdef MFC_MPI
                    block
                        integer :: nv_max, face_max
                        nv_max = sys_size
                        face_max = max((n_lso_ds + 1)*(p_lso_ds + 1), (m_lso_ds + 1)*(p_lso_ds + 1), (m_lso_ds + 1)*(n_lso_ds + 1))
                        allocate (buff_crs_send(0:nv_max*lso_crs_gw*face_max - 1))
                        allocate (buff_crs_recv(0:nv_max*lso_crs_gw*face_max - 1))
                    end block
#endif
                end if
            end if
        end if

    end subroutine s_initialize_lso_filter_module

    impure subroutine s_finalize_lso_filter_module()

        integer :: i

        @:DEALLOCATE(lso_tmp)

        if (lso_filter_wrt) then
            do i = 1, sys_size
                @:DEALLOCATE(q_filt_vf(i)%sf)
            end do
            @:DEALLOCATE(q_filt_vf)

            if (ib) then
                @:DEALLOCATE(q_lso_mask_vf(1)%sf)
                @:DEALLOCATE(q_lso_mask_vf)
            end if

            if (lso_down_sample_factor > 1) then
                do i = 1, sys_size
                    @:DEALLOCATE(q_filt_ds_vf(i)%sf)
                end do
                @:DEALLOCATE(q_filt_ds_vf)
                if (ib) then
                    @:DEALLOCATE(q_lso_mask_ds_vf(1)%sf)
                    @:DEALLOCATE(q_lso_mask_ds_vf)
                end if
                if (allocated(lso2_tmp)) deallocate (lso2_tmp)
#ifdef MFC_MPI
                if (allocated(buff_crs_send)) deallocate (buff_crs_send, buff_crs_recv)
#endif
            end if
        end if

    end subroutine s_finalize_lso_filter_module

    !> Copy q_cons_vf into q_filt_vf on the device and filter in place, leaving the original conserved array untouched for the
    !! primary write.
    impure subroutine s_copy_and_apply_lso_filter(q_cons_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer                           :: i, j, k, l

        call nvtxStartRange("LSO-FILTER")

        call nvtxStartRange("LSO-FILTER-COPY")
        do i = 1, sys_size
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = idwbuff(3)%beg, idwbuff(3)%end
                do k = idwbuff(2)%beg, idwbuff(2)%end
                    do j = idwbuff(1)%beg, idwbuff(1)%end
                        q_filt_vf(i)%sf(j, k, l) = q_cons_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
        end do
        call nvtxEndRange

        if (ib) then
            ! Normalized convolution q_filt = filter(m*q)/filter(m) with the gas mask m (1 in fluid, 0 in solid).
            call nvtxStartRange("LSO-FILTER-MASK")
            $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        if (ib_markers%sf(j, k, l) > 0) then
                            q_lso_mask_vf(1)%sf(j, k, l) = 0._stp
                        else
                            q_lso_mask_vf(1)%sf(j, k, l) = 1._stp
                        end if
                    end do
                end do
            end do
            $:END_GPU_PARALLEL_LOOP()
            do i = 1, sys_size
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_filt_vf(i)%sf(j, k, l) = real(real(q_filt_vf(i)%sf(j, k, l), wp)*real(q_lso_mask_vf(1)%sf(j, k, l), &
                                      & wp), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            call nvtxEndRange

            call s_apply_lso_filter(q_filt_vf)
            call s_apply_lso_filter(q_lso_mask_vf)

            ! Normalize everywhere (the solid interior takes the fluid average). Under the two-stage
            ! pyramid the numerator and denominator stay unnormalized until s_lso_filter_stage2.
            if (lso2_n_passes_x <= 0) then
                do i = 1, sys_size
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_filt_vf(i)%sf(j, k, l) = real(real(q_filt_vf(i)%sf(j, k, l), &
                                          & wp)/max(real(q_lso_mask_vf(1)%sf(j, k, l), wp), lso_w_floor), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
        else
            call s_apply_lso_filter(q_filt_vf)
        end if

        call nvtxEndRange

    end subroutine s_copy_and_apply_lso_filter

    !> Apply the LSO filter in place to q_cons_vf, direction by direction. Passes loop outside variables so a halo exchange can run
    !! between passes.
    !!
    !! @note With forward Euler (stor == 1) this modifies the live state; with any RK
    !! scheme (stor == 2) only the save copy is touched. IBM runs use RK.
    impure subroutine s_apply_lso_filter(q_cons_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer                           :: i, ipass, j, k, l, nv
        real(wp)                          :: c0, c1, c2, c3, c4

        nv = size(q_cons_vf)

        call nvtxStartRange("LSO-FILTER-X")
        call s_lso_filter_ghost_refresh(q_cons_vf, 1)
        do ipass = 1, lso_n_passes_x
            c0 = lso_a_x(1, ipass)
            c1 = lso_a_x(2, ipass)
            c2 = lso_a_x(3, ipass)
            c3 = lso_a_x(4, ipass)
            c4 = lso_a_x(5, ipass)
            do i = 1, nv
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j - 1, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_cons_vf(i)%sf(j - 2, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 2, k, l), wp)) + c3*(real(q_cons_vf(i)%sf(j - 3, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 3, k, l), wp)) + c4*(real(q_cons_vf(i)%sf(j - 4, k, l), &
                                    & wp) + real(q_cons_vf(i)%sf(j + 4, k, l), wp))
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
                $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end do
            if (ipass < lso_n_passes_x) then
                call s_lso_filter_ghost_refresh(q_cons_vf, 1)
            end if
        end do
        call nvtxEndRange

        ! y-direction (2D/3D)
        if (n > 0) then
            call nvtxStartRange("LSO-FILTER-Y")
            call s_lso_filter_ghost_refresh(q_cons_vf, 2)
            do ipass = 1, lso_n_passes_y
                c0 = lso_a_y(1, ipass)
                c1 = lso_a_y(2, ipass)
                c2 = lso_a_y(3, ipass)
                c3 = lso_a_y(4, ipass)
                c4 = lso_a_y(5, ipass)
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k - 1, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 1, l), wp)) + c2*(real(q_cons_vf(i)%sf(j, k - 2, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 2, l), wp)) + c3*(real(q_cons_vf(i)%sf(j, k - 3, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 3, l), wp)) + c4*(real(q_cons_vf(i)%sf(j, k - 4, l), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k + 4, l), wp))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                if (ipass < lso_n_passes_y) then
                    call s_lso_filter_ghost_refresh(q_cons_vf, 2)
                end if
            end do
            call nvtxEndRange
        end if

        ! z-direction (3D)
        if (p > 0) then
            call nvtxStartRange("LSO-FILTER-Z")
            call s_lso_filter_ghost_refresh(q_cons_vf, 3)
            do ipass = 1, lso_n_passes_z
                c0 = lso_a_z(1, ipass)
                c1 = lso_a_z(2, ipass)
                c2 = lso_a_z(3, ipass)
                c3 = lso_a_z(4, ipass)
                c4 = lso_a_z(5, ipass)
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k, l - 1), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 1), wp)) + c2*(real(q_cons_vf(i)%sf(j, k, l - 2), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 2), wp)) + c3*(real(q_cons_vf(i)%sf(j, k, l - 3), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 3), wp)) + c4*(real(q_cons_vf(i)%sf(j, k, l - 4), &
                                        & wp) + real(q_cons_vf(i)%sf(j, k, l + 4), wp))
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
                if (ipass < lso_n_passes_z) then
                    call s_lso_filter_ghost_refresh(q_cons_vf, 3)
                end if
            end do
            call nvtxEndRange
        end if

    end subroutine s_apply_lso_filter

    !> Resample q_src_vf into the (pre-allocated) q_dst_vf onto the coarsened grid using trilinear interpolation. Output cell j maps
    !! to source position j*m/m_lso_ds (uniformly spaced from 0 to m), so the first and last cells are always exact and the domain
    !! is fully preserved regardless of divisibility by the stride factor.
    impure subroutine s_lso_stride_sample(q_src_vf, q_dst_vf)

        type(scalar_field), intent(in)    :: q_src_vf(:)
        type(scalar_field), intent(inout) :: q_dst_vf(:)
        integer                           :: i, j, k, l, nv
        integer                           :: j0, j1, k0, k1, l0, l1
        integer                           :: sidx(3)
        real(wp)                          :: alpha, beta, gamma_pos, wj, wk, wl

        ! Global offset of this rank's block (start_idx is only allocated with parallel_io).

        sidx = 0
        if (allocated(start_idx)) sidx(1:size(start_idx)) = start_idx

        nv = size(q_src_vf)
        do i = 1, nv
            do l = 0, p_lso_ds
                ! Sample at the coarse-cell centres: fine position (J + 1/2)(m_glb + 1)/(m_glb_ds + 1) - 1/2 for the
                ! global coarse index J = start_idx/factor + j, so the sample spacing equals the coarse cell size
                ! (interior-only reads for factor >= 2).
                if (p_lso_ds > 0) then
                    gamma_pos = (real(sidx(3)/lso_down_sample_factor + l, wp) + 0.5_wp)*real(p_glb + 1, &
                                 & wp)/real(p_glb_lso_ds + 1, wp) - 0.5_wp - real(sidx(3), wp)
                    l0 = floor(gamma_pos); l1 = l0 + 1; wl = gamma_pos - real(l0, wp)
                else
                    l0 = 0; l1 = 0; wl = 0._wp
                end if

                do k = 0, n_lso_ds
                    if (n_lso_ds > 0) then
                        beta = (real(sidx(2)/lso_down_sample_factor + k, wp) + 0.5_wp)*real(n_glb + 1, wp)/real(n_glb_lso_ds + 1, &
                                & wp) - 0.5_wp - real(sidx(2), wp)
                        k0 = floor(beta); k1 = k0 + 1; wk = beta - real(k0, wp)
                    else
                        k0 = 0; k1 = 0; wk = 0._wp
                    end if

                    do j = 0, m_lso_ds
                        if (m_lso_ds > 0) then
                            alpha = (real(sidx(1)/lso_down_sample_factor + j, wp) + 0.5_wp)*real(m_glb + 1, &
                                     & wp)/real(m_glb_lso_ds + 1, wp) - 0.5_wp - real(sidx(1), wp)
                            j0 = floor(alpha); j1 = j0 + 1; wj = alpha - real(j0, wp)
                        else
                            j0 = 0; j1 = 0; wj = 0._wp
                        end if

                        q_dst_vf(i)%sf(j, k, l) = real((1._wp - wj)*(1._wp - wk)*(1._wp - wl)*real(q_src_vf(i)%sf(j0, k0, l0), &
                                 & wp) + wj*(1._wp - wk)*(1._wp - wl)*real(q_src_vf(i)%sf(j1, k0, l0), &
                                 & wp) + (1._wp - wj)*wk*(1._wp - wl)*real(q_src_vf(i)%sf(j0, k1, l0), &
                                 & wp) + wj*wk*(1._wp - wl)*real(q_src_vf(i)%sf(j1, k1, l0), &
                                 & wp) + (1._wp - wj)*(1._wp - wk)*wl*real(q_src_vf(i)%sf(j0, k0, l1), &
                                 & wp) + wj*(1._wp - wk)*wl*real(q_src_vf(i)%sf(j1, k0, l1), &
                                 & wp) + (1._wp - wj)*wk*wl*real(q_src_vf(i)%sf(j0, k1, l1), &
                                 & wp) + wj*wk*wl*real(q_src_vf(i)%sf(j1, k1, l1), wp), stp)
                    end do
                end do
            end do
        end do

    end subroutine s_lso_stride_sample

    !> Host halo refresh for the coarse (decimated) grid used by the stage-2 pyramid cascade: MPI exchange of lso_crs_gw layers on
    !! decomposed directions (sequential counter packing, identical nesting on pack and unpack), then periodic wrap or edge clamp
    !! for physical boundaries.
    impure subroutine s_lso_coarse_ghost_refresh(q_vf, mpi_dir)

        type(scalar_field), intent(inout) :: q_vf(:)
        integer, intent(in)               :: mpi_dir
        integer                           :: nv, i, j, k, l, g
        integer                           :: beg_bc, end_bc, grid_dim

#ifdef MFC_MPI
        integer :: ierr, r, cnt, pack_offset, unpack_offset, pbc_loc, side
        integer :: dst_proc, src_proc, send_tag, recv_tag
        integer :: beg_end(2)
        logical :: beg_end_geq_0
#endif

        nv = size(q_vf)
        select case (mpi_dir)
        case (1)
            beg_bc = bc_x%beg; end_bc = bc_x%end; grid_dim = m_lso_ds
        case (2)
            beg_bc = bc_y%beg; end_bc = bc_y%end; grid_dim = n_lso_ds
        case (3)
            beg_bc = bc_z%beg; end_bc = bc_z%end; grid_dim = p_lso_ds
        end select

#ifdef MFC_MPI
        beg_end = (/beg_bc, end_bc/)
        select case (mpi_dir)
        case (1)
            cnt = nv*lso_crs_gw*(n_lso_ds + 1)*(p_lso_ds + 1)
        case (2)
            cnt = nv*lso_crs_gw*(m_lso_ds + 1)*(p_lso_ds + 1)
        case (3)
            cnt = nv*lso_crs_gw*(m_lso_ds + 1)*(n_lso_ds + 1)
        end select

        do side = 1, 2
            pbc_loc = 2*side - 3  ! -1 (beg) then +1 (end)
            if (beg_end(side) < 0) cycle
            beg_end_geq_0 = beg_end(max(pbc_loc, 0) - pbc_loc + 1) >= 0
            send_tag = f_logical_to_int(.not. f_xor(beg_end_geq_0, pbc_loc == 1))
            recv_tag = f_logical_to_int(pbc_loc == 1)
            dst_proc = beg_end(1 + f_logical_to_int(f_xor(pbc_loc == 1, beg_end_geq_0)))
            src_proc = beg_end(1 + f_logical_to_int(pbc_loc == 1))
            pack_offset = 0
            if (f_xor(pbc_loc == 1, beg_end_geq_0)) pack_offset = grid_dim - lso_crs_gw + 1
            unpack_offset = 0
            if (pbc_loc == 1) unpack_offset = grid_dim + lso_crs_gw + 1

            r = -1
            select case (mpi_dir)
            case (1)
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do g = 0, lso_crs_gw - 1
                            do i = 1, nv
                                r = r + 1
                                buff_crs_send(r) = real(q_vf(i)%sf(g + pack_offset, k, l), wp)
                            end do
                        end do
                    end do
                end do
            case (2)
                do l = 0, p_lso_ds
                    do g = 0, lso_crs_gw - 1
                        do j = 0, m_lso_ds
                            do i = 1, nv
                                r = r + 1
                                buff_crs_send(r) = real(q_vf(i)%sf(j, g + pack_offset, l), wp)
                            end do
                        end do
                    end do
                end do
            case (3)
                do g = 0, lso_crs_gw - 1
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            do i = 1, nv
                                r = r + 1
                                buff_crs_send(r) = real(q_vf(i)%sf(j, k, g + pack_offset), wp)
                            end do
                        end do
                    end do
                end do
            end select

            call MPI_SENDRECV(buff_crs_send, cnt, mpi_p, dst_proc, send_tag, buff_crs_recv, cnt, mpi_p, src_proc, recv_tag, &
                              & MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

            r = -1
            select case (mpi_dir)
            case (1)
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do g = -lso_crs_gw, -1
                            do i = 1, nv
                                r = r + 1
                                q_vf(i)%sf(g + unpack_offset, k, l) = real(buff_crs_recv(r), stp)
                            end do
                        end do
                    end do
                end do
            case (2)
                do l = 0, p_lso_ds
                    do g = -lso_crs_gw, -1
                        do j = 0, m_lso_ds
                            do i = 1, nv
                                r = r + 1
                                q_vf(i)%sf(j, g + unpack_offset, l) = real(buff_crs_recv(r), stp)
                            end do
                        end do
                    end do
                end do
            case (3)
                do g = -lso_crs_gw, -1
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            do i = 1, nv
                                r = r + 1
                                q_vf(i)%sf(j, k, g + unpack_offset) = real(buff_crs_recv(r), stp)
                            end do
                        end do
                    end do
                end do
            end select
        end do
#endif

        ! Physical boundaries: periodic wrap, else clamp to the edge cell.
        do i = 1, nv
            select case (mpi_dir)
            case (1)
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do g = 1, lso_crs_gw
                            if (beg_bc == BC_PERIODIC) then
                                q_vf(i)%sf(-g, k, l) = q_vf(i)%sf(grid_dim - g + 1, k, l)
                            else if (beg_bc < 0) then
                                q_vf(i)%sf(-g, k, l) = q_vf(i)%sf(0, k, l)
                            end if
                            if (end_bc == BC_PERIODIC) then
                                q_vf(i)%sf(grid_dim + g, k, l) = q_vf(i)%sf(g - 1, k, l)
                            else if (end_bc < 0) then
                                q_vf(i)%sf(grid_dim + g, k, l) = q_vf(i)%sf(grid_dim, k, l)
                            end if
                        end do
                    end do
                end do
            case (2)
                do l = 0, p_lso_ds
                    do g = 1, lso_crs_gw
                        do j = 0, m_lso_ds
                            if (beg_bc == BC_PERIODIC) then
                                q_vf(i)%sf(j, -g, l) = q_vf(i)%sf(j, grid_dim - g + 1, l)
                            else if (beg_bc < 0) then
                                q_vf(i)%sf(j, -g, l) = q_vf(i)%sf(j, 0, l)
                            end if
                            if (end_bc == BC_PERIODIC) then
                                q_vf(i)%sf(j, grid_dim + g, l) = q_vf(i)%sf(j, g - 1, l)
                            else if (end_bc < 0) then
                                q_vf(i)%sf(j, grid_dim + g, l) = q_vf(i)%sf(j, grid_dim, l)
                            end if
                        end do
                    end do
                end do
            case (3)
                do g = 1, lso_crs_gw
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            if (beg_bc == BC_PERIODIC) then
                                q_vf(i)%sf(j, k, -g) = q_vf(i)%sf(j, k, grid_dim - g + 1)
                            else if (beg_bc < 0) then
                                q_vf(i)%sf(j, k, -g) = q_vf(i)%sf(j, k, 0)
                            end if
                            if (end_bc == BC_PERIODIC) then
                                q_vf(i)%sf(j, k, grid_dim + g) = q_vf(i)%sf(j, k, g - 1)
                            else if (end_bc < 0) then
                                q_vf(i)%sf(j, k, grid_dim + g) = q_vf(i)%sf(j, k, grid_dim)
                            end if
                        end do
                    end do
                end do
            end select
        end do

    end subroutine s_lso_coarse_ghost_refresh

    !> Stage-2 cascade of the two-stage in-situ pyramid: apply the lso2_a_* passes to the coarse-grid fields in place (host).
    impure subroutine s_apply_lso_filter_coarse(q_vf)

        type(scalar_field), intent(inout) :: q_vf(:)
        integer                           :: i, ipass, j, k, l, nv
        real(wp)                          :: c0, c1, c2, c3, c4

        nv = size(q_vf)

        call s_lso_coarse_ghost_refresh(q_vf, 1)
        do ipass = 1, lso2_n_passes_x
            c0 = lso2_a_x(1, ipass); c1 = lso2_a_x(2, ipass); c2 = lso2_a_x(3, ipass)
            c3 = lso2_a_x(4, ipass); c4 = lso2_a_x(5, ipass)
            do i = 1, nv
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            lso2_tmp(j, k, l) = c0*real(q_vf(i)%sf(j, k, l), wp) + c1*(real(q_vf(i)%sf(j - 1, k, l), &
                                     & wp) + real(q_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_vf(i)%sf(j - 2, k, l), &
                                     & wp) + real(q_vf(i)%sf(j + 2, k, l), wp)) + c3*(real(q_vf(i)%sf(j - 3, k, l), &
                                     & wp) + real(q_vf(i)%sf(j + 3, k, l), wp)) + c4*(real(q_vf(i)%sf(j - 4, k, l), &
                                     & wp) + real(q_vf(i)%sf(j + 4, k, l), wp))
                        end do
                    end do
                end do
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            q_vf(i)%sf(j, k, l) = real(lso2_tmp(j, k, l), stp)
                        end do
                    end do
                end do
            end do
            if (ipass < lso2_n_passes_x) call s_lso_coarse_ghost_refresh(q_vf, 1)
        end do

        if (n_lso_ds > 0) then
            call s_lso_coarse_ghost_refresh(q_vf, 2)
            do ipass = 1, lso2_n_passes_y
                c0 = lso2_a_y(1, ipass); c1 = lso2_a_y(2, ipass); c2 = lso2_a_y(3, ipass)
                c3 = lso2_a_y(4, ipass); c4 = lso2_a_y(5, ipass)
                do i = 1, nv
                    do l = 0, p_lso_ds
                        do k = 0, n_lso_ds
                            do j = 0, m_lso_ds
                                lso2_tmp(j, k, l) = c0*real(q_vf(i)%sf(j, k, l), wp) + c1*(real(q_vf(i)%sf(j, k - 1, l), &
                                         & wp) + real(q_vf(i)%sf(j, k + 1, l), wp)) + c2*(real(q_vf(i)%sf(j, k - 2, l), &
                                         & wp) + real(q_vf(i)%sf(j, k + 2, l), wp)) + c3*(real(q_vf(i)%sf(j, k - 3, l), &
                                         & wp) + real(q_vf(i)%sf(j, k + 3, l), wp)) + c4*(real(q_vf(i)%sf(j, k - 4, l), &
                                         & wp) + real(q_vf(i)%sf(j, k + 4, l), wp))
                            end do
                        end do
                    end do
                    do l = 0, p_lso_ds
                        do k = 0, n_lso_ds
                            do j = 0, m_lso_ds
                                q_vf(i)%sf(j, k, l) = real(lso2_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < lso2_n_passes_y) call s_lso_coarse_ghost_refresh(q_vf, 2)
            end do
        end if

        if (p_lso_ds > 0) then
            call s_lso_coarse_ghost_refresh(q_vf, 3)
            do ipass = 1, lso2_n_passes_z
                c0 = lso2_a_z(1, ipass); c1 = lso2_a_z(2, ipass); c2 = lso2_a_z(3, ipass)
                c3 = lso2_a_z(4, ipass); c4 = lso2_a_z(5, ipass)
                do i = 1, nv
                    do l = 0, p_lso_ds
                        do k = 0, n_lso_ds
                            do j = 0, m_lso_ds
                                lso2_tmp(j, k, l) = c0*real(q_vf(i)%sf(j, k, l), wp) + c1*(real(q_vf(i)%sf(j, k, l - 1), &
                                         & wp) + real(q_vf(i)%sf(j, k, l + 1), wp)) + c2*(real(q_vf(i)%sf(j, k, l - 2), &
                                         & wp) + real(q_vf(i)%sf(j, k, l + 2), wp)) + c3*(real(q_vf(i)%sf(j, k, l - 3), &
                                         & wp) + real(q_vf(i)%sf(j, k, l + 3), wp)) + c4*(real(q_vf(i)%sf(j, k, l - 4), &
                                         & wp) + real(q_vf(i)%sf(j, k, l + 4), wp))
                            end do
                        end do
                    end do
                    do l = 0, p_lso_ds
                        do k = 0, n_lso_ds
                            do j = 0, m_lso_ds
                                q_vf(i)%sf(j, k, l) = real(lso2_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < lso2_n_passes_z) call s_lso_coarse_ghost_refresh(q_vf, 3)
            end do
        end if

    end subroutine s_apply_lso_filter_coarse

    !> Complete the two-stage in-situ pyramid on the downsampled arrays: stage-2 filter the decimated numerator (and mask
    !! denominator when ib), then normalize. No-op unless lso2 passes are configured.
    impure subroutine s_lso_filter_stage2()

        integer :: i, j, k, l

        if (lso2_n_passes_x <= 0) return

        call s_apply_lso_filter_coarse(q_filt_ds_vf)
        if (ib) then
            call s_apply_lso_filter_coarse(q_lso_mask_ds_vf)
            do i = 1, sys_size
                do l = 0, p_lso_ds
                    do k = 0, n_lso_ds
                        do j = 0, m_lso_ds
                            q_filt_ds_vf(i)%sf(j, k, l) = real(real(q_filt_ds_vf(i)%sf(j, k, l), &
                                         & wp)/max(real(q_lso_mask_ds_vf(1)%sf(j, k, l), wp), lso_w_floor), stp)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine s_lso_filter_stage2

    !> Refresh ghost cells between filter passes for direction mpi_dir (1=x, 2=y, 3=z). MPI ghosts come via sendrecv, then
    !! IBM-flagged ghosts are zero-extrapolated to keep particle velocity out of the fluid stencil. BC_GHOST_EXTRAP faces are
    !! re-extrapolated from the current edge cell; other physical BCs are left as is.
    impure subroutine s_lso_filter_ghost_refresh(q_cons_vf, mpi_dir)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer, intent(in)               :: mpi_dir
        integer                           :: i, j, k, l, beg_bc, end_bc, nv

        nv = size(q_cons_vf)

        select case (mpi_dir)
        case (1)
            beg_bc = bc_x%beg
            end_bc = bc_x%end
        case (2)
            beg_bc = bc_y%beg
            end_bc = bc_y%end
        case (3)
            beg_bc = bc_z%beg
            end_bc = bc_z%end
        end select

        ! MPI rank boundaries
#ifdef MFC_MPI
        if (beg_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, -1, nv)
        if (end_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, 1, nv)
#endif

        ! BC_GHOST_EXTRAP: re-extrapolate from the filtered edge cell each pass.
        select case (mpi_dir)
        case (1)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(0, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            else if (beg_bc == BC_PERIODIC) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(m - j + 1, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(m, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            else if (end_bc == BC_PERIODIC) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 0, n
                            do j = 1, buff_size
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(j - 1, k, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
#ifdef MFC_MPI
#endif
        case (2)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, 0, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            else if (beg_bc == BC_PERIODIC) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, n - k + 1, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, n, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            else if (end_bc == BC_PERIODIC) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 0, p
                        do k = 1, buff_size
                            do j = 0, m
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, k - 1, l)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
#ifdef MFC_MPI
#endif
        case (3)
            if (beg_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, 0)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            else if (beg_bc == BC_PERIODIC) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, p - l + 1)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
            if (end_bc == BC_GHOST_EXTRAP) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, p)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            else if (end_bc == BC_PERIODIC) then
                do i = 1, nv
                    $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l]')
                    do l = 1, buff_size
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, l - 1)
                            end do
                        end do
                    end do
                    $:END_GPU_PARALLEL_LOOP()
                end do
            end if
#ifdef MFC_MPI
#endif
        end select

    end subroutine s_lso_filter_ghost_refresh

end module m_lso_filter
