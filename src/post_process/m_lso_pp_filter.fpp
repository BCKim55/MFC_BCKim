!>
!! @file
!! @brief Contains module m_lso_pp_filter

#:include 'macros.fpp'

!> @brief Least-squares optimized (LSO) Gaussian filter applied by post_process.
!!
!! Applies repeated symmetric 9-point FIR passes whose composed transfer function
!! approximates a Gaussian of standard deviation lso_filter_sigma_target. The per-pass
!! stencil weights lso_pp_a_* come from the Python block-coordinate-descent design
!! (toolchain/mfc/lso_filter.py) via post_process.inp.
!!
!! With immersed boundaries the filter is the mask-normalized convolution
!! q <- filter(w*q)/max(filter(w), floor) with w the binary gas mask built from
!! ib_markers, so the non-physical solid-interior state never enters the fluid average.
module m_lso_pp_filter

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_constants
    use m_data_input, only: ib_markers

    implicit none

    private

    public :: s_initialize_lso_pp_filter_module, s_finalize_lso_pp_filter_module, s_apply_lso_pp_filter, &
        & s_apply_lso_pp_filter_masked, s_lso_pp_mask_from_ib

    !> Floor on the normalized-convolution denominator filter(w)
    real(wp), parameter :: lso_w_floor = 1.0e-3_wp

    !> Scratch buffer for one directional pass
    real(wp), allocatable :: lso_pp_tmp(:,:,:)

contains

    impure subroutine s_initialize_lso_pp_filter_module()

        allocate (lso_pp_tmp(0:m,0:n,0:p))

    end subroutine s_initialize_lso_pp_filter_module

    impure subroutine s_finalize_lso_pp_filter_module()

        if (allocated(lso_pp_tmp)) deallocate (lso_pp_tmp)

    end subroutine s_finalize_lso_pp_filter_module

    !> Apply the LSO filter to each field of q_cons_vf in place (interior cells), sweeping direction by direction with a ghost
    !! refresh between passes.
    !! @param q_cons_vf Fields to filter
    impure subroutine s_apply_lso_pp_filter(q_cons_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:)
        integer                           :: i, ipass, j, k, l, nv
        real(wp)                          :: c0, c1, c2, c3, c4

        nv = size(q_cons_vf)

        ! x-direction
        call s_lso_pp_filter_ghost_refresh(q_cons_vf, 1)
        do ipass = 1, lso_pp_n_passes_x
            c0 = lso_pp_a_x(1, ipass)
            c1 = lso_pp_a_x(2, ipass)
            c2 = lso_pp_a_x(3, ipass)
            c3 = lso_pp_a_x(4, ipass)
            c4 = lso_pp_a_x(5, ipass)
            do i = 1, nv
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j - 1, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 1, k, l), wp)) + c2*(real(q_cons_vf(i)%sf(j - 2, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 2, k, l), wp)) + c3*(real(q_cons_vf(i)%sf(j - 3, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 3, k, l), wp)) + c4*(real(q_cons_vf(i)%sf(j - 4, k, l), &
                                       & wp) + real(q_cons_vf(i)%sf(j + 4, k, l), wp))
                        end do
                    end do
                end do
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                        end do
                    end do
                end do
            end do
            if (ipass < lso_pp_n_passes_x) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 1)
        end do

        ! y-direction (2D/3D)
        if (n > 0) then
            call s_lso_pp_filter_ghost_refresh(q_cons_vf, 2)
            do ipass = 1, lso_pp_n_passes_y
                c0 = lso_pp_a_y(1, ipass)
                c1 = lso_pp_a_y(2, ipass)
                c2 = lso_pp_a_y(3, ipass)
                c3 = lso_pp_a_y(4, ipass)
                c4 = lso_pp_a_y(5, ipass)
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k - 1, &
                                           & l), wp) + real(q_cons_vf(i)%sf(j, k + 1, l), wp)) + c2*(real(q_cons_vf(i)%sf(j, &
                                           & k - 2, l), wp) + real(q_cons_vf(i)%sf(j, k + 2, l), &
                                           & wp)) + c3*(real(q_cons_vf(i)%sf(j, k - 3, l), wp) + real(q_cons_vf(i)%sf(j, k + 3, &
                                           & l), wp)) + c4*(real(q_cons_vf(i)%sf(j, k - 4, l), wp) + real(q_cons_vf(i)%sf(j, &
                                           & k + 4, l), wp))
                            end do
                        end do
                    end do
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < lso_pp_n_passes_y) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 2)
            end do
        end if

        ! z-direction (3D)
        if (p > 0) then
            call s_lso_pp_filter_ghost_refresh(q_cons_vf, 3)
            do ipass = 1, lso_pp_n_passes_z
                c0 = lso_pp_a_z(1, ipass)
                c1 = lso_pp_a_z(2, ipass)
                c2 = lso_pp_a_z(3, ipass)
                c3 = lso_pp_a_z(4, ipass)
                c4 = lso_pp_a_z(5, ipass)
                do i = 1, nv
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                lso_pp_tmp(j, k, l) = c0*real(q_cons_vf(i)%sf(j, k, l), wp) + c1*(real(q_cons_vf(i)%sf(j, k, &
                                           & l - 1), wp) + real(q_cons_vf(i)%sf(j, k, l + 1), wp)) + c2*(real(q_cons_vf(i)%sf(j, &
                                           & k, l - 2), wp) + real(q_cons_vf(i)%sf(j, k, l + 2), &
                                           & wp)) + c3*(real(q_cons_vf(i)%sf(j, k, l - 3), wp) + real(q_cons_vf(i)%sf(j, k, &
                                           & l + 3), wp)) + c4*(real(q_cons_vf(i)%sf(j, k, l - 4), wp) + real(q_cons_vf(i)%sf(j, &
                                           & k, l + 4), wp))
                            end do
                        end do
                    end do
                    do l = 0, p
                        do k = 0, n
                            do j = 0, m
                                q_cons_vf(i)%sf(j, k, l) = real(lso_pp_tmp(j, k, l), stp)
                            end do
                        end do
                    end do
                end do
                if (ipass < lso_pp_n_passes_z) call s_lso_pp_filter_ghost_refresh(q_cons_vf, 3)
            end do
        end if

    end subroutine s_apply_lso_pp_filter

    !> Fill the interior of w_vf with the binary gas mask from ib_markers: 1 in fluid, 0 inside an immersed body.
    !! @param w_vf Mask field to fill
    impure subroutine s_lso_pp_mask_from_ib(w_vf)

        type(scalar_field), intent(inout) :: w_vf(1:1)
        integer                           :: j, k, l

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    if (ib_markers%sf(j, k, l) > 0) then
                        w_vf(1)%sf(j, k, l) = 0._stp
                    else
                        w_vf(1)%sf(j, k, l) = 1._stp
                    end if
                end do
            end do
        end do

    end subroutine s_lso_pp_mask_from_ib

    !> Mask-normalized filtering: q <- filter(w*q)/max(filter(w), floor), applied in place (the solid interior takes the normalized
    !! fluid average).
    !! @param q_cons_vf Fields to filter
    !! @param w_vf Binary gas mask; holds filter(w) on return
    impure subroutine s_apply_lso_pp_filter_masked(q_cons_vf, w_vf)

        type(scalar_field), intent(inout) :: q_cons_vf(:), w_vf(1:1)
        integer                           :: i, j, k, l, nv

        nv = size(q_cons_vf)
        do i = 1, nv
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = real(real(q_cons_vf(i)%sf(j, k, l), wp)*real(w_vf(1)%sf(j, k, l), wp), stp)
                    end do
                end do
            end do
        end do

        call s_apply_lso_pp_filter(q_cons_vf)
        call s_apply_lso_pp_filter(w_vf)

        do i = 1, nv
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_vf(i)%sf(j, k, l) = real(real(q_cons_vf(i)%sf(j, k, l), wp)/max(real(w_vf(1)%sf(j, k, l), wp), &
                                  & lso_w_floor), stp)
                    end do
                end do
            end do
        end do

    end subroutine s_apply_lso_pp_filter_masked

    !> Refresh ghosts for direction mpi_dir (1=x, 2=y, 3=z): MPI faces via s_mpi_sendrecv_variables_buffers, BC_GHOST_EXTRAP faces
    !! re-extrapolated from the edge cell, BC_PERIODIC faces of an undecomposed direction wrapped locally.
    !! @param q_cons_vf Fields whose ghosts are refreshed
    !! @param mpi_dir Sweep direction
    impure subroutine s_lso_pp_filter_ghost_refresh(q_cons_vf, mpi_dir)

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

#ifdef MFC_MPI
        if (beg_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, -1, nv)
        if (end_bc >= 0) call s_mpi_sendrecv_variables_buffers(q_cons_vf, mpi_dir, 1, nv)
#endif

        select case (mpi_dir)
        case (1)
            do i = 1, nv
                do l = 0, p
                    do k = 0, n
                        do j = 1, buff_size
                            if (beg_bc == BC_GHOST_EXTRAP) then
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(0, k, l)
                            else if (beg_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(-j, k, l) = q_cons_vf(i)%sf(m - j + 1, k, l)
                            end if
                            if (end_bc == BC_GHOST_EXTRAP) then
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(m, k, l)
                            else if (end_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(m + j, k, l) = q_cons_vf(i)%sf(j - 1, k, l)
                            end if
                        end do
                    end do
                end do
            end do
        case (2)
            do i = 1, nv
                do l = 0, p
                    do k = 1, buff_size
                        do j = 0, m
                            if (beg_bc == BC_GHOST_EXTRAP) then
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, 0, l)
                            else if (beg_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(j, -k, l) = q_cons_vf(i)%sf(j, n - k + 1, l)
                            end if
                            if (end_bc == BC_GHOST_EXTRAP) then
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, n, l)
                            else if (end_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(j, n + k, l) = q_cons_vf(i)%sf(j, k - 1, l)
                            end if
                        end do
                    end do
                end do
            end do
        case (3)
            do i = 1, nv
                do l = 1, buff_size
                    do k = 0, n
                        do j = 0, m
                            if (beg_bc == BC_GHOST_EXTRAP) then
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, 0)
                            else if (beg_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(j, k, -l) = q_cons_vf(i)%sf(j, k, p - l + 1)
                            end if
                            if (end_bc == BC_GHOST_EXTRAP) then
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, p)
                            else if (end_bc == BC_PERIODIC) then
                                q_cons_vf(i)%sf(j, k, p + l) = q_cons_vf(i)%sf(j, k, l - 1)
                            end if
                        end do
                    end do
                end do
            end do
        end select

    end subroutine s_lso_pp_filter_ghost_refresh

end module m_lso_pp_filter
