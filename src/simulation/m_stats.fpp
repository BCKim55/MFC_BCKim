!>
!! @file
!! @brief Contains module m_stats

#:include 'macros.fpp'

!> @brief Running time integrals of primitive moments for Reynolds and Favre statistics.
!!
!! Per interior cell, stats_sum holds int dt * f over t >= t_stats_start for
!! f = rho, p, p^2, u_i, rho u_i, u_i u_j, rho u_i u_j (i <= j), in that order; stats_time = int dt.
!! Sums, not means, are stored so a restart continues the integral. File restart_data/stats_<save>.dat:
!! stats_time, then each moment as a global (m_glb+1)(n_glb+1)(p_glb+1) block, x fastest, all in wp.
module m_stats

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy

    implicit none

    private; public :: s_initialize_stats_module, s_accumulate_stats, s_write_stats_file, s_finalize_stats_module

    ! wp, not stp: half-precision sums over 1e5 steps lose the mean
    real(wp), allocatable, dimension(:,:,:,:) :: stats_sum
    $:GPU_DECLARE(create='[stats_sum]')

    real(wp) :: stats_time
    integer  :: n_stats

contains

    !> Allocate the sums and, on a restart, continue from the saved ones
    impure subroutine s_initialize_stats_module(restart_save)

        integer, intent(in) :: restart_save  !< save index restarted from, 0 for a fresh run

        n_stats = 3 + 2*num_vels + num_vels*(num_vels + 1)
        @:ALLOCATE(stats_sum(0:m, 0:n, 0:p, 1:n_stats))
        stats_sum = 0._wp
        stats_time = 0._wp
        if (restart_save > 0) call s_read_stats_file(restart_save)
        $:GPU_UPDATE(device='[stats_sum]')

    end subroutine s_initialize_stats_module

    !> Add dt times each moment of the time-level-n primitives
    impure subroutine s_accumulate_stats(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        real(wp)                                            :: w, rho, pres
        real(wp), dimension(3)                              :: vel
        integer                                             :: i, j, k, l, q, c

        if (mytime < t_stats_start) return
        w = dt
        stats_time = stats_time + w

        $:GPU_PARALLEL_LOOP(collapse=3, private='[i, q, c, rho, pres, vel]', copyin='[w]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    rho = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = eqn_idx%cont%beg, eqn_idx%cont%end
                        rho = rho + q_prim_vf(i)%sf(j, k, l)
                    end do
                    pres = q_prim_vf(eqn_idx%E)%sf(j, k, l)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel(i) = q_prim_vf(eqn_idx%mom%beg + i - 1)%sf(j, k, l)
                    end do

                    stats_sum(j, k, l, 1) = stats_sum(j, k, l, 1) + w*rho
                    stats_sum(j, k, l, 2) = stats_sum(j, k, l, 2) + w*pres
                    stats_sum(j, k, l, 3) = stats_sum(j, k, l, 3) + w*pres*pres
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        stats_sum(j, k, l, 3 + i) = stats_sum(j, k, l, 3 + i) + w*vel(i)
                        stats_sum(j, k, l, 3 + num_vels + i) = stats_sum(j, k, l, 3 + num_vels + i) + w*rho*vel(i)
                    end do
                    c = 3 + 2*num_vels
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        $:GPU_LOOP(parallelism='[seq]')
                        do q = i, num_vels
                            c = c + 1
                            stats_sum(j, k, l, c) = stats_sum(j, k, l, c) + w*vel(i)*vel(q)
                            stats_sum(j, k, l, c + num_vels*(num_vels + 1)/2) = stats_sum(j, k, l, &
                                      & c + num_vels*(num_vels + 1)/2) + w*rho*vel(i)*vel(q)
                        end do
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_accumulate_stats

    !> Write stats_time and the sums to restart_data/stats_<save>.dat
    impure subroutine s_write_stats_file(save_count)

        integer, intent(in) :: save_count

#ifdef MFC_MPI
        integer                              :: ifile, ierr, view, i
        integer, dimension(MPI_STATUS_SIZE)  :: status
        integer(kind=MPI_OFFSET_KIND)        :: disp, wp_bytes, n_glb_cells
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist

        $:GPU_UPDATE(host='[stats_sum]')
        call s_stats_view(view, wp_bytes, n_glb_cells)

        write (file_loc, '(A,I0,A)') '/restart_data/stats_', save_count, '.dat'
        file_loc = trim(case_dir) // trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (file_exist .and. proc_rank == 0) call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
        call s_mpi_barrier()
        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)
        call s_check_mpi_file_open(ierr, file_loc)

        if (proc_rank == 0) call MPI_FILE_WRITE_AT(ifile, 0_MPI_OFFSET_KIND, stats_time, 1, mpi_p, status, ierr)
        do i = 1, n_stats
            disp = wp_bytes*(1 + n_glb_cells*int(i - 1, MPI_OFFSET_KIND))
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)
            call MPI_FILE_WRITE_ALL(ifile, stats_sum(:,:,:,i), (m + 1)*(n + 1)*(p + 1), mpi_p, status, ierr)
        end do

        call MPI_FILE_CLOSE(ifile, ierr)
        call MPI_TYPE_FREE(view, ierr)
#endif

    end subroutine s_write_stats_file

    !> Read stats_time and the sums from restart_data/stats_<save>.dat; start from zero if it is absent
    impure subroutine s_read_stats_file(save_count)

        integer, intent(in) :: save_count

#ifdef MFC_MPI
        integer                              :: ifile, ierr, view, i
        integer, dimension(MPI_STATUS_SIZE)  :: status
        integer(kind=MPI_OFFSET_KIND)        :: disp, wp_bytes, n_glb_cells
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist

        write (file_loc, '(A,I0,A)') '/restart_data/stats_', save_count, '.dat'
        file_loc = trim(case_dir) // trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            if (proc_rank == 0) print '(A)', 'stats: ' // trim(file_loc) // ' not found, statistics start from zero'
            return
        end if

        call s_stats_view(view, wp_bytes, n_glb_cells)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, MPI_MODE_RDONLY, mpi_info_int, ifile, ierr)
        call s_check_mpi_file_open(ierr, file_loc)

        call MPI_FILE_READ_AT_ALL(ifile, 0_MPI_OFFSET_KIND, stats_time, 1, mpi_p, status, ierr)
        do i = 1, n_stats
            disp = wp_bytes*(1 + n_glb_cells*int(i - 1, MPI_OFFSET_KIND))
            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, view, 'native', mpi_info_int, ierr)
            call MPI_FILE_READ_ALL(ifile, stats_sum(:,:,:,i), (m + 1)*(n + 1)*(p + 1), mpi_p, status, ierr)
        end do

        call MPI_FILE_CLOSE(ifile, ierr)
        call MPI_TYPE_FREE(view, ierr)
#endif

    end subroutine s_read_stats_file

#ifdef MFC_MPI
    !> Subarray type placing this rank's interior in the global array
    impure subroutine s_stats_view(view, wp_bytes, n_glb_cells)

        integer, intent(out)                       :: view
        integer(kind=MPI_OFFSET_KIND), intent(out) :: wp_bytes, n_glb_cells
        integer, dimension(num_dims)               :: sizes_glb, sizes_loc
        integer                                    :: ierr

        sizes_glb(1) = m_glb + 1; sizes_loc(1) = m + 1
        if (n > 0) then
            sizes_glb(2) = n_glb + 1; sizes_loc(2) = n + 1
            if (p > 0) then
                sizes_glb(3) = p_glb + 1; sizes_loc(3) = p + 1
            end if
        end if
        call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_idx, MPI_ORDER_FORTRAN, mpi_p, view, ierr)
        call MPI_TYPE_COMMIT(view, ierr)

        wp_bytes = int(storage_size(0._wp)/8, MPI_OFFSET_KIND)
        n_glb_cells = int(m_glb + 1, MPI_OFFSET_KIND)*int(n_glb + 1, MPI_OFFSET_KIND)*int(p_glb + 1, MPI_OFFSET_KIND)

    end subroutine s_stats_view
#endif

    impure subroutine s_finalize_stats_module

        @:DEALLOCATE(stats_sum)

    end subroutine s_finalize_stats_module

end module m_stats
