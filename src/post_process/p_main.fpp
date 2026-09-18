!>
!! @file
!! @brief Contains program p_main

!> Post-process raw simulation data into formatted database files (Silo-HDF5 or Binary)
program p_main

    use m_global_parameters
    use m_start_up
    use m_data_output, only: s_switch_output_dirs
    use m_data_input, only: q_cons_vf
    use m_derived_types, only: scalar_field
    use m_lso_pp_filter, only: s_apply_lso_pp_filter, s_compute_lso_pp_stat_fields, s_apply_lso_pp_filter_masked, &
        & s_lso_pp_mask_from_ib

    implicit none

    integer :: t_step, nt_step  !< Iterator for the main time-stepping loop
    !> Generic storage for the name(s) of the flow variable(s) that will be added to the formatted database file(s)
    character(LEN=name_len) :: varname
    real(wp)                :: pres
    real(wp)                :: c
    real(wp)                :: H
    real(wp)                :: start, finish

    call s_initialize_mpi_domain()

    call s_initialize_modules()

    if (cfl_dt) then
        t_step = n_start
        n_save = int(t_stop/t_save) + 1
    else
        ! Setting the time-step iterator to the first time step to be post-processed
        t_step = t_step_start
    end if

    ! Time-Marching Loop
    do
        ! If all time-steps are not ready to be post-processed and one rank is faster than another, the slower rank processing the
        ! last available step might be killed when the faster rank attempts to process the first missing step, before the slower
        ! rank finishes writing the last available step. To avoid this, we force synchronization here.
        call s_mpi_barrier()

        call cpu_time(start)

        ! Primary pass: read (LSO-filtered when lso_filter_wrt=T) and write.
        call s_perform_time_step(t_step)

        ! Post-process filter: filter q_cons_vf in place (pre-filtered coarse data or original data), reconvert
        ! to primitive, and compute the stat products while q_cons_vf still holds the filtered state.
        if (lso_pp_filter) then
            ! Normalization weight: the simulation-written filtered mask, the ib_markers gas mask, or none.
            block
                type(scalar_field) :: w_vf(1:1)
                logical            :: have_w
                have_w = .false.
                allocate (w_vf(1)%sf(lbound(q_cons_vf(1)%sf, 1):ubound(q_cons_vf(1)%sf, 1),lbound(q_cons_vf(1)%sf, &
                          & 2):ubound(q_cons_vf(1)%sf, 2),lbound(q_cons_vf(1)%sf, 3):ubound(q_cons_vf(1)%sf, 3)))
                if (lso_filter_wrt) then
                    call s_read_lso_mask(w_vf, t_step, have_w)
                    if (.not. have_w .and. ib .and. proc_rank == 0) then
                        print '(A)', 'Warning: lso_mask file not found; post_process filter runs without IB normalization.'
                    end if
                else if (ib) then
                    call s_lso_pp_mask_from_ib(w_vf)
                    have_w = .true.
                end if
                ! Original data: the stat products are formed from the unfiltered state (phase-weighted by the gas mask)
                ! and filtered afterwards, so the closures below are F[rho u u] - F[rho u]F[rho u]/F[rho] and not zero.
                if (lso_stat_wrt .and. .not. lso_filter_wrt) then
                    if (have_w) then
                        call s_compute_lso_pp_stat_fields(q_cons_vf, w_vf)
                    else
                        call s_compute_lso_pp_stat_fields(q_cons_vf)
                    end if
                end if
                if (have_w) then
                    call s_apply_lso_pp_filter_masked(q_cons_vf, w_vf)
                else
                    call s_apply_lso_pp_filter(q_cons_vf)
                end if
                deallocate (w_vf(1)%sf)
            end block
            call s_reconvert_filtered_to_primitive()
            if (lso_stat_wrt) then
                if (lso_filter_wrt) then
                    call s_compute_lso_pp_stat_fields(q_cons_vf)  ! widening of pre-filtered data: products of the widened state
                else
                    call s_filter_lso_pp_stat_fields()
                end if
            end if
        end if

        call s_save_data(t_step, varname, pres, c, H)

        ! Secondary pass: when lso_filter_wrt=T, also write the unfiltered data to silo_hdf5/. Grid files are NOT re-opened
        ! (grid_loaded flag skips them), avoiding a known MPI-IO hang from re-opening the same file with
        ! MPI_FILE_OPEN(MPI_COMM_WORLD, fp) a second time. NOTE: When lso_down_sample_factor > 1 the LSO files are on a coarser grid
        ! than the unfiltered data. The secondary pass cannot switch grid dimensions at runtime, so it is skipped. Run post_process
        ! again with lso_filter_wrt=.false. to process the full-resolution unfiltered data.
        if (lso_filter_wrt .and. lso_down_sample_factor <= 1) then
            lso_filter_wrt = .false.
            call s_reload_data(t_step)
            call s_switch_output_dirs(.false.)
            call s_save_data(t_step, varname, pres, c, H)
            call s_switch_output_dirs(.true.)
            lso_filter_wrt = .true.
        end if

        ! Third pass: LSO statistical product fields to silo_hdf5_lso_stat/ (computed above in post_process
        ! filter mode, otherwise read from the simulation-written lso_stat binary).
        if (lso_stat_wrt) then
            if (lso_pp_filter) then
                call s_save_lso_pp_stat_data(t_step)
            else if (lso_filter_wrt) then
                call s_save_lso_stat_data(t_step)
            end if
        end if

        ! Fourth pass: Euler-Lagrange closure fields to silo_hdf5_lso_closure/, at the in-situ width or
        ! (lso_pp_filter with downsampled input) at the widened width.
        if (lso_closure_wrt) then
            if (lso_filter_wrt .and. lso_stat_wrt .and. .not. lso_pp_filter) then
                call s_save_lso_closure_data(t_step)
            else if (lso_filter_wrt .and. lso_stat_wrt .and. lso_pp_filter .and. lso_down_sample_factor > 1) then
                call s_save_lso_pp_closure_data(t_step)
            else if (lso_pp_filter .and. lso_stat_wrt .and. .not. lso_filter_wrt) then
                call s_save_lso_pp_raw_closure_data(t_step)
            else if (proc_rank == 0 .and. t_step == t_step_start) then
                print '(A)', &
                    & 'Warning: lso_closure_wrt requires lso_stat_wrt=T and lso_filter_wrt=T or lso_pp_filter=T ' &
                    & // '(both: lso_down_sample_factor>1); closure output skipped.'
            end if
        end if

        call cpu_time(finish)

        wall_time = abs(finish - start)

        if (cfl_dt) then
            nt_step = t_step - n_start + 1
        else
            nt_step = (t_step - t_step_start)/t_step_save + 1
        end if
        if (nt_step >= 2) then
            wall_time_avg = (wall_time + (nt_step - 2)*wall_time_avg)/(nt_step - 1)
        else
            wall_time_avg = 0._wp
        end if

        if (cfl_dt) then
            if (t_step == n_save - 1) then
                exit
            end if
        else
            ! Adjust time-step iterator to reach final step if needed, else exit
            if ((t_step_stop - t_step) < t_step_save .and. t_step_stop /= t_step) then
                t_step = t_step_stop - t_step_save
            else if (t_step == t_step_stop) then
                exit
            end if
        end if

        if (cfl_dt) then
            t_step = t_step + 1
        else
            ! Incrementing time-step iterator to next time-step to be post-processed
            t_step = t_step + t_step_save
        end if
    end do
    ! END: Time-Marching Loop

    close (11)

    call s_finalize_modules()
end program p_main
