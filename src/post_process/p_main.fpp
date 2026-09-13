!>
!! @file
!! @brief Contains program p_main

!> Post-process raw simulation data into formatted database files (Silo-HDF5 or Binary)
program p_main

    use m_global_parameters
    use m_start_up
    use m_derived_types, only: scalar_field
    use m_data_input, only: q_cons_vf
    use m_lso_pp_filter, only: s_apply_lso_pp_filter, s_apply_lso_pp_filter_masked, s_lso_pp_mask_from_ib

    implicit none

    integer :: t_step, nt_step  !< Iterator for the main time-stepping loop
    !> Generic storage for the name(s) of the flow variable(s) that will be added to the formatted database file(s)
    character(LEN=name_len) :: varname
    real(wp)                :: pres
    real(wp)                :: c
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

        call s_perform_time_step(t_step)

        ! LSO filter: filter q_cons_vf in place (mask-normalized when immersed boundaries
        ! are present), then rebuild the primitive state from the filtered fields.
        if (lso_pp_filter) then
            block
                type(scalar_field) :: w_vf(1:1)

                if (ib) then
                    allocate (w_vf(1)%sf(lbound(q_cons_vf(1)%sf, 1):ubound(q_cons_vf(1)%sf, 1),lbound(q_cons_vf(1)%sf, &
                              & 2):ubound(q_cons_vf(1)%sf, 2),lbound(q_cons_vf(1)%sf, 3):ubound(q_cons_vf(1)%sf, 3)))
                    call s_lso_pp_mask_from_ib(w_vf)
                    call s_apply_lso_pp_filter_masked(q_cons_vf, w_vf)
                    deallocate (w_vf(1)%sf)
                else
                    call s_apply_lso_pp_filter(q_cons_vf)
                end if
            end block
            call s_reconvert_filtered_to_primitive()
        end if

        call s_save_data(t_step, varname, pres, c)

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
