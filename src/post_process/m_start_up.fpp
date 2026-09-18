#:include 'macros.fpp'

!>
!! @file
!! @brief  Contains module m_start_up

!> @brief Reads and validates user inputs, allocates variables, and configures MPI decomposition and I/O for post-processing

module m_start_up

    use, intrinsic :: iso_c_binding

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_boundary_common
    use m_variables_conversion
    use m_data_input
    use m_data_output
    use m_derived_variables
    use m_helper
    use m_compile_specific
    use m_checker_common
    use m_checker
    use m_thermochem, only: num_species, species_names
    use m_finite_differences
    use m_chemistry
    use m_lso_pp_filter

#ifdef MFC_MPI
    use mpi
#endif

    implicit none

    include 'fftw3.f03'

    type(c_ptr)                             :: fwd_plan_x, fwd_plan_y, fwd_plan_z
    complex(c_double_complex), allocatable  :: data_in(:), data_out(:)
    complex(c_double_complex), allocatable  :: data_cmplx(:,:,:), data_cmplx_y(:,:,:), data_cmplx_z(:,:,:)
    real(wp), allocatable, dimension(:,:,:) :: En_real
    real(wp), allocatable, dimension(:)     :: En
    integer                                 :: num_procs_x, num_procs_y, num_procs_z
    integer                                 :: Nx, Ny, Nz, Nxloc, Nyloc, Nyloc2, Nzloc, Nf
    integer                                 :: ierr
    integer                                 :: MPI_COMM_CART, MPI_COMM_CART12, MPI_COMM_CART13
    integer, dimension(3)                   :: cart3d_coords
    integer, dimension(2)                   :: cart2d12_coords, cart2d13_coords
    integer                                 :: proc_rank12, proc_rank13

contains

    !> Reads the configuration file post_process.inp, in order to populate parameters in module m_global_parameters.f90 with the
    !! user provided inputs
    impure subroutine s_read_input_file

        character(LEN=name_len) :: file_loc
        logical                 :: file_check
        integer                 :: iostatus
        character(len=1000)     :: line

        namelist /user_inputs/ case_dir, m, n, p, t_step_start, t_step_stop, t_step_save, model_eqns, num_fluids, mpp_lim, &
            & weno_order, bc_x, bc_y, bc_z, fluid_pp, bub_pp, format, precision, output_partial_domain, x_output, y_output, &
            & z_output, hypoelasticity, G, mhd, chem_wrt_Y, chem_wrt_T, avg_state, alpha_rho_wrt, rho_wrt, mom_wrt, vel_wrt, &
            & E_wrt, fft_wrt, pres_wrt, alpha_wrt, gamma_wrt, heat_ratio_wrt, pi_inf_wrt, pres_inf_wrt, cons_vars_wrt, &
            & prim_vars_wrt, c_wrt, omega_wrt, qm_wrt, liutex_wrt, schlieren_wrt, schlieren_alpha, fd_order, mixture_err, &
            & alt_soundspeed, flux_lim, flux_wrt, cyl_coord, parallel_io, rhoref, pref, bubbles_euler, qbmm, sigR, R0ref, nb, &
            & polytropic, thermal, Ca, Web, Re_inv, polydisperse, poly_sigma, file_per_process, relax, relax_model, cf_wrt, &
            & sigma, adv_n, ib, num_ibs, cfl_adap_dt, cfl_const_dt, t_save, t_stop, n_start, cfl_target, surface_tension, &
            & bubbles_lagrange, sim_data, hyperelasticity, Bx0, relativity, cont_damage, hyper_cleaning, num_bc_patches, igr, &
            & igr_order, down_sample, recon_type, muscl_order, lag_header, lag_txt_wrt, lag_db_wrt, lag_id_wrt, lag_pos_wrt, &
            & lag_pos_prev_wrt, lag_vel_wrt, lag_rad_wrt, lag_rvel_wrt, lag_r0_wrt, lag_rmax_wrt, lag_rmin_wrt, lag_dphidt_wrt, &
            & lag_pres_wrt, lag_mv_wrt, lag_mg_wrt, lag_betaT_wrt, lag_betaC_wrt, alpha_rho_e_wrt, ib_state_wrt, lso_filter_wrt, &
            & lso_down_sample_factor, lso_stat_wrt, lso_pp_filter, lso_closure_wrt, lso_pp_n_passes_x, lso_pp_n_passes_y, &
            & lso_pp_n_passes_z, lso_pp_a_x, lso_pp_a_y, lso_pp_a_z, lso_R_gas, lso_mu, lso_conductivity

        file_loc = 'post_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        if (file_check) then
            open (1, FILE=trim(file_loc), form='formatted', STATUS='old', ACTION='read')
            read (1, NML=user_inputs, iostat=iostatus)

            if (iostatus /= 0) then
                backspace (1)
                read (1, fmt='(A)') line
                print *, 'Invalid line in namelist: ' // trim(line)
                call s_mpi_abort('Invalid line in post_process.inp. It is ' // 'likely due to a datatype mismatch. Exiting.')
            end if

            close (1)

            call s_update_cell_bounds(cells_bounds, m, n, p)

            if (down_sample) then
                m = int((m + 1)/3) - 1
                n = int((n + 1)/3) - 1
                p = int((p + 1)/3) - 1
            end if

            ! LSO-filtered output saved on a coarser grid: shrink m/n/p to match the file.
            if (lso_filter_wrt .and. lso_down_sample_factor > 1) then
                m = int((m + 1)/lso_down_sample_factor) - 1
                if (n > 0) n = int((n + 1)/lso_down_sample_factor) - 1
                if (p > 0) p = int((p + 1)/lso_down_sample_factor) - 1
            end if

            ! LSO stat index layout (same as simulation). num_dims is not set yet, so derive a local dim count from m/n/p.
            if (lso_filter_wrt .and. lso_stat_wrt) then
                block
                    integer :: loc_num_dims
                    loc_num_dims = 1 + min(1, n) + min(1, p)
                    lso_stat_phi_p_beg = 1; lso_stat_phi_p_end = 1
                    lso_stat_rho_beg = lso_stat_phi_p_end + 1
                    lso_stat_rho_end = lso_stat_rho_beg
                    lso_stat_rhoke_beg = lso_stat_rho_end + 1
                    lso_stat_rhoke_end = lso_stat_rhoke_beg
                    lso_stat_up_beg = lso_stat_rhoke_end + 1
                    lso_stat_up_end = lso_stat_up_beg + loc_num_dims - 1
                    lso_stat_rhou_beg = lso_stat_up_end + 1
                    lso_stat_rhou_end = lso_stat_rhou_beg + loc_num_dims - 1
                    lso_stat_rhouu_beg = lso_stat_rhou_end + 1
                    if (loc_num_dims == 1) then
                        lso_stat_rhouu_end = lso_stat_rhouu_beg
                    else if (loc_num_dims == 2) then
                        lso_stat_rhouu_end = lso_stat_rhouu_beg + 2
                    else
                        lso_stat_rhouu_end = lso_stat_rhouu_beg + 5
                    end if
                    lso_stat_rhouke_beg = lso_stat_rhouu_end + 1
                    lso_stat_rhouke_end = lso_stat_rhouke_beg + loc_num_dims - 1
                    lso_stat_rhouT_beg = lso_stat_rhouke_end + 1
                    lso_stat_rhouT_end = lso_stat_rhouT_beg + loc_num_dims - 1
                    lso_stat_tau_beg = lso_stat_rhouT_end + 1
                    lso_stat_tau_end = lso_stat_tau_beg + (lso_stat_rhouu_end - lso_stat_rhouu_beg)
                    lso_stat_q_beg = lso_stat_tau_end + 1
                    lso_stat_q_end = lso_stat_q_beg + loc_num_dims - 1
                    lso_stat_rhotau_u_beg = lso_stat_q_end + 1
                    lso_stat_rhotau_u_end = lso_stat_rhotau_u_beg + loc_num_dims - 1
                    n_lso_stat = lso_stat_rhotau_u_end
                end block
            end if

            ! Same layout for the post_process filter path.
            if (lso_pp_filter .and. lso_stat_wrt) then
                block
                    integer :: loc_num_dims
                    loc_num_dims = 1 + min(1, n) + min(1, p)
                    lso_stat_phi_p_beg = 1; lso_stat_phi_p_end = 1
                    lso_stat_rho_beg = lso_stat_phi_p_end + 1
                    lso_stat_rho_end = lso_stat_rho_beg
                    lso_stat_rhoke_beg = lso_stat_rho_end + 1
                    lso_stat_rhoke_end = lso_stat_rhoke_beg
                    lso_stat_up_beg = lso_stat_rhoke_end + 1
                    lso_stat_up_end = lso_stat_up_beg + loc_num_dims - 1
                    lso_stat_rhou_beg = lso_stat_up_end + 1
                    lso_stat_rhou_end = lso_stat_rhou_beg + loc_num_dims - 1
                    lso_stat_rhouu_beg = lso_stat_rhou_end + 1
                    if (loc_num_dims == 1) then
                        lso_stat_rhouu_end = lso_stat_rhouu_beg
                    else if (loc_num_dims == 2) then
                        lso_stat_rhouu_end = lso_stat_rhouu_beg + 2
                    else
                        lso_stat_rhouu_end = lso_stat_rhouu_beg + 5
                    end if
                    lso_stat_rhouke_beg = lso_stat_rhouu_end + 1
                    lso_stat_rhouke_end = lso_stat_rhouke_beg + loc_num_dims - 1
                    lso_stat_rhouT_beg = lso_stat_rhouke_end + 1
                    lso_stat_rhouT_end = lso_stat_rhouT_beg + loc_num_dims - 1
                    lso_stat_tau_beg = lso_stat_rhouT_end + 1
                    lso_stat_tau_end = lso_stat_tau_beg + (lso_stat_rhouu_end - lso_stat_rhouu_beg)
                    lso_stat_q_beg = lso_stat_tau_end + 1
                    lso_stat_q_end = lso_stat_q_beg + loc_num_dims - 1
                    lso_stat_rhotau_u_beg = lso_stat_q_end + 1
                    lso_stat_rhotau_u_end = lso_stat_rhotau_u_beg + loc_num_dims - 1
                    n_lso_stat = lso_stat_rhotau_u_end
                end block
            end if

            m_glb = m
            n_glb = n
            p_glb = p

            nGlobal = int(m_glb + 1, kind=8)*int(n_glb + 1, kind=8)*int(p_glb + 1, kind=8)

            if (cfl_adap_dt .or. cfl_const_dt) cfl_dt = .true.

            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -17) .or. num_bc_patches > 0) then
                bc_io = .true.
            end if
        else
            call s_mpi_abort('File post_process.inp is missing. Exiting.')
        end if

    end subroutine s_read_input_file

    !> Checking that the user inputs make sense, i.e. that the individual choices are compatible with the code's options and that
    !! the combination of these choices results into a valid configuration for the post-process
    impure subroutine s_check_input_file

        character(LEN=len_trim(case_dir)) :: file_loc
        logical                           :: dir_check

        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir) // '/.'

        call my_inquire(file_loc, dir_check)

        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Unsupported choice for the value of ' // 'case_dir. Exiting.')
        end if

        call s_check_inputs_common()
        call s_check_inputs()

    end subroutine s_check_input_file

    !> Load grid and conservative data for a time step, fill ghost-cell buffers, and convert to primitive variables.
    impure subroutine s_perform_time_step(t_step)

        integer, intent(inout) :: t_step
        integer                :: eta_hh, eta_mm, eta_ss
        real(wp)               :: eta_sec

        if (proc_rank == 0) then
            if (cfl_dt) then
                eta_sec = wall_time_avg*real(n_save - 1 - t_step, wp)
                eta_hh = int(eta_sec)/3600
                eta_mm = mod(int(eta_sec), 3600)/60
                eta_ss = mod(int(eta_sec), 60)
                print '(" [", I3, "%]  Saving ", I8, " of ", I0, " Time Avg = ", ES16.6,  " Time/step = ", ES12.6, " ETA (HH:MM:SS)  = ", I0, ":", I2.2, ":", I2.2)', &
                    & int(ceiling(100._wp*(real(t_step - n_start)/(n_save)))), t_step, n_save, wall_time_avg, wall_time, eta_hh, &
                    & eta_mm, eta_ss
            else
                eta_sec = wall_time_avg*real((t_step_stop - t_step)/t_step_save, wp)
                eta_hh = int(eta_sec)/3600
                eta_mm = mod(int(eta_sec), 3600)/60
                eta_ss = mod(int(eta_sec), 60)
                print '(" [", I3, "%]  Saving ", I8, " of ", I0, " @ t_step = ", I8, " Time Avg = ", ES16.6,  " Time/step = ", ES12.6, " ETA (HH:MM:SS) = ", I0, ":", I2.2, ":", I2.2)', &
                    & int(ceiling(100._wp*(real(t_step - t_step_start)/(t_step_stop - t_step_start + 1)))), &
                    & (t_step - t_step_start)/t_step_save + 1, (t_step_stop - t_step_start)/t_step_save + 1, t_step, &
                    & wall_time_avg, wall_time, eta_hh, eta_mm, eta_ss
            end if
        end if

        call s_read_data_files(t_step)

        if (chemistry) call s_compute_q_T_sf(q_T_sf, q_cons_vf, idwbuff)

        if (buff_size > 0) then
            call s_populate_grid_variables_buffers()
            call s_populate_variables_buffers(bc_type, q_cons_vf, q_T_sf=q_T_sf)
        end if

        call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_vf, idwbuff)

    end subroutine s_perform_time_step

    !> Reload conservative variable data for a time step and reconvert to primitive variables, WITHOUT printing the progress bar.
    !! Used for the LSO two-pass write: after the normal filtered pass the caller toggles lso_filter_wrt=.false., calls
    !! s_reload_data to read the unfiltered conservative data, then writes a second Silo database to silo_hdf5/. Grid file reads are
    !! skipped automatically (grid_loaded flag in m_data_input).
    impure subroutine s_reload_data(t_step)

        integer, intent(in) :: t_step

        call s_read_data_files(t_step)

        if (chemistry) call s_compute_q_T_sf(q_T_sf, q_cons_vf, idwbuff)

        if (buff_size > 0) then
            call s_populate_grid_variables_buffers()
            call s_populate_variables_buffers(bc_type, q_cons_vf, q_T_sf=q_T_sf)
        end if

        call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_vf, idwbuff)

    end subroutine s_reload_data

    !> Reconvert q_cons_vf to q_prim_vf after the post_process LSO filter modified it in place (same sequence as s_reload_data).
    impure subroutine s_reconvert_filtered_to_primitive()

        if (chemistry) call s_compute_q_T_sf(q_T_sf, q_cons_vf, idwbuff)

        if (buff_size > 0) then
            call s_populate_variables_buffers(bc_type, q_cons_vf, q_T_sf=q_T_sf)
        end if

        call s_convert_conservative_to_primitive_variables(q_cons_vf, q_T_sf, q_prim_vf, idwbuff)

    end subroutine s_reconvert_filtered_to_primitive

    !> Derive requested flow quantities from primitive variables and write them to the formatted database files.
    impure subroutine s_save_data(t_step, varname, pres, c, H)

        integer, intent(inout)                 :: t_step
        character(LEN=name_len), intent(inout) :: varname
        real(wp), intent(inout)                :: pres, c, H

        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end, &
             & -offset_z%beg:p + offset_z%end) :: liutex_mag
        real(wp), dimension(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end, &
             & 3) :: liutex_axis
        integer       :: i, j, k, l, kx, ky, kz, kf, j_glb, k_glb, l_glb
        character(50) :: filename
        logical       :: file_exists
        integer       :: x_beg, x_end, y_beg, y_end, z_beg, z_end

        if (output_partial_domain) then
            call s_define_output_region
            x_beg = -offset_x%beg + x_output_idx%beg
            x_end = offset_x%end + x_output_idx%end
            y_beg = -offset_y%beg + y_output_idx%beg
            y_end = offset_y%end + y_output_idx%end
            z_beg = -offset_z%beg + z_output_idx%beg
            z_end = offset_z%end + z_output_idx%end
        else
            x_beg = -offset_x%beg
            x_end = offset_x%end + m
            y_beg = -offset_y%beg
            y_end = offset_y%end + n
            z_beg = -offset_z%beg
            z_end = offset_z%end + p
        end if

        call s_open_formatted_database_file(t_step)

        if (sim_data .and. proc_rank == 0) then
            call s_open_intf_data_file()
            call s_open_energy_data_file()
        end if

        if (sim_data) then
            call s_write_intf_data_file(q_prim_vf)
            call s_write_energy_data_file(q_prim_vf, q_cons_vf)
        end if

        call s_write_grid_to_formatted_database_file(t_step)

        if (omega_wrt(2) .or. omega_wrt(3) .or. qm_wrt .or. liutex_wrt .or. schlieren_wrt) then
            call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, fd_number, fd_order, offset_x)
        end if

        if (omega_wrt(1) .or. omega_wrt(3) .or. qm_wrt .or. liutex_wrt .or. (n > 0 .and. schlieren_wrt)) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, fd_number, fd_order, offset_y)
        end if

        if (omega_wrt(1) .or. omega_wrt(2) .or. qm_wrt .or. liutex_wrt .or. (p > 0 .and. schlieren_wrt)) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z, buff_size, fd_number, fd_order, offset_z)
        end if

        if ((model_eqns == 2) .or. (model_eqns == 3) .or. (model_eqns == 4)) then
            do i = 1, num_fluids
                if (alpha_rho_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    q_sf(:,:,:) = q_cons_vf(i)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    if (model_eqns /= 4) then
                        write (varname, '(A,I0)') 'alpha_rho', i
                    else
                        write (varname, '(A,I0)') 'rho', i
                    end if
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '
                end if
            end do
        end if

        if ((rho_wrt .or. (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) .and. (.not. relativity)) then
            q_sf(:,:,:) = rho_sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'rho'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (relativity .and. (rho_wrt .or. prim_vars_wrt)) then
            q_sf(:,:,:) = q_prim_vf(1)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'rho'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (relativity .and. (rho_wrt .or. cons_vars_wrt)) then
            ! For relativistic flow, conservative and primitive densities are different Hard-coded single-component for now
            q_sf(:,:,:) = q_cons_vf(1)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'D'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        do i = 1, eqn_idx%E - eqn_idx%mom%beg
            if (mom_wrt(i) .or. cons_vars_wrt) then
                q_sf(:,:,:) = q_cons_vf(i + eqn_idx%cont%end)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A,I0)') 'mom', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        do i = 1, eqn_idx%E - eqn_idx%mom%beg
            if (vel_wrt(i) .or. prim_vars_wrt) then
                q_sf(:,:,:) = q_prim_vf(i + eqn_idx%cont%end)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A,I0)') 'vel', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        if (chemistry) then
            do i = 1, num_species
                if (chem_wrt_Y(i) .or. prim_vars_wrt) then
                    q_sf(:,:,:) = q_prim_vf(eqn_idx%species%beg + i - 1)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,A)') 'Y_', trim(species_names(i))
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '
                end if
            end do

            if (chem_wrt_T) then
                q_sf(:,:,:) = q_T_sf%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A)') 'T'
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end if

        do i = 1, eqn_idx%E - eqn_idx%mom%beg
            if (flux_wrt(i)) then
                call s_derive_flux_limiter(i, q_prim_vf, q_sf)

                write (varname, '(A,I0)') 'flux', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        if (E_wrt .or. cons_vars_wrt) then
            q_sf(:,:,:) = q_cons_vf(eqn_idx%E)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'E'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (model_eqns == 3) then
            do i = 1, num_fluids
                if (alpha_rho_e_wrt(i) .or. cons_vars_wrt) then
                    q_sf = q_cons_vf(i + eqn_idx%int_en%beg - 1)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,I0)') 'alpha_rho_e', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '
                end if
            end do
        end if

        if (fft_wrt) then
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        data_cmplx(j + 1, k + 1, l + 1) = cmplx(q_cons_vf(eqn_idx%mom%beg)%sf(j, k, l)/q_cons_vf(1)%sf(j, k, l), &
                                   & 0._wp)
                    end do
                end do
            end do

            call s_mpi_FFT_fwd()

            En_real = 0.5_wp*abs(data_cmplx_z)**2._wp/(1._wp*Nx*Ny*Nz)**2._wp

            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        data_cmplx(j + 1, k + 1, l + 1) = cmplx(q_cons_vf(eqn_idx%mom%beg + 1)%sf(j, k, l)/q_cons_vf(1)%sf(j, k, &
                                   & l), 0._wp)
                    end do
                end do
            end do

            call s_mpi_FFT_fwd()

            En_real = En_real + 0.5_wp*abs(data_cmplx_z)**2._wp/(1._wp*Nx*Ny*Nz)**2._wp

            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        data_cmplx(j + 1, k + 1, l + 1) = cmplx(q_cons_vf(eqn_idx%mom%beg + 2)%sf(j, k, l)/q_cons_vf(1)%sf(j, k, &
                                   & l), 0._wp)
                    end do
                end do
            end do

            call s_mpi_FFT_fwd()

            En_real = En_real + 0.5_wp*abs(data_cmplx_z)**2._wp/(1._wp*Nx*Ny*Nz)**2._wp

            do kf = 1, Nf
                En(kf) = 0._wp
            end do

            do l = 1, Nz
                do k = 1, Nyloc2
                    do j = 1, Nxloc
                        j_glb = j + cart3d_coords(2)*Nxloc
                        k_glb = k + cart3d_coords(3)*Nyloc2
                        l_glb = l

                        if (j_glb >= (m_glb + 1)/2) then
                            kx = (j_glb - 1) - (m_glb + 1)
                        else
                            kx = j_glb - 1
                        end if

                        if (k_glb >= (n_glb + 1)/2) then
                            ky = (k_glb - 1) - (n_glb + 1)
                        else
                            ky = k_glb - 1
                        end if

                        if (l_glb >= (p_glb + 1)/2) then
                            kz = (l_glb - 1) - (p_glb + 1)
                        else
                            kz = l_glb - 1
                        end if

                        kf = nint(sqrt(kx**2._wp + ky**2._wp + kz**2._wp)) + 1

                        En(kf) = En(kf) + En_real(j, k, l)
                    end do
                end do
            end do

#ifdef MFC_MPI
            call MPI_ALLREDUCE(MPI_IN_PLACE, En, Nf, mpi_p, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

            if (proc_rank == 0) then
                call s_create_directory('En_FFT_DATA')
                write (filename, '(a,i0,a)') 'En_FFT_DATA/En_tot', t_step, '.dat'
                inquire (FILE=filename, EXIST=file_exists)
                if (file_exists) then
                    call s_delete_file(trim(filename))
                end if
            end if

            do kf = 1, Nf
                if (proc_rank == 0) then
                    write (filename, '(a,i0,a)') 'En_FFT_DATA/En_tot', t_step, '.dat'
                    inquire (FILE=filename, EXIST=file_exists)
                    if (file_exists) then
                        open (1, file=filename, position='append', status='old')
                        write (1, *) En(kf), t_step
                        close (1)
                    else
                        open (1, file=filename, status='new')
                        write (1, *) En(kf), t_step
                        close (1)
                    end if
                end if
            end do
        end if

        if (mhd .and. prim_vars_wrt) then
            do i = eqn_idx%B%beg, eqn_idx%B%end
                q_sf(:,:,:) = q_prim_vf(i)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)

                ! 1D: output By, Bz
                if (n == 0) then
                    if (i == eqn_idx%B%beg) then
                        write (varname, '(A)') 'By'
                    else
                        write (varname, '(A)') 'Bz'
                    end if
                    ! 2D/3D: output Bx, By, Bz
                else
                    if (i == eqn_idx%B%beg) then
                        write (varname, '(A)') 'Bx'
                    else if (i == eqn_idx%B%beg + 1) then
                        write (varname, '(A)') 'By'
                    else
                        write (varname, '(A)') 'Bz'
                    end if
                end if

                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end do
        end if

        if (elasticity) then
            do i = 1, eqn_idx%stress%end - eqn_idx%stress%beg + 1
                if (prim_vars_wrt) then
                    q_sf(:,:,:) = q_prim_vf(i - 1 + eqn_idx%stress%beg)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,I0)') 'tau', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)
                end if
                varname(:) = ' '
            end do
        end if

        if (hyperelasticity) then
            do i = 1, eqn_idx%xi%end - eqn_idx%xi%beg + 1
                if (prim_vars_wrt) then
                    q_sf(:,:,:) = q_prim_vf(i - 1 + eqn_idx%xi%beg)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,I0)') 'xi', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)
                end if
                varname(:) = ' '
            end do
        end if

        if (cont_damage) then
            q_sf(:,:,:) = q_cons_vf(eqn_idx%damage)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'damage_state'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (hyper_cleaning) then
            q_sf = q_cons_vf(eqn_idx%psi)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'psi'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (pres_wrt .or. prim_vars_wrt) then
            q_sf(:,:,:) = q_prim_vf(eqn_idx%E)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'pres'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (((model_eqns == 2) .and. (bubbles_euler .neqv. .true.)) .or. (model_eqns == 3)) then
            do i = 1, num_fluids - 1
                if (alpha_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    q_sf(:,:,:) = q_cons_vf(i + eqn_idx%E)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,I0)') 'alpha', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '
                end if
            end do

            if (alpha_wrt(num_fluids) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                if (igr) then
                    do k = z_beg, z_end
                        do j = y_beg, y_end
                            do i = x_beg, x_end
                                q_sf(i, j, k) = 1._wp
                                do l = 1, num_fluids - 1
                                    q_sf(i, j, k) = q_sf(i, j, k) - q_cons_vf(eqn_idx%E + l)%sf(i, j, k)
                                end do
                            end do
                        end do
                    end do
                else
                    q_sf(:,:,:) = q_cons_vf(eqn_idx%adv%end)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                end if
                write (varname, '(A,I0)') 'alpha', num_fluids
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end if

        if (gamma_wrt .or. (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) then
            q_sf(:,:,:) = gamma_sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'gamma'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (heat_ratio_wrt) then
            call s_derive_specific_heat_ratio(q_sf)

            write (varname, '(A)') 'heat_ratio'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (pi_inf_wrt .or. (model_eqns == 1 .and. (cons_vars_wrt .or. prim_vars_wrt))) then
            q_sf(:,:,:) = pi_inf_sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'pi_inf'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (pres_inf_wrt) then
            call s_derive_liquid_stiffness(q_sf)

            write (varname, '(A)') 'pres_inf'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (c_wrt) then
            do k = -offset_z%beg, p + offset_z%end
                do j = -offset_y%beg, n + offset_y%end
                    do i = -offset_x%beg, m + offset_x%end
                        do l = 1, eqn_idx%adv%end - eqn_idx%E
                            adv(l) = q_prim_vf(eqn_idx%E + l)%sf(i, j, k)
                        end do

                        pres = q_prim_vf(eqn_idx%E)%sf(i, j, k)

                        H = ((gamma_sf(i, j, k) + 1._wp)*pres + pi_inf_sf(i, j, k) + qv_sf(i, j, k))/rho_sf(i, j, k)

                        call s_compute_speed_of_sound(pres, rho_sf(i, j, k), gamma_sf(i, j, k), pi_inf_sf(i, j, k), H, adv, &
                                                      & 0._wp, 0._wp, c, qv_sf(i, j, k))

                        q_sf(i, j, k) = c
                    end do
                end do
            end do

            write (varname, '(A)') 'c'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        do i = 1, 3
            if (omega_wrt(i)) then
                call s_derive_vorticity_component(i, q_prim_vf, q_sf)

                write (varname, '(A,I0)') 'omega', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        if (ib) then
            q_sf(:,:,:) = real(ib_markers%sf(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end, &
                 & -offset_z%beg:p + offset_z%end))
            varname = 'ib_markers'
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if

        if (p > 0 .and. qm_wrt) then
            call s_derive_qm(q_prim_vf, q_sf)

            write (varname, '(A)') 'qm'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (liutex_wrt) then
            call s_derive_liutex(q_prim_vf, liutex_mag, liutex_axis)

            q_sf = liutex_mag

            write (varname, '(A)') 'liutex_mag'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '

            do i = 1, 3
                q_sf = liutex_axis(:,:,:,i)

                write (varname, '(A,I0)') 'liutex_axis', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end do
        end if

        if (schlieren_wrt) then
            call s_derive_numerical_schlieren_function(q_cons_vf, q_sf)

            write (varname, '(A)') 'schlieren'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (cf_wrt) then
            q_sf(:,:,:) = q_cons_vf(eqn_idx%c)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A,I0)') 'color_function'
            call s_write_variable_to_formatted_database_file(varname, t_step)
            varname(:) = ' '
        end if

        if (bubbles_euler) then
            do i = eqn_idx%adv%beg, eqn_idx%adv%end
                q_sf(:,:,:) = q_cons_vf(i)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A,I0)') 'alpha', i - eqn_idx%E
                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end do
        end if

        if (bubbles_euler) then
            ! nR
            do i = 1, nb
                q_sf(:,:,:) = q_cons_vf(qbmm_idx%rs(i))%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A,I3.3)') 'nR', i
                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end do

            ! nRdot
            do i = 1, nb
                q_sf(:,:,:) = q_cons_vf(qbmm_idx%vs(i))%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A,I3.3)') 'nV', i
                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end do
            if ((polytropic .neqv. .true.) .and. (.not. qbmm)) then
                ! nP
                do i = 1, nb
                    q_sf(:,:,:) = q_cons_vf(qbmm_idx%ps(i))%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,I3.3)') 'nP', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)
                    varname(:) = ' '
                end do

                ! nM
                do i = 1, nb
                    q_sf(:,:,:) = q_cons_vf(qbmm_idx%ms(i))%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,I3.3)') 'nM', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)
                    varname(:) = ' '
                end do
            end if

            ! number density
            if (adv_n) then
                q_sf(:,:,:) = q_cons_vf(eqn_idx%n)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A)') 'n'
                call s_write_variable_to_formatted_database_file(varname, t_step)
                varname(:) = ' '
            end if
        end if

        if (bubbles_lagrange) then
            ! Void fraction field
            q_sf(:,:,:) = 1._wp - q_cons_vf(beta_idx)%sf(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end, &
                 & -offset_z%beg:p + offset_z%end)
            write (varname, '(A)') 'voidFraction'
            call s_write_variable_to_formatted_database_file(varname, t_step)
            varname(:) = ' '

            if (lag_txt_wrt) call s_write_lag_bubbles_results_to_text(t_step)  ! text output
            if (lag_db_wrt) call s_write_lag_bubbles_to_formatted_database_file(t_step)  ! silo file output
        end if

        if (ib_state_wrt) call s_write_ib_bodies_to_formatted_database_file(t_step)

        if (sim_data .and. proc_rank == 0) then
            call s_close_intf_data_file()
            call s_close_energy_data_file()
        end if

        call s_close_formatted_database_file()

    end subroutine s_save_data

    !> Read lso_stat_<t_step>.dat (n_stat MPI-IO subarrays back-to-back) into q_stat_vf. Sets found = .false. when the file is
    !! missing.
    impure subroutine s_read_lso_stat_file(q_stat_vf, n_stat, t_step, found, fname)

        type(scalar_field), intent(inout)      :: q_stat_vf(:)
        integer, intent(in)                    :: n_stat, t_step
        logical, intent(out)                   :: found
        character(LEN=*), intent(in), optional :: fname  !< file basename prefix (default 'lso_stat_')

#ifdef MFC_MPI
        integer                              :: ifile, ierr, data_size, i
        integer, dimension(MPI_STATUS_SIZE)  :: status
        integer(kind=MPI_OFFSET_KIND)        :: disp
        integer(kind=MPI_OFFSET_KIND)        :: m_MOK, n_MOK, p_MOK
        integer(kind=MPI_OFFSET_KIND)        :: WP_MOK, var_MOK, MOK
        integer, dimension(num_dims)         :: sizes_glb, sizes_loc, start_stat
        integer                              :: mpi_view
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist

        found = .true.

        sizes_glb(1) = m_glb + 1; sizes_loc(1) = m + 1; start_stat(1) = start_idx(1)
        if (num_dims >= 2) then
            sizes_glb(2) = n_glb + 1; sizes_loc(2) = n + 1; start_stat(2) = start_idx(2)
        end if
        if (num_dims == 3) then
            sizes_glb(3) = p_glb + 1; sizes_loc(3) = p + 1; start_stat(3) = start_idx(3)
        end if
        data_size = (m + 1)*(n + 1)*(p + 1)
        m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
        n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
        p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
        WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
        MOK = int(1._wp, MPI_OFFSET_KIND)

        if (present(fname)) then
            write (file_loc, '(A,I0,A)') trim(fname), t_step, '.dat'
        else
            write (file_loc, '(A,I0,A)') 'lso_stat_', t_step, '.dat'
        end if
        file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
        inquire (FILE=trim(file_loc), EXIST=file_exist)
        if (.not. file_exist) then
            found = .false.
            return
        end if
        ! MPI_COMM_SELF so ranks open independently (collective open would deadlock if some ranks skip this path).
        call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, MPI_MODE_RDONLY, MPI_INFO_NULL, ifile, ierr)

        do i = 1, n_stat
            call MPI_TYPE_CREATE_SUBARRAY(num_dims, sizes_glb, sizes_loc, start_stat, MPI_ORDER_FORTRAN, mpi_p, mpi_view, ierr)
            call MPI_TYPE_COMMIT(mpi_view, ierr)

            var_MOK = int(i, MPI_OFFSET_KIND)
            disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

            call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, mpi_view, 'native', MPI_INFO_NULL, ierr)
            call MPI_FILE_READ(ifile, q_stat_vf(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)

            call MPI_TYPE_FREE(mpi_view, ierr)
        end do

        call MPI_FILE_CLOSE(ifile, ierr)
#else
        found = .true.
#endif

    end subroutine s_read_lso_stat_file

    !> Read the simulation-written stage-1 filtered gas-mask (lso_mask_<t>.dat) into the INTERIOR of w_vf (caller allocates
    !! w_vf(1)%sf with ghost bounds). found is made collectively consistent so all ranks agree on whether to take the masked path.
    impure subroutine s_lso_sync_found(found)

        logical, intent(inout) :: found

#ifdef MFC_MPI
        integer :: found_int, found_min, ierr

        found_int = merge(1, 0, found)
        call MPI_ALLREDUCE(found_int, found_min, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
        found = (found_min == 1)
#endif

    end subroutine s_lso_sync_found

    impure subroutine s_read_lso_mask(w_vf, t_step, found)

        type(scalar_field), intent(inout) :: w_vf(1:1)
        integer, intent(in)               :: t_step
        logical, intent(out)              :: found
        type(scalar_field)                :: w_io(1:1)
        integer                           :: j, k, l

        allocate (w_io(1)%sf(0:m,0:n,0:p))
        call s_read_lso_stat_file(w_io, 1, t_step, found, 'lso_mask_')

        call s_lso_sync_found(found)

        if (found) then
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        w_vf(1)%sf(j, k, l) = w_io(1)%sf(j, k, l)
                    end do
                end do
            end do
        end if
        deallocate (w_io(1)%sf)

    end subroutine s_read_lso_mask

    !> Read the lso_stat binary for t_step and emit each variable into silo_hdf5_lso_stat/. No-op when the file is missing (e.g.
    !! step 0 in cfl_dt mode).
    impure subroutine s_save_lso_stat_data(t_step)

        integer, intent(in)             :: t_step
        type(scalar_field), allocatable :: q_stat_vf(:)
        logical                         :: found
        integer                         :: i

        if (n_lso_stat <= 0) return

        allocate (q_stat_vf(1:n_lso_stat))
        do i = 1, n_lso_stat
            allocate (q_stat_vf(i)%sf(0:m,0:n,0:p))
        end do

        call s_read_lso_stat_file(q_stat_vf, n_lso_stat, t_step, found)

        call s_lso_sync_found(found)

        if (found) call s_write_lso_stat_fields(q_stat_vf, t_step)

        do i = 1, n_lso_stat
            deallocate (q_stat_vf(i)%sf)
        end do
        deallocate (q_stat_vf)

    end subroutine s_save_lso_stat_data

    !> Compute the Euler-Lagrange closure fields from the simulation-written LSO stat binary (+ filtered conserved data and mask)
    !! and write them to silo_hdf5_lso_closure/. No-op when the stat file is missing.
    impure subroutine s_save_lso_closure_data(t_step)

        integer, intent(in)             :: t_step
        type(scalar_field), allocatable :: q_stat_vf(:), q_cls_vf(:)
        type(scalar_field)              :: w_vf(1:1)
        logical                         :: found, found_w
        integer                         :: i, n_cls

        if (n_lso_stat <= 0) return

        allocate (q_stat_vf(1:n_lso_stat))
        do i = 1, n_lso_stat
            allocate (q_stat_vf(i)%sf(0:m,0:n,0:p))
        end do
        call s_read_lso_stat_file(q_stat_vf, n_lso_stat, t_step, found)
        call s_lso_sync_found(found)

        if (found) then
            allocate (w_vf(1)%sf(0:m,0:n,0:p))
            call s_read_lso_mask(w_vf, t_step, found_w)
            if (.not. found_w) w_vf(1)%sf = 1._stp

            n_cls = f_lso_n_closure()
            allocate (q_cls_vf(1:n_cls))
            do i = 1, n_cls
                allocate (q_cls_vf(i)%sf(0:m,0:n,0:p))
            end do

            call s_compute_lso_closure_fields(q_stat_vf, q_cons_vf, w_vf, q_cls_vf)
            call s_write_lso_closure_fields(q_cls_vf, t_step)

            do i = 1, n_cls
                deallocate (q_cls_vf(i)%sf)
            end do
            deallocate (q_cls_vf)
            deallocate (w_vf(1)%sf)
        end if

        do i = 1, n_lso_stat
            deallocate (q_stat_vf(i)%sf)
        end do
        deallocate (q_stat_vf)

    end subroutine s_save_lso_closure_data

    !> Closure fields for the post-widened data (lso_pp_filter = T): apply the sigma2 pass to every stat product, then compute the
    !! closures against the post-filtered conserved state and weight w2 = filter2(w1).
    impure subroutine s_save_lso_pp_closure_data(t_step)

        integer, intent(in)             :: t_step
        type(scalar_field), allocatable :: q_io_vf(:), q_stat_vf(:), q_cls_vf(:)
        logical                         :: found
        integer                         :: i, c, j, k, l, n_cls

        if (n_lso_stat <= 0) return

        ! Read into interior-sized temporaries, then move into ghost-extended fields (the pp filter needs halo cells).
        allocate (q_io_vf(1:n_lso_stat), q_stat_vf(1:n_lso_stat))
        do i = 1, n_lso_stat
            allocate (q_io_vf(i)%sf(0:m,0:n,0:p))
            allocate (q_stat_vf(i)%sf(lbound(q_cons_vf(1)%sf, 1):ubound(q_cons_vf(1)%sf, 1),lbound(q_cons_vf(1)%sf, &
                      & 2):ubound(q_cons_vf(1)%sf, 2),lbound(q_cons_vf(1)%sf, 3):ubound(q_cons_vf(1)%sf, 3)))
        end do
        call s_read_lso_stat_file(q_io_vf, n_lso_stat, t_step, found)
        call s_lso_sync_found(found)

        if (found) then
            do i = 1, n_lso_stat
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            q_stat_vf(i)%sf(j, k, l) = q_io_vf(i)%sf(j, k, l)
                        end do
                    end do
                end do
            end do

            ! sigma2 pass over all stat products, in chunks of sys_size fields so the
            ! shared MPI halo buffers (sized for sys_size) are never exceeded.
            do c = 1, n_lso_stat, sys_size
                call s_apply_lso_pp_filter(q_stat_vf(c:min(c + sys_size - 1, n_lso_stat)))
            end do

            n_cls = f_lso_n_closure()
            allocate (q_cls_vf(1:n_cls))
            do i = 1, n_cls
                allocate (q_cls_vf(i)%sf(0:m,0:n,0:p))
            end do
            call s_compute_lso_closure_fields(q_stat_vf, q_cons_vf, q_lso_pp_w_vf, q_cls_vf)
            call s_write_lso_closure_fields(q_cls_vf, t_step)
            do i = 1, n_cls
                deallocate (q_cls_vf(i)%sf)
            end do
            deallocate (q_cls_vf)
        end if

        do i = 1, n_lso_stat
            deallocate (q_io_vf(i)%sf, q_stat_vf(i)%sf)
        end do
        deallocate (q_io_vf, q_stat_vf)

    end subroutine s_save_lso_pp_closure_data

    !> Closure fields for original data filtered in post_process (lso_pp_filter = T, lso_filter_wrt = F): the stat products in
    !! q_lso_pp_stat_vf were formed from the unfiltered state and filtered (p_main), the weight is the filtered gas mask.
    impure subroutine s_save_lso_pp_raw_closure_data(t_step)

        integer, intent(in)             :: t_step
        type(scalar_field), allocatable :: q_stat_vf(:), q_cls_vf(:)
        integer                         :: i, j, k, l, n_cls

        if (n_lso_stat <= 0) return

        allocate (q_stat_vf(1:n_lso_stat))
        do i = 1, n_lso_stat
            allocate (q_stat_vf(i)%sf(lbound(q_cons_vf(1)%sf, 1):ubound(q_cons_vf(1)%sf, 1),lbound(q_cons_vf(1)%sf, &
                      & 2):ubound(q_cons_vf(1)%sf, 2),lbound(q_cons_vf(1)%sf, 3):ubound(q_cons_vf(1)%sf, 3)))
            q_stat_vf(i)%sf = 0._stp
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_stat_vf(i)%sf(j, k, l) = q_lso_pp_stat_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
        n_cls = f_lso_n_closure()
        allocate (q_cls_vf(1:n_cls))
        do i = 1, n_cls
            allocate (q_cls_vf(i)%sf(0:m,0:n,0:p))
        end do
        call s_compute_lso_closure_fields(q_stat_vf, q_cons_vf, q_lso_pp_w_vf, q_cls_vf)
        call s_write_lso_closure_fields(q_cls_vf, t_step)
        do i = 1, n_cls
            deallocate (q_cls_vf(i)%sf)
        end do
        do i = 1, n_lso_stat
            deallocate (q_stat_vf(i)%sf)
        end do
        deallocate (q_cls_vf, q_stat_vf)

    end subroutine s_save_lso_pp_raw_closure_data

    !> Emit the closure fields into silo_hdf5_lso_closure/ with ParaView-friendly names. Field order must match
    !! s_compute_lso_closure_fields / f_lso_n_closure.
    impure subroutine s_write_lso_closure_fields(q_cls_vf, t_step)

        type(scalar_field), intent(in) :: q_cls_vf(:)
        integer, intent(in)            :: t_step
        character(LEN=name_len)        :: varname
        character(LEN=2)               :: tc(6)
        character(LEN=1)               :: dc(3)
        integer                        :: nd, nt, c, a, idx

        nd = num_dims
        nt = nd*(nd + 1)/2
        dc = (/'x', 'y', 'z'/)
        if (nd == 3) then
            tc(1:6) = (/'11', '12', '13', '22', '23', '33'/)
        else
            tc(1:3) = (/'11', '12', '22'/)
        end if

        call s_switch_to_lso_stat_dir('silo_hdf5_lso_closure')
        call s_open_formatted_database_file(t_step)
        call s_write_grid_to_formatted_database_file(t_step)

        idx = 0
        do c = 1, nt
            call put_next('R_sg_' // tc(c))
        end do
        do a = 1, nd
            call put_next('Q_T_' // dc(a))
        end do
        do a = 1, nd
            call put_next('E_ku_' // dc(a))
        end do
        do a = 1, nd
            call put_next('W_tau_u_' // dc(a))
        end do
        do c = 1, nt
            call put_next('R_mu_sg_' // tc(c))
        end do
        do a = 1, nd
            call put_next('R_lam_sg_' // dc(a))
        end do
        call put_next('T_tilde')
        do a = 1, nd
            call put_next('u_favre_' // dc(a))
        end do

        call s_close_formatted_database_file()
        call s_switch_output_dirs(.true.)

    contains

        impure subroutine put_next(name)

            character(LEN=*), intent(in) :: name

            idx = idx + 1
            varname = name
            call s_put_lso_var(q_cls_vf(idx))
            call s_write_variable_to_formatted_database_file(varname, t_step)

        end subroutine put_next

    end subroutine s_write_lso_closure_fields

    !> Write the post_process-computed LSO stat fields (q_lso_pp_stat_vf from m_lso_pp_filter) to silo_hdf5_lso_stat/. The fields
    !! are computed in-process by s_compute_lso_pp_stat_fields from the post_process-filtered conserved state, so no binary file
    !! read is required (unlike s_save_lso_stat_data).
    impure subroutine s_save_lso_pp_stat_data(t_step)

        integer, intent(in) :: t_step

        if (n_lso_stat <= 0) return

        call s_write_lso_stat_fields(q_lso_pp_stat_vf, t_step)

    end subroutine s_save_lso_pp_stat_data

    !> Copy a ghost-less LSO field into q_sf, replicating edge cells into the Silo pad region.
    impure subroutine s_put_lso_var(fld)

        type(scalar_field), intent(in) :: fld
        integer                        :: j, k, l

        do l = lbound(q_sf, 3), ubound(q_sf, 3)
            do k = lbound(q_sf, 2), ubound(q_sf, 2)
                do j = lbound(q_sf, 1), ubound(q_sf, 1)
                    q_sf(j, k, l) = real(fld%sf(min(max(j, 0), m), min(max(k, 0), n), min(max(l, 0), p)), wp)
                end do
            end do
        end do

    end subroutine s_put_lso_var

    !> Write the LSO statistical product fields in q_stat_vf to silo_hdf5_lso_stat/ and restore the LSO output directory.
    impure subroutine s_write_lso_stat_fields(q_stat_vf, t_step)

        type(scalar_field), intent(in) :: q_stat_vf(:)
        integer, intent(in)            :: t_step
        character(LEN=name_len)        :: varname
        character(LEN=2)               :: tc(6)
        integer                        :: c, nt

        nt = num_dims*(num_dims + 1)/2
        if (num_dims == 3) then
            tc(1:6) = (/'11', '12', '13', '22', '23', '33'/)
        else
            tc(1:3) = (/'11', '12', '22'/)
        end if

        call s_switch_to_lso_stat_dir()
        call s_open_formatted_database_file(t_step)
        call s_write_grid_to_formatted_database_file(t_step)

        ! phi_p
        varname = 'phi_p'
        call s_put_lso_var(q_stat_vf(lso_stat_phi_p_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)

        ! rho scalar (gas_mask * rho)
        varname = 'rho'
        call s_put_lso_var(q_stat_vf(lso_stat_rho_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)

        ! rhoke scalar (gas_mask * (mom1^2+mom2^2+mom3^2)/rho)
        varname = 'rho_ke'
        call s_put_lso_var(q_stat_vf(lso_stat_rhoke_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)

        ! phi_p * u_p
        varname = 'phi_p_up_x'
        call s_put_lso_var(q_stat_vf(lso_stat_up_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)
        if (num_dims >= 2) then
            varname = 'phi_p_up_y'
            call s_put_lso_var(q_stat_vf(lso_stat_up_beg + 1))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if
        if (num_dims == 3) then
            varname = 'phi_p_up_z'
            call s_put_lso_var(q_stat_vf(lso_stat_up_beg + 2))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if

        ! rho*u
        varname = 'rho_u_x'
        call s_put_lso_var(q_stat_vf(lso_stat_rhou_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)
        if (num_dims >= 2) then
            varname = 'rho_u_y'
            call s_put_lso_var(q_stat_vf(lso_stat_rhou_beg + 1))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if
        if (num_dims == 3) then
            varname = 'rho_u_z'
            call s_put_lso_var(q_stat_vf(lso_stat_rhou_beg + 2))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if

        ! rho*u*u upper triangle
        do c = 1, nt
            varname = 'rho_uu_' // tc(c)
            call s_put_lso_var(q_stat_vf(lso_stat_rhouu_beg + c - 1))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end do

        ! rho*u*|u|^2
        varname = 'rho_uke_x'
        call s_put_lso_var(q_stat_vf(lso_stat_rhouke_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)
        if (num_dims >= 2) then
            varname = 'rho_uke_y'
            call s_put_lso_var(q_stat_vf(lso_stat_rhouke_beg + 1))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if
        if (num_dims == 3) then
            varname = 'rho_uke_z'
            call s_put_lso_var(q_stat_vf(lso_stat_rhouke_beg + 2))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if

        ! rho*u*T
        varname = 'rho_uT_x'
        call s_put_lso_var(q_stat_vf(lso_stat_rhouT_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)
        if (num_dims >= 2) then
            varname = 'rho_uT_y'
            call s_put_lso_var(q_stat_vf(lso_stat_rhouT_beg + 1))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if
        if (num_dims == 3) then
            varname = 'rho_uT_z'
            call s_put_lso_var(q_stat_vf(lso_stat_rhouT_beg + 2))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if

        ! tau (same layout as rho_uu)
        do c = 1, nt
            varname = 'tau_' // tc(c)
            call s_put_lso_var(q_stat_vf(lso_stat_tau_beg + c - 1))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end do

        ! q_i
        varname = 'q_x'
        call s_put_lso_var(q_stat_vf(lso_stat_q_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)
        if (num_dims >= 2) then
            varname = 'q_y'
            call s_put_lso_var(q_stat_vf(lso_stat_q_beg + 1))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if
        if (num_dims == 3) then
            varname = 'q_z'
            call s_put_lso_var(q_stat_vf(lso_stat_q_beg + 2))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if

        ! (tau u)_i
        varname = 'rho_tau_u_x'
        call s_put_lso_var(q_stat_vf(lso_stat_rhotau_u_beg))
        call s_write_variable_to_formatted_database_file(varname, t_step)
        if (num_dims >= 2) then
            varname = 'rho_tau_u_y'
            call s_put_lso_var(q_stat_vf(lso_stat_rhotau_u_beg + 1))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if
        if (num_dims == 3) then
            varname = 'rho_tau_u_z'
            call s_put_lso_var(q_stat_vf(lso_stat_rhotau_u_beg + 2))
            call s_write_variable_to_formatted_database_file(varname, t_step)
        end if

        call s_close_formatted_database_file()

        ! Restore LSO directory for subsequent passes
        call s_switch_output_dirs(.true.)

    end subroutine s_write_lso_stat_fields

    !> Transpose 3-D complex data from x-pencil to y-pencil layout via MPI_Alltoall.
    subroutine s_mpi_transpose_x2y

        complex(c_double_complex), allocatable :: sendbuf(:), recvbuf(:)
        integer                                :: dest_rank, src_rank
        integer                                :: i, j, k, l

#ifdef MFC_MPI
        allocate (sendbuf(Nx*Nyloc*Nzloc))
        allocate (recvbuf(Nx*Nyloc*Nzloc))

        do dest_rank = 0, num_procs_y - 1
            do l = 1, Nzloc
                do k = 1, Nyloc
                    do j = 1, Nxloc
                        sendbuf(j + (k - 1)*Nxloc + (l - 1)*Nxloc*Nyloc + dest_rank*Nxloc*Nyloc*Nzloc) = data_cmplx(j &
                                & + dest_rank*Nxloc, k, l)
                    end do
                end do
            end do
        end do

        call MPI_Alltoall(sendbuf, Nxloc*Nyloc*Nzloc, MPI_C_DOUBLE_COMPLEX, recvbuf, Nxloc*Nyloc*Nzloc, MPI_C_DOUBLE_COMPLEX, &
                          & MPI_COMM_CART12, ierr)

        do src_rank = 0, num_procs_y - 1
            do l = 1, Nzloc
                do k = 1, Nyloc
                    do j = 1, Nxloc
                        data_cmplx_y(j, k + src_rank*Nyloc, &
                                     & l) = recvbuf(j + (k - 1)*Nxloc + (l - 1)*Nxloc*Nyloc + src_rank*Nxloc*Nyloc*Nzloc)
                    end do
                end do
            end do
        end do

        deallocate (sendbuf)
        deallocate (recvbuf)
#endif

    end subroutine s_mpi_transpose_x2y

    !> Transpose 3-D complex data from y-pencil to z-pencil layout via MPI_Alltoall.
    subroutine s_mpi_transpose_y2z

        complex(c_double_complex), allocatable :: sendbuf(:), recvbuf(:)
        integer                                :: dest_rank, src_rank
        integer                                :: j, k, l

#ifdef MFC_MPI
        allocate (sendbuf(Ny*Nxloc*Nzloc))
        allocate (recvbuf(Ny*Nxloc*Nzloc))

        do dest_rank = 0, num_procs_z - 1
            do l = 1, Nzloc
                do j = 1, Nxloc
                    do k = 1, Nyloc2
                        sendbuf(k + (j - 1)*Nyloc2 + (l - 1)*(Nyloc2*Nxloc) + dest_rank*Nyloc2*Nxloc*Nzloc) = data_cmplx_y(j, &
                                & k + dest_rank*Nyloc2, l)
                    end do
                end do
            end do
        end do

        call MPI_Alltoall(sendbuf, Nyloc2*Nxloc*Nzloc, MPI_C_DOUBLE_COMPLEX, recvbuf, Nyloc2*Nxloc*Nzloc, MPI_C_DOUBLE_COMPLEX, &
                          & MPI_COMM_CART13, ierr)

        do src_rank = 0, num_procs_z - 1
            do l = 1, Nzloc
                do j = 1, Nxloc
                    do k = 1, Nyloc2
                        data_cmplx_z(j, k, &
                                     & l + src_rank*Nzloc) = recvbuf(k + (j - 1)*Nyloc2 + (l - 1)*(Nyloc2*Nxloc) &
                                     & + src_rank*Nyloc2*Nxloc*Nzloc)
                    end do
                end do
            end do
        end do

        deallocate (sendbuf)
        deallocate (recvbuf)
#endif

    end subroutine s_mpi_transpose_y2z

    !> Initialize all post-process sub-modules, set up I/O pointers, and prepare FFTW plans and MPI communicators.
    impure subroutine s_initialize_modules

        integer :: size_n(1), inembed(1), onembed(1)

        call s_initialize_global_parameters_module()
        if (bubbles_euler .or. bubbles_lagrange) then
            call s_initialize_bubbles_model()
        end if
        if (num_procs > 1) then
            call s_initialize_mpi_proxy_module()
            call s_initialize_mpi_common_module()
        end if
        call s_initialize_boundary_common_module()
        call s_initialize_variables_conversion_module()
        call s_initialize_data_input_module()
        call s_initialize_derived_variables_module()
        call s_initialize_data_output_module()
        if (lso_pp_filter) call s_initialize_lso_pp_filter_module()

        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
        end if

#ifdef MFC_MPI
        if (fft_wrt) then
            num_procs_x = (m_glb + 1)/(m + 1)
            num_procs_y = (n_glb + 1)/(n + 1)
            num_procs_z = (p_glb + 1)/(p + 1)

            Nx = m_glb + 1
            Ny = n_glb + 1
            Nz = p_glb + 1

            Nxloc = (m_glb + 1)/num_procs_y
            Nyloc = n + 1
            Nyloc2 = (n_glb + 1)/num_procs_z
            Nzloc = p + 1

            Nf = max(Nx, Ny, Nz)

            @:ALLOCATE(data_in(Nx*Nyloc*Nzloc))
            @:ALLOCATE(data_out(Nx*Nyloc*Nzloc))

            @:ALLOCATE(data_cmplx(Nx, Nyloc, Nzloc))
            @:ALLOCATE(data_cmplx_y(Nxloc, Ny, Nzloc))
            @:ALLOCATE(data_cmplx_z(Nxloc, Nyloc2, Nz))

            @:ALLOCATE(En_real(Nxloc, Nyloc2, Nz))
            @:ALLOCATE(En(Nf))

            size_n(1) = Nx
            inembed(1) = Nx
            onembed(1) = Nx

            fwd_plan_x = fftw_plan_many_dft(1, size_n, Nyloc*Nzloc, data_in, inembed, 1, Nx, data_out, onembed, 1, Nx, &
                                            & FFTW_FORWARD, FFTW_MEASURE)

            size_n(1) = Ny
            inembed(1) = Ny
            onembed(1) = Ny

            fwd_plan_y = fftw_plan_many_dft(1, size_n, Nxloc*Nzloc, data_out, inembed, 1, Ny, data_in, onembed, 1, Ny, &
                                            & FFTW_FORWARD, FFTW_MEASURE)

            size_n(1) = Nz
            inembed(1) = Nz
            onembed(1) = Nz

            fwd_plan_z = fftw_plan_many_dft(1, size_n, Nxloc*Nyloc2, data_in, inembed, 1, Nz, data_out, onembed, 1, Nz, &
                                            & FFTW_FORWARD, FFTW_MEASURE)

            call MPI_CART_CREATE(MPI_COMM_WORLD, 3, (/num_procs_x, num_procs_y, num_procs_z/), (/.true., .true., .true./), &
                                 & .false., MPI_COMM_CART, ierr)
            call MPI_CART_COORDS(MPI_COMM_CART, proc_rank, 3, cart3d_coords, ierr)

            call MPI_Cart_SUB(MPI_COMM_CART, (/.true., .true., .false./), MPI_COMM_CART12, ierr)
            call MPI_COMM_RANK(MPI_COMM_CART12, proc_rank12, ierr)
            call MPI_CART_COORDS(MPI_COMM_CART12, proc_rank12, 2, cart2d12_coords, ierr)

            call MPI_Cart_SUB(MPI_COMM_CART, (/.true., .false., .true./), MPI_COMM_CART13, ierr)
            call MPI_COMM_RANK(MPI_COMM_CART13, proc_rank13, ierr)
            call MPI_CART_COORDS(MPI_COMM_CART13, proc_rank13, 2, cart2d13_coords, ierr)
        end if
#endif

    end subroutine s_initialize_modules

    !> Perform a distributed forward 3-D FFT using pencil decomposition with FFTW and MPI transposes.
    subroutine s_mpi_FFT_fwd

        integer :: j, k, l

#ifdef MFC_MPI
        do l = 1, Nzloc
            do k = 1, Nyloc
                do j = 1, Nx
                    data_in(j + (k - 1)*Nx + (l - 1)*Nx*Nyloc) = data_cmplx(j, k, l)
                end do
            end do
        end do

        call fftw_execute_dft(fwd_plan_x, data_in, data_out)

        do l = 1, Nzloc
            do k = 1, Nyloc
                do j = 1, Nx
                    data_cmplx(j, k, l) = data_out(j + (k - 1)*Nx + (l - 1)*Nx*Nyloc)
                end do
            end do
        end do

        call s_mpi_transpose_x2y !!Change Pencil from data_cmplx to data_cmpx_y

        do l = 1, Nzloc
            do k = 1, Nxloc
                do j = 1, Ny
                    data_out(j + (k - 1)*Ny + (l - 1)*Ny*Nxloc) = data_cmplx_y(k, j, l)
                end do
            end do
        end do

        call fftw_execute_dft(fwd_plan_y, data_out, data_in)

        do l = 1, Nzloc
            do k = 1, Nxloc
                do j = 1, Ny
                    data_cmplx_y(k, j, l) = data_in(j + (k - 1)*Ny + (l - 1)*Ny*Nxloc)
                end do
            end do
        end do

        call s_mpi_transpose_y2z !!Change Pencil from data_cmplx_y to data_cmpx_z

        do l = 1, Nyloc2
            do k = 1, Nxloc
                do j = 1, Nz
                    data_in(j + (k - 1)*Nz + (l - 1)*Nz*Nxloc) = data_cmplx_z(k, l, j)
                end do
            end do
        end do

        call fftw_execute_dft(fwd_plan_z, data_in, data_out)

        do l = 1, Nyloc2
            do k = 1, Nxloc
                do j = 1, Nz
                    data_cmplx_z(k, l, j) = data_out(j + (k - 1)*Nz + (l - 1)*Nz*Nxloc)
                end do
            end do
        end do
#endif

    end subroutine s_mpi_FFT_fwd

    !> Set up the MPI environment, read and broadcast user inputs, and decompose the computational domain.
    impure subroutine s_initialize_mpi_domain

        num_dims = 1 + min(1, n) + min(1, p)

        call s_mpi_initialize()

        if (proc_rank == 0) then
            call s_assign_default_values_to_user_inputs()
            call s_read_input_file()
            call s_check_input_file()

            print '(" Post-processing a ", I0, "x", I0, "x", I0, " case on ", I0, " rank(s)")', m, n, p, num_procs
        end if

        call s_mpi_bcast_user_inputs()
        call s_initialize_parallel_io()
        call s_mpi_decompose_computational_domain()
        call s_check_inputs_fft()

    end subroutine s_initialize_mpi_domain

    !> Destroy FFTW plans, free MPI communicators, and finalize all post-process sub-modules.
    impure subroutine s_finalize_modules

        s_read_data_files => null()

        if (fft_wrt) then
            if (c_associated(fwd_plan_x)) call fftw_destroy_plan(fwd_plan_x)
            if (c_associated(fwd_plan_y)) call fftw_destroy_plan(fwd_plan_y)
            if (c_associated(fwd_plan_z)) call fftw_destroy_plan(fwd_plan_z)
            if (allocated(data_in)) deallocate (data_in)
            if (allocated(data_out)) deallocate (data_out)
            if (allocated(data_cmplx)) deallocate (data_cmplx)
            if (allocated(data_cmplx_y)) deallocate (data_cmplx_y)
            if (allocated(data_cmplx_z)) deallocate (data_cmplx_z)
            if (allocated(En_real)) deallocate (En_real)
            if (allocated(En)) deallocate (En)
            call fftw_cleanup()
        end if

#ifdef MFC_MPI
        if (fft_wrt) then
            if (MPI_COMM_CART12 /= MPI_COMM_NULL) call MPI_Comm_free(MPI_COMM_CART12, ierr)
            if (MPI_COMM_CART13 /= MPI_COMM_NULL) call MPI_Comm_free(MPI_COMM_CART13, ierr)
            if (MPI_COMM_CART /= MPI_COMM_NULL) call MPI_Comm_free(MPI_COMM_CART, ierr)
        end if
#endif

        if (lso_pp_filter) call s_finalize_lso_pp_filter_module()
        call s_finalize_data_output_module()
        call s_finalize_derived_variables_module()
        call s_finalize_data_input_module()
        call s_finalize_variables_conversion_module()
        if (num_procs > 1) then
            call s_finalize_mpi_proxy_module()
            call s_finalize_mpi_common_module()
        end if
        call s_finalize_global_parameters_module()

        call s_mpi_finalize()

    end subroutine s_finalize_modules

end module m_start_up
