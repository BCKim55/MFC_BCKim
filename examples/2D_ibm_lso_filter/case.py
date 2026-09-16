import json
import math

# 2D compressible flow over a single immersed cylinder, with the LSO Gaussian
# filter coarse-graining the conserved variables in situ at every save step.
# See README.md for how to run and post-process this case.

# Fluid (air) and flow conditions
gam_a = 1.4
R_gas = 287.0
mu = 1.81e-5
T0 = 293.0
p0 = 101325.0
rho0 = p0 / (R_gas * T0)
c0 = math.sqrt(gam_a * p0 / rho0)

# Particle: one cylinder of diameter d_p resolved by 10 cells
d_p = 5.0e-4
Nx, Ny = 200, 120
Lx, Ly = 1.0e-2, 6.0e-3
dx = Lx / (Nx + 1)

# Streamwise body force sized to reach Ma ~ 0.1 within a few flow-through times
u_target = 0.1 * c0
t_flow = Lx / u_target
g_x = u_target / (0.5 * t_flow)

# LSO filter width: 1.5 particle diameters (= 15 cells), written on a grid
# coarsened 5x. The pass counts and stencil weights are derived by the
# toolchain; only the four lso_* switches below are user inputs.
filter_sigma = 1.5 * d_p

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": -1.5e-3,
            "x_domain%end": -1.5e-3 + Lx,
            "y_domain%beg": -0.5 * Ly,
            "y_domain%end": 0.5 * Ly,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "cyl_coord": "F",
            "cfl_adap_dt": "T",
            "cfl_target": 0.5,
            "n_start": 0,
            "t_stop": 3.0 * t_flow,
            "t_save": 0.6 * t_flow,
            # Simulation Algorithm Parameters
            "num_patches": 1,
            "num_fluids": 1,
            "model_eqns": 2,
            "alt_soundspeed": "F",
            "mpp_lim": "F",
            "mixture_err": "T",
            "time_stepper": 3,
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "T",
            "weno_avg": "F",
            "avg_state": 2,
            "mapped_weno": "T",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "viscous": "T",
            "bc_x%beg": -1,
            "bc_x%end": -1,
            "bc_y%beg": -1,
            "bc_y%end": -1,
            # Body force driving the flow through the periodic box
            "bf_x": "T",
            "g_x": g_x,
            "k_x": 0.0,
            "w_x": 0.0,
            "p_x": 0.0,
            # Formatted Database Files Structure Parameters
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "fd_order": 1,
            # Patch: quiescent air filling the whole domain
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": -1.5e-3 + 0.5 * Lx,
            "patch_icpp(1)%y_centroid": 0.0,
            "patch_icpp(1)%length_x": Lx,
            "patch_icpp(1)%length_y": Ly,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p0,
            "patch_icpp(1)%alpha_rho(1)": rho0,
            "patch_icpp(1)%alpha(1)": 1.0,
            # Immersed cylinder
            "ib": "T",
            "num_ibs": 1,
            "patch_ib(1)%geometry": 2,
            "patch_ib(1)%x_centroid": 1.0e-3,
            "patch_ib(1)%y_centroid": 0.0,
            "patch_ib(1)%radius": 0.5 * d_p,
            "patch_ib(1)%slip": "F",
            # LSO Gaussian filter (in-situ): filter the conserved variables at
            # every save step and write them as lustre_lso_*.dat on a 5x
            # coarser grid, next to the unfiltered restart data
            "lso_filter": "T",
            "lso_filter_wrt": "T",
            "filter_sigma": filter_sigma,
            "lso_down_sample_factor": 5,
            # Fluid Physical Parameters
            "fluid_pp(1)%gamma": 1.0e00 / (gam_a - 1.0e00),
            "fluid_pp(1)%pi_inf": 0.0,
            "fluid_pp(1)%Re(1)": 1.0 / mu,
        }
    )
)
