# 2D IBM cylinder with the LSO Gaussian filter

Compressible flow over a single immersed cylinder (2D, fully periodic, body-force
driven), coarse-grained in situ by the LSO (least-squares optimized) Gaussian
filter. This is the minimal working setup for extracting mesoscale
(LES / Euler-Lagrange) fields from a particle-resolved simulation.

## What the filter does

At every save step the simulation applies a Gaussian low-pass filter of standard
deviation `filter_sigma` (here 1.5 particle diameters = 15 cells) to the conserved
variables and writes the result next to the unfiltered restart data, on a grid
coarsened by `lso_down_sample_factor` (here 5x). With immersed boundaries the
filter is a mask-normalized convolution, so the non-physical state inside the
cylinder never contaminates the fluid average.

The filter itself is a cascade of 9-point FIR passes whose composed transfer
function matches the target Gaussian to a frequency-domain RMS error below 1e-3.
The pass counts and stencil weights are designed automatically by the toolchain
from `filter_sigma` and the grid; for widths beyond ~45 cells the toolchain splits
the filter into a two-stage pyramid around the decimation (which is why
`parallel_io` is required here).

## Running

```shell
./mfc.sh run examples/2D_ibm_lso_filter/case.py -n 4
```

This runs pre_process and simulation (about two minutes on four ranks) and leaves,
under `examples/2D_ibm_lso_filter/restart_data/`:

| file | content |
| ---: | :--- |
| `lustre_<t>.dat` | unfiltered conserved variables (as always) |
| `lustre_lso_<t>.dat` | filtered conserved variables on the 40x24 coarse grid |
| `lustre_lso_mask_<t>.dat` | filtered gas mask `w` (immersed-boundary cases) |
| `lustre_lso_{x,y}_cb.dat` | coarse-grid cell-boundary coordinates |

Then convert the filtered data to Silo for visualization:

```shell
./mfc.sh run examples/2D_ibm_lso_filter/case.py -n 4 -t post_process
```

which writes `silo_hdf5_lso/` (open `silo_hdf5_lso/root/collection_*.silo` in
ParaView or VisIt). The initial condition is written unfiltered by pre_process, so
post_process skips step 0 with a warning -- filtered output starts at the first
simulation save.

## The two offline modes

**Widen the filter after the fact.** The saved mask makes a further Gaussian pass
compose exactly with the in-situ one, so the final width can be chosen at
post-processing time. To re-filter this run's 15-cell data to 30 cells, add

```python
"lso_pp_filter": "T",
"lso_filter_sigma_target": 2 * filter_sigma,
```

and rerun the post_process command above.

**Filter original data from scratch.** Without the in-situ filter (`lso_filter_wrt`
absent or `F`), the same two parameters make post_process filter the full-resolution
restart data directly -- useful for runs that were made without the filter. Widths
up to ~45 cells per direction are supported on a single grid.

## Parameter summary

| parameter | meaning |
| ---: | :--- |
| `lso_filter` | enable the in-situ LSO filter |
| `lso_filter_wrt` | write the filtered fields at every save step |
| `filter_sigma` | Gaussian standard deviation, physical units (default `patch_ib(1)%radius`) |
| `lso_down_sample_factor` | coarsening factor of the filtered output grid |
| `lso_pp_filter` | apply a further LSO filter in post_process |
| `lso_filter_sigma_target` | total width the post_process filter brings the data to |
| `lso_filter_sigma_in` | width of the input data when widening (defaults to `filter_sigma`) |

All other `lso_*` namelist entries (pass counts, stencil weights) are derived by
the toolchain and injected into the input files; do not set them by hand.

## Reading the filtered binaries directly

The MPI-IO files are headerless little-endian float64, one contiguous block per
conserved variable, x-fastest:

```python
import numpy as np

nxc, nyc, nvars = 40, 24, 5   # (m+1)//5, (n+1)//5, sys_size
q = np.fromfile("restart_data/lustre_lso_5.dat").reshape(nvars, nyc, nxc)
w = np.fromfile("restart_data/lustre_lso_mask_5.dat").reshape(nyc, nxc)
x = np.fromfile("restart_data/lustre_lso_x_cb.dat")   # nxc+1 cell boundaries
```

`q` holds the mask-normalized filtered state (finite everywhere, including inside
the cylinder); `w` is the filtered gas mask, needed to compose further filtering
offline or to phase-weight averages.
