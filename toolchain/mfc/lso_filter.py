"""LSO (Least-Squares Optimized) variable-weight 9-point Gaussian filter.

LSO denotes the least-squares optimization used to fit repeated 9-point stencils to a
target Gaussian transfer function.  Each pass has a symmetric FIR response

    H(xi) = a0 + 2*a1*cos(xi) + 2*a2*cos(2*xi) + 2*a3*cos(3*xi) + 2*a4*cos(4*xi),

normalized so H(0) = 1.  Per-pass weights are fitted by block coordinate descent (BCD)
against G(xi) = exp(-0.5 * sigma^2 * xi^2) with sigma = filter_sigma / dx.
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

# Must stay in sync with lso_max_passes in the Fortran side.
LSO_MAX_PASSES: int = 60
# A single cascade is only well conditioned up to ~40 cells (finite-precision wall at ~45),
# so post_process splits a wider target across two same-grid stages.
PP_SPLIT_CELLS: float = 30.0
PP_STAGE_MAX_CELLS: float = 40.0
_CONV_TOL: float = 1e-3
_N_XI: int = 600


def _filter_tf(coeffs: np.ndarray, xi: np.ndarray) -> np.ndarray:
    a0, a1, a2, a3, a4 = coeffs
    return a0 + 2.0 * a1 * np.cos(xi) + 2.0 * a2 * np.cos(2.0 * xi) + 2.0 * a3 * np.cos(3.0 * xi) + 2.0 * a4 * np.cos(4.0 * xi)


def _taylor_init(sigma_per_pass: float) -> np.ndarray:
    """Taylor-matched initial coefficients for one pass: matches G to O(xi^8)
    by equating H at xi = 1, 2, 3, 4 with normalization H(0) = 1."""
    r0 = min(sigma_per_pass, 1.38)
    r02 = r0**2

    rows = []
    for m in [1, 2, 3, 4]:
        row = [
            1.0,
            -(m * r02) / 2.0,
            (m * r02) ** 2 / 8.0,
            -((m * r02) ** 3) / 48.0,
            (m * r02) ** 4 / 384.0,
        ]
        rows.append(row)
    rows.append([1.0, 2.0, 2.0, 2.0, 2.0])  # normalization H(0)=1

    A = np.array(rows, dtype=float)
    b = np.array(
        [np.exp(-0.5 * m**2 * r02) for m in [1, 2, 3, 4]] + [1.0],
        dtype=float,
    )
    try:
        return np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        # Identity pass: H(xi) = 1 everywhere.  Satisfies H(0) = 1.
        return np.array([1.0, 0.0, 0.0, 0.0, 0.0])


def design_lso_filter(
    sigma_target: float,
    n_passes: int,
    tol_inner: float = 1e-8,
    n_max: int = 500,
    n_xi: int = _N_XI,
) -> List[Tuple[float, ...]]:
    """Block coordinate descent fit of n_passes 9-point filters against
    the target Gaussian G. Returns n_passes tuples (a0..a4)."""
    xi = np.linspace(0.0, np.pi, n_xi)
    G = np.exp(-0.5 * sigma_target**2 * xi**2)
    Phi = np.column_stack(
        [
            np.ones(n_xi),
            2.0 * np.cos(xi),
            2.0 * np.cos(2.0 * xi),
            2.0 * np.cos(3.0 * xi),
            2.0 * np.cos(4.0 * xi),
        ]
    )
    c = np.array([1.0, 2.0, 2.0, 2.0, 2.0])  # H(0) = 1

    sigma_per_pass = sigma_target / np.sqrt(max(n_passes, 1))
    betas = [_taylor_init(sigma_per_pass).copy() for _ in range(n_passes)]

    for _ in range(n_max):
        prev = [b.copy() for b in betas]

        for k in range(n_passes):
            H_rest = np.ones(n_xi)
            for j, b in enumerate(betas):
                if j != k:
                    H_rest *= Phi @ b

            Phi_w = H_rest[:, None] * Phi
            A_w = Phi_w.T @ Phi_w + 1e-9 * np.eye(5)
            rhs = Phi_w.T @ G
            A_inv = np.linalg.inv(A_w)
            beta_u = A_inv @ rhs
            Ainv_c = A_inv @ c
            lam = (c @ beta_u - 1.0) / (c @ Ainv_c + 1e-12)
            betas[k] = beta_u - lam * Ainv_c

        delta = max(np.linalg.norm(betas[k] - prev[k]) for k in range(n_passes))
        if delta < tol_inner:
            break

    return [tuple(float(x) for x in b) for b in betas]


def find_min_lso_passes(
    sigma_target: float,
    conv_tol: float = _CONV_TOL,
    max_passes: int = LSO_MAX_PASSES,
    n_xi: int = _N_XI,
) -> Tuple[int, float, List[Tuple[float, ...]]]:
    """Smallest n_passes such that the frequency-domain L2 error
    sqrt(mean((H_cum - G)^2)) on [0, pi] is below conv_tol."""
    xi = np.linspace(0.0, np.pi, n_xi)
    G = np.exp(-0.5 * sigma_target**2 * xi**2)

    err = float("inf")
    coeffs: List[Tuple[float, ...]] = []

    for n_passes in range(1, max_passes + 1):
        coeffs = design_lso_filter(sigma_target, n_passes, n_xi=n_xi)
        H_cum = np.ones(n_xi)
        for ck in coeffs:
            H_cum *= _filter_tf(np.array(ck), xi)
        err = float(np.sqrt(np.mean((H_cum - G) ** 2)))
        if err < conv_tol:
            return n_passes, err, coeffs

    return max_passes, err, coeffs


def _worst_partial_product(coeffs_list: List[Tuple[float, ...]], xi: np.ndarray) -> float:
    """Peak magnitude of the running partial product over the pass sequence (1.0 = never amplifies)."""
    P = np.ones_like(xi)
    worst = 1.0
    for c in coeffs_list:
        P = P * _filter_tf(np.asarray(c, dtype=float), xi)
        worst = max(worst, float(np.max(np.abs(P))))
    return worst


def reorder_passes_for_conditioning(coeffs_list: List[Tuple[float, ...]], xi: np.ndarray) -> List[Tuple[float, ...]]:
    """Greedily reorder passes to minimize the running partial product. The composed
    transfer function is unchanged; only the conditioning of the sequential application improves."""
    if len(coeffs_list) <= 1:
        return list(coeffs_list)

    Hs = [_filter_tf(np.asarray(c, dtype=float), xi) for c in coeffs_list]
    remaining = list(range(len(coeffs_list)))
    order: List[int] = []
    P = np.ones_like(xi)
    while remaining:
        best = min(remaining, key=lambda i: float(np.max(np.abs(P * Hs[i]))))
        order.append(best)
        P = P * Hs[best]
        remaining.remove(best)

    return [coeffs_list[i] for i in order]


def compute_lso_params(
    d_p: float,
    dx: float,
    dy: float,
    dz: float,
    filter_sigma: float,
    conv_tol: float = _CONV_TOL,
    max_passes: int = LSO_MAX_PASSES,
) -> dict:
    """Entry point used by case.py. Pass dy = 0 (1D) or dz = 0 (2D) to skip
    those directions. Returns {lso_n_passes_x/y/z, lso_a_x/y/z}."""
    directions = [("x", dx)]
    if dy > 0.0:
        directions.append(("y", dy))
    if dz > 0.0:
        directions.append(("z", dz))

    result: dict = {}
    xi = np.linspace(0.0, np.pi, _N_XI)
    for tag, d in directions:
        sigma_target = filter_sigma / d
        # Finite-precision stability limit of the repeated 9-point cascade.
        if sigma_target > 45.0:
            raise ValueError(
                f"LSO filter: sigma = {sigma_target:.1f} cells in {tag} exceeds the "
                f"stability limit (~45) of the repeated 9-point cascade. Filter at a "
                f"smaller in-situ width and compose to the target in post_process "
                f"(lso_filter_sigma_in/lso_filter_sigma_target), or coarsen the grid."
            )
        n_passes, err, coeffs = find_min_lso_passes(sigma_target, conv_tol=conv_tol, max_passes=max_passes)
        coeffs = reorder_passes_for_conditioning(coeffs, xi)
        gain = _worst_partial_product(coeffs, xi)
        result[f"lso_n_passes_{tag}"] = n_passes
        result[f"lso_a_{tag}"] = coeffs
        n_res_note = f"N_res={d_p / d:.2f}, " if d_p > 0.0 else ""
        print(f"  [LSO] {tag}-dir: {n_res_note}sigma_target={sigma_target:.2f} cells, N_passes={n_passes}, L2_err={err:.2e}, peak_partial_gain={gain:.1e}")

    if dy <= 0.0:
        result["lso_n_passes_y"] = 0
        result["lso_a_y"] = []
    if dz <= 0.0:
        result["lso_n_passes_z"] = 0
        result["lso_a_z"] = []

    return result


def lso_namelist_lines(lso_params: dict, prefix: str = "lso_") -> str:
    """Format the dict from compute_lso_params() as &user_inputs namelist lines.
    lso_a_x is shape (5, LSO_MAX_PASSES), column-major, zero-padded.

    The dict keys from compute_lso_params() are always in lso_ form (lso_n_passes_x,
    lso_a_x). ``prefix`` sets the namelist variable prefix in the emitted lines: "lso_" for
    the in-situ filter, "lso_pp_" for the post_process additional filter."""
    lines = []

    for tag in ["x", "y", "z"]:
        n_passes = lso_params.get(f"lso_n_passes_{tag}", 0)
        coeffs = lso_params.get(f"lso_a_{tag}", [])

        lines.append(f"{prefix}n_passes_{tag} = {n_passes}")

        if n_passes > 0 and coeffs:
            # Fortran column-major: coeff index varies fastest.
            flat = []
            for i_pass in range(LSO_MAX_PASSES):
                for j_coeff in range(5):
                    if i_pass < len(coeffs):
                        flat.append(f"{coeffs[i_pass][j_coeff]:.17e}")
                    else:
                        flat.append("0.0e0")
            lines.append(f"{prefix}a_{tag} = {' '.join(flat)}")

    return "\n".join(lines) + "\n"
