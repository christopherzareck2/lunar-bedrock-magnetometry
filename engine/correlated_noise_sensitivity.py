"""Along-track correlated noise for the correlated-noise and drift tests.

Per axis, with the total variance fixed at (0.50 nT)^2 for every correlation length L:

    n(s) = w(s) + c(s)
    w ~ N(0, 0.05^2) white sensor noise
    c   stationary Ornstein-Uhlenbeck (AR(1)) process, variance 0.50^2 - 0.05^2, correlation exp(-ds / L)

L -> 0 gives the i.i.d. case; L -> infinity gives a constant offset per survey.
"""
import numpy as np

from noise_budget import ROVER_ADOPTED_NT, ROVER_INSTRUMENT_NT   # Table 1 (engine/noise_budget.py)

SPACING = 0.01            # m, along-track sampling
SIGMA_TOT = ROVER_ADOPTED_NT      # nT, total per-axis std (0.50)
SIGMA_W = ROVER_INSTRUMENT_NT     # nT, white sensor noise (0.05)
SIGMA_C = float(np.sqrt(SIGMA_TOT ** 2 - SIGMA_W ** 2))   # correlated std (~0.4975 nT)


def along_track_ds(sxyz):
    """Positive along-track step lengths between consecutive stations (m)."""
    s = np.asarray(sxyz, float)
    return np.linalg.norm(np.diff(s[:, :2], axis=0), axis=1)   # (nR-1,), xy path


def make_noise(sxyz, L, n_draws, rng, n_axes=3):
    """Correlated noise, shape (n_draws, nR*3), station-major [bx, by, bz per station], variance SIGMA_TOT^2.

    L = None gives white noise (SIGMA_TOT); L = np.inf gives a constant offset per survey.
    """
    nR = len(sxyz)
    if L is None:
        n = rng.normal(0.0, SIGMA_TOT, size=(n_draws, nR, n_axes))
        return n.reshape(n_draws, nR * n_axes)
    ds = along_track_ds(sxyz)                                   # (nR-1,)
    rho = np.ones_like(ds) if np.isinf(L) else np.exp(-ds / L)  # lag-1 correlation per step
    w = rng.normal(0.0, SIGMA_W, size=(n_draws, nR, n_axes))
    c = np.empty((n_draws, nR, n_axes))
    c[:, 0, :] = SIGMA_C * rng.normal(size=(n_draws, n_axes))
    innov = SIGMA_C * rng.normal(size=(n_draws, nR, n_axes))
    for k in range(1, nR):
        rk = rho[k - 1]
        c[:, k, :] = rk * c[:, k - 1, :] + np.sqrt(1.0 - rk ** 2) * innov[:, k, :]
    return (w + c).reshape(n_draws, nR * n_axes)


def effective_sample_count(sxyz, L):
    """Effective number of independent samples per axis for the full noise autocovariance.

    tau = 1 + 2 * sum_{k>=1} (1-k/N) * gamma_k/gamma_0, N_eff = N/tau, with gamma_0 = SIGMA_W^2 + SIGMA_C^2
    and gamma_k = SIGMA_C^2 * exp(-|dist_k|/L) (white noise at lag 0 only). Uses the median along-track step.
    """
    N = len(sxyz)
    g0 = SIGMA_W ** 2 + SIGMA_C ** 2
    if L is None:
        return {"N": N, "N_eff": float(N), "N_eff_over_N": 1.0, "rho_lag1": 0.0}
    ds = along_track_ds(sxyz)
    d = float(np.median(ds[ds > 0])) if np.any(ds > 0) else SPACING
    if np.isinf(L):
        rho1, gk = 1.0, np.full(N - 1, SIGMA_C ** 2)
    else:
        rho1 = float(np.exp(-d / L))
        k = np.arange(1, N)
        gk = SIGMA_C ** 2 * np.exp(-k * d / L)
    kk = np.arange(1, N)
    tau = 1.0 + 2.0 * np.sum((1.0 - kk / N) * gk / g0)
    return {"N": N, "N_eff": float(N / tau), "N_eff_over_N": float(1.0 / tau), "rho_lag1": rho1}
