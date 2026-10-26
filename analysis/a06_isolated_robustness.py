"""a06 - Isolated outcrop: robustness tests (paper Table 3, Sections 4.2 and 5.1).

Geometry, noise and solver as in a02 (4 x 4 x 1 m prism on the ground, sensor 0.5 m above the ground, 1 m
standoff, 1 cm sampling, sigma_B = 0.50 nT per axis, closed-form weighted least squares). The body is at its
engine position (z = -0.5..0.5) with the sensor at z = 0, as in a02.

Parts
  A. Rover-frame and world-frame DC offsets and drift fitted as nuisance parameters. Models: none, rover_dc
     (3 rover-frame offsets, rotated with heading), rover_world_dc (+3 world-frame offsets),
     rover_world_dc_drift (+3 world-frame linear drifts). Linearly dependent columns are dropped.
     A1  i.i.d. noise: precision lost by fitting the nuisance terms.
     A2  i.i.d. noise plus an unknown constant rover-frame offset (5 nT) and world-frame offset (2 nT).
     A3  correlated noise (a05 model: 0.05 nT white + OU, total 0.50 nT per axis), L = 0.5 cm to 100 m and
         L = inf (a constant offset per survey).
  B. Covariance scan: mean direction error over magnetization inclination and declination, sensor heights
     0.25-1.5 m and standoffs 0.5-2 m (i.i.d. noise, sampled from the estimator covariance).
  C. Buried extent: data from an outcrop extending b m below the surface, inverted with the exposed-only
     (0..1 m) model.
  D. Internal heterogeneity: the prism split into 4 x 4 x 2 cells (1 x 1 x 0.5 m), each magnetized in a
     Fisher-distributed direction about -z; the recovered direction is compared with the volume-mean vector
     and the Fisher mean direction.
  E. Along-track spacing 1-25 cm at fixed per-sample sigma and at fixed noise amplitude spectral density (ASD;
     sigma scaled by sqrt(1 cm / spacing)): accuracy is set by the ASD, not the spacing.
  F. Wavenumber band of the estimator weights, and a prediction of the correlated-noise errors from the ASD
     in that band.
  G. Regional field gradient (fitted or not) and slow drift about a constant (OU, L = 100 m), with rover- and
     local-frame offsets fitted.
  H. Local-frame offsets only fitted (direction scan), tolerance to an unfitted constant offset, and buried
     extent down to 5 m.
  I. Survey designs (approach leg, second loop) that separate a fitted site-field offset from the outcrop.
  J. Figure 8 panel B at 500 inversions per point (the check in the Figure 8 caption): the G2 drift
     requirement (drive-around, offsets fitted) rerun with n = 500, once with the random-number stream of G2
     (seed label a06:G2) and once with an independent one (a06:G2N500), compared with the n = 300 points of G2.

Run:  python analysis/a06_isolated_robustness.py [--parts ABCDEFGHIJ]
Out:  results/analysis/a06_isolated_robustness.json
"""
import os, sys, json, time, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "engine")); os.chdir(os.path.join(ROOT, "engine"))

import discretize
import isolated_outcrop_1cm as ISO
import correlated_noise_sensitivity as CN
from inversion_helpers import (prism_cell_operator, _sim_for, dir_err_deg, seed_for,
                               M_NRM_MIN, SIGMA_B, STANDOFF, PRISM_DIMS, PAPER_SURVEYS, PAPER_SPACING)
from core import generate_traverse

SPACING = PAPER_SPACING
SURVEYS = PAPER_SURVEYS
M_TRUE = M_NRM_MIN * np.array([0.0, 0.0, -1.0])
BND = tuple(float(v) for v in ISO.bnd)
L_GRID = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, np.inf]
NUIS = ["none", "rover_dc", "rover_world_dc", "rover_world_dc_drift", "rover_dc_world_poly3", "rover_dc_world_poly9"]
ALL_PARTS = "ABCDEFGHIJ"
OUT = os.path.join(ROOT, "results", "analysis", "a06_isolated_robustness.json")


def survey(name, height=0.5, standoff=STANDOFF, spacing=SPACING):
    kind, approach = SURVEYS[name]
    sp = {"cx": 0.0, "cy": 0.0, "sx": PRISM_DIMS[0], "sy": PRISM_DIMS[1], "sz": PRISM_DIMS[2]}
    return np.asarray(generate_traverse(traverse_type=kind, shape_params=sp, sensor_height=height - 0.5,
                                        spacing=spacing, orbit_offset=standoff, approach_distance=approach,
                                        shape="rect"), float)


def stats(a):
    a = np.asarray(a, float)
    return {"median": float(np.median(a)), "mean": float(np.mean(a)), "std": float(np.std(a)),
            "p95": float(np.percentile(a, 95)), "max": float(np.max(a)), "N": int(a.size)}


def metrics(Mrec, mt):
    th = np.array([dir_err_deg(m, mt) for m in Mrec])
    dm = 100 * (np.linalg.norm(Mrec, axis=1) - np.linalg.norm(mt)) / np.linalg.norm(mt)
    return {"dir_err_deg": stats(th), "dM_abs_pct": stats(np.abs(dm)), "dM_signed_pct": stats(dm)}


# ---------------------------------------------------------------- nuisance columns
def headings(sxyz):
    s = sxyz[:, :2]
    d = np.gradient(s, axis=0)
    return np.arctan2(d[:, 1], d[:, 0])


def nuisance_columns(sxyz, model):
    """Station-major (nD x q) columns for the DC offsets and drift of the given nuisance model."""
    n = len(sxyz)
    cols = []
    if model in ("rover_dc", "rover_world_dc", "rover_world_dc_drift"):
        psi = headings(sxyz)
        c, s = np.cos(psi), np.sin(psi)
        for v in (np.c_[c, s, 0 * c], np.c_[-s, c, 0 * c], np.c_[0 * c, 0 * c, 1 + 0 * c]):
            cols.append(v.reshape(-1))
    if model == "world_dc":
        for j in range(3):
            v = np.zeros((n, 3)); v[:, j] = 1.0; cols.append(v.reshape(-1))
    if model in ("rover_world_dc", "rover_world_dc_drift"):
        for j in range(3):
            v = np.zeros((n, 3)); v[:, j] = 1.0; cols.append(v.reshape(-1))
    if model.startswith("rover_dc_world_poly"):
        psi = headings(sxyz); c, s = np.cos(psi), np.sin(psi)
        for v in (np.c_[c, s, 0 * c], np.c_[-s, c, 0 * c]):
            cols.append(v.reshape(-1))
        t = np.linspace(-1.0, 1.0, n)
        for deg in range(int(model[len("rover_dc_world_poly"):]) + 1):
            P = np.polynomial.legendre.Legendre.basis(deg)(t)
            for j in range(3):
                v = np.zeros((n, 3)); v[:, j] = P; cols.append(v.reshape(-1))
    if model == "rover_world_dc_drift":
        t = np.linspace(-1.0, 1.0, n)
        for j in range(3):
            v = np.zeros((n, 3)); v[:, j] = t; cols.append(v.reshape(-1))
    if not cols:
        return np.zeros((3 * n, 0))
    return np.column_stack(cols)


def augmented(G, N):
    """Append the nuisance columns that are not linearly dependent on the columns already present."""
    A = G.copy()
    kept = 0
    for k in range(N.shape[1]):
        trial = np.column_stack([A, N[:, k]])
        if np.linalg.matrix_rank(trial, tol=1e-8 * np.linalg.norm(trial)) == trial.shape[1]:
            A = trial; kept += 1
    return A, kept


def lsK(A):
    return np.linalg.solve(A.T @ A, A.T)


# ---------------------------------------------------------------- Part A
def part_A(out):
    res = {}
    for name in SURVEYS:
        sxyz = survey(name)
        G = prism_cell_operator(BND, sxyz)
        d0 = G @ M_TRUE
        res[name] = {"n_stations": int(len(sxyz))}
        ops = {}
        for nm in NUIS:
            A, kept = augmented(G, nuisance_columns(sxyz, nm))
            ops[nm] = (A, lsK(A), kept)
            res[name].setdefault("n_nuisance_kept", {})[nm] = kept
        # A1: i.i.d. noise
        rng = np.random.default_rng(seed_for(f"a06:A1:{name}"))
        noise = rng.normal(0.0, SIGMA_B, size=(1000, d0.size))
        res[name]["A1_iid"] = {nm: metrics((d0[None] + noise) @ ops[nm][1][:3].T, M_TRUE) for nm in NUIS}
        # A2: i.i.d. noise plus unknown offsets (random directions, fixed amplitudes)
        rng = np.random.default_rng(seed_for(f"a06:A2:{name}"))
        psi = headings(sxyz); c, s = np.cos(psi), np.sin(psi)
        n2 = 500
        br = rng.normal(size=(n2, 3)); br *= 5.0 / np.linalg.norm(br, axis=1, keepdims=True)
        bw = rng.normal(size=(n2, 3)); bw *= 2.0 / np.linalg.norm(bw, axis=1, keepdims=True)
        off = np.empty((n2, len(sxyz), 3))
        off[:, :, 0] = br[:, [0]] * c - br[:, [1]] * s + bw[:, [0]]
        off[:, :, 1] = br[:, [0]] * s + br[:, [1]] * c + bw[:, [1]]
        off[:, :, 2] = br[:, [2]] + bw[:, [2]]
        D = d0[None] + off.reshape(n2, -1) + rng.normal(0.0, SIGMA_B, size=(n2, d0.size))
        res[name]["A2_offsets_5nT_rover_2nT_world"] = {nm: metrics(D @ ops[nm][1][:3].T, M_TRUE) for nm in NUIS}
        # A3: correlated noise
        res[name]["A3_correlated"] = []
        for L in L_GRID:
            rng = np.random.default_rng(seed_for(f"a06:A3:{name}:{L}"))
            nz = CN.make_noise(sxyz, L, 300, rng)
            row = {"L_m": ("inf" if np.isinf(L) else L)}
            for nm in NUIS:
                row[nm] = metrics((d0[None] + nz) @ ops[nm][1][:3].T, M_TRUE)
            res[name]["A3_correlated"].append(row)
            print(f"[A3 {name}] L={row['L_m']}: " + "  ".join(
                f"{nm} {row[nm]['dir_err_deg']['mean']:.1f}" for nm in NUIS), flush=True)
        for key in ("A1_iid", "A2_offsets_5nT_rover_2nT_world"):
            print(f"[{key} {name}] " + "  ".join(
                f"{nm} {res[name][key][nm]['dir_err_deg']['mean']:.2f}/{res[name][key][nm]['dM_abs_pct']['median']:.1f}%"
                for nm in NUIS), flush=True)
    out["A_dc_offsets"] = res


# ---------------------------------------------------------------- Part B
INC = [-90, -60, -30, 0, 30, 60, 90]
DEC = [0, 45, 90, 135, 180, 225, 270, 315]


def dir_scan(G, n=4000, seed=0, A=None):
    A = G if A is None else A
    C = SIGMA_B ** 2 * np.linalg.inv(A.T @ A)[:3, :3]
    Lc = np.linalg.cholesky(C)
    z = np.random.default_rng(seed).normal(size=(n, 3))
    rows = []
    for I in INC:
        for D in (DEC if abs(I) != 90 else [0]):
            Ir, Dr = np.radians(I), np.radians(D)
            u = np.array([np.cos(Ir) * np.cos(Dr), np.cos(Ir) * np.sin(Dr), -np.sin(Ir)])  # +I points down (-z)
            mt = M_NRM_MIN * u
            Mrec = mt[None] + z @ Lc.T
            th = np.degrees(np.arccos(np.clip(Mrec @ u / np.linalg.norm(Mrec, axis=1), -1, 1)))
            dm = 100 * np.abs(np.linalg.norm(Mrec, axis=1) - M_NRM_MIN) / M_NRM_MIN
            rows.append({"inc_deg": I, "dec_deg": D, "dir_err_mean": float(th.mean()),
                         "dir_err_median": float(np.median(th)), "dM_abs_median": float(np.median(dm))})
    m = np.array([r["dir_err_mean"] for r in rows])
    return {"rows": rows, "dir_err_mean_min": float(m.min()), "dir_err_mean_max": float(m.max()),
            "dir_err_mean_avg": float(m.mean()),
            "sigma_M_Am": np.sqrt(np.diag(C)).tolist()}


def part_B(out):
    res = {}
    for name in SURVEYS:
        res[name] = {"height": [], "standoff": []}
        for h in [0.25, 0.5, 0.75, 1.0, 1.5]:
            G = prism_cell_operator(BND, survey(name, height=h))
            r = dir_scan(G, seed=seed_for(f"a06:B:{name}:h{h}")); r["height_m"] = h
            res[name]["height"].append(r)
            print(f"[B {name}] h={h}: mean dir err over directions {r['dir_err_mean_min']:.2f}-{r['dir_err_mean_max']:.2f} deg", flush=True)
        sx = survey(name); G = prism_cell_operator(BND, sx)
        res[name]["with_nuisance"] = {}
        for nm in NUIS[1:]:
            A, kept = augmented(G, nuisance_columns(sx, nm))
            r = dir_scan(G, seed=seed_for(f"a06:B:{name}:nuis"), A=A); r["n_nuisance"] = kept
            res[name]["with_nuisance"][nm] = r
            print(f"[B {name}] nuisance {nm}: {r['dir_err_mean_min']:.2f}-{r['dir_err_mean_max']:.2f} deg; "
                  f"|dM| med {min(x['dM_abs_median'] for x in r['rows']):.1f}-{max(x['dM_abs_median'] for x in r['rows']):.1f}%", flush=True)
        for so in [0.5, 1.0, 2.0]:
            G = prism_cell_operator(BND, survey(name, standoff=so))
            r = dir_scan(G, seed=seed_for(f"a06:B:{name}:s{so}")); r["standoff_m"] = so
            res[name]["standoff"].append(r)
            print(f"[B {name}] standoff={so}: {r['dir_err_mean_min']:.2f}-{r['dir_err_mean_max']:.2f} deg", flush=True)
    out["B_direction_height_scan"] = res


# ---------------------------------------------------------------- Part C
def part_C(out):
    res = {}
    dirs = {"-z": np.array([0, 0, -1.0]), "+x": np.array([1.0, 0, 0]), "inc45": np.array([0.7071, 0, -0.7071])}
    for name in SURVEYS:
        sxyz = survey(name)
        Gi = prism_cell_operator(BND, sxyz)
        K = lsK(Gi)
        res[name] = []
        for b in [0.0, 0.1, 0.25, 0.5, 1.0]:
            bt = (BND[0], BND[1], BND[2], BND[3], BND[4] - b, BND[5])
            Gt = prism_cell_operator(bt, sxyz)
            row = {"buried_depth_m": b}
            for dn, u in dirs.items():
                mt = M_NRM_MIN * u
                rng = np.random.default_rng(seed_for(f"a06:C:{name}:{b}:{dn}"))
                D = (Gt @ mt)[None] + rng.normal(0.0, SIGMA_B, size=(1000, Gt.shape[0]))
                row[dn] = metrics(D @ K.T, mt)
                row[dn]["noise_free_dir_err_deg"] = dir_err_deg(K @ (Gt @ mt), mt)
                row[dn]["noise_free_dM_pct"] = float(100 * (np.linalg.norm(K @ (Gt @ mt)) / M_NRM_MIN - 1))
            res[name].append(row)
            print(f"[C {name}] b={b}: " + "  ".join(
                f"{dn} {row[dn]['dir_err_deg']['median']:.1f} deg (nf {row[dn]['noise_free_dir_err_deg']:.1f}, "
                f"dM {row[dn]['noise_free_dM_pct']:+.0f}%)" for dn in dirs), flush=True)
    out["C_buried_extent"] = res


# ---------------------------------------------------------------- Part D
def fisher(mean, kappa, n, rng):
    u = rng.uniform(size=n)
    w = 1 + np.log(u + (1 - u) * np.exp(-2 * kappa)) / kappa
    phi = rng.uniform(0, 2 * np.pi, n)
    v = np.c_[np.sqrt(1 - w ** 2) * np.cos(phi), np.sqrt(1 - w ** 2) * np.sin(phi), w]
    mean = mean / np.linalg.norm(mean)
    a = np.array([0, 0, 1.0])
    if np.allclose(mean, a):
        return v
    if np.allclose(mean, -a):
        return v * np.array([1, -1, -1])
    k = np.cross(a, mean); s = np.linalg.norm(k); c = a @ mean; k /= s
    Kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + s * Kx + (1 - c) * Kx @ Kx
    return v @ R.T


def part_D(out):
    x0, x1, y0, y1, z0, z1 = BND
    mesh = discretize.TensorMesh([np.full(4, 1.0), np.full(4, 1.0), np.full(2, 0.5)], origin=[x0, y0, z0])
    nC = mesh.n_cells
    res = {}
    for name in SURVEYS:
        sxyz = survey(name)
        Gc = _sim_for(mesh, sxyz, nC).G.copy()
        Gu = np.column_stack([Gc[:, 0:nC].sum(1), Gc[:, nC:2 * nC].sum(1), Gc[:, 2 * nC:].sum(1)])
        K = lsK(Gu)
        res[name] = []
        for kappa in [100.0, 30.0, 10.0, 5.0]:
            rng = np.random.default_rng(seed_for(f"a06:D:{name}:{kappa}"))
            th_mean, th_fish, dm_mean, th_nf = [], [], [], []
            for _ in range(500):
                U = fisher(np.array([0, 0, -1.0]), kappa, nC, rng)
                Mc = M_NRM_MIN * U
                m_full = np.concatenate([Mc[:, 0], Mc[:, 1], Mc[:, 2]])
                d = Gc @ m_full
                mrec_nf = K @ d
                mrec = K @ (d + rng.normal(0, SIGMA_B, d.size))
                vm = Mc.mean(0)
                th_mean.append(dir_err_deg(mrec, vm)); th_fish.append(dir_err_deg(mrec, [0, 0, -1]))
                th_nf.append(dir_err_deg(mrec_nf, vm))
                dm_mean.append(100 * (np.linalg.norm(mrec) / np.linalg.norm(vm) - 1))
            cell_scatter = float(np.degrees(np.arccos(np.clip(fisher(np.array([0, 0, -1.0]), kappa, 20000, rng)[:, 2] * -1, -1, 1))).mean())
            row = {"kappa": kappa, "mean_cell_angle_from_mean_deg": cell_scatter,
                   "err_vs_volume_mean": stats(th_mean), "err_vs_fisher_mean": stats(th_fish),
                   "noise_free_err_vs_volume_mean": stats(th_nf), "dM_vs_volume_mean_pct": stats(dm_mean)}
            res[name].append(row)
            print(f"[D {name}] kappa={kappa} (cells {cell_scatter:.0f} deg): vs vol-mean {row['err_vs_volume_mean']['mean']:.1f} "
                  f"(noise-free {row['noise_free_err_vs_volume_mean']['mean']:.1f}), vs Fisher mean {row['err_vs_fisher_mean']['mean']:.1f} deg", flush=True)
    out["D_heterogeneity"] = res


def part_E(out):
    """Recovery at different along-track spacings with the same noise ASD."""
    res = {}
    for name in SURVEYS:
        res[name] = []
        for sp in [0.01, 0.02, 0.05, 0.10, 0.25]:
            sx = survey(name, spacing=sp)
            G = prism_cell_operator(BND, sx)
            for label, sig in (("fixed_sigma", SIGMA_B), ("fixed_ASD", SIGMA_B * np.sqrt(0.01 / sp))):
                C = sig ** 2 * np.linalg.inv(G.T @ G)
                z = np.random.default_rng(seed_for(f"a06:E:{name}:{sp}:{label}")).normal(size=(4000, 3))
                Mrec = M_TRUE[None] + z @ np.linalg.cholesky(C).T
                r = metrics(Mrec, M_TRUE); r.update(spacing_m=sp, noise=label, sigma_nT=float(sig), n_stations=len(sx))
                res[name].append(r)
                print(f"[E {name}] spacing {sp*100:.0f} cm {label:11s} sigma {sig:.2f}: dir mean {r['dir_err_deg']['mean']:.2f} deg", flush=True)
    out["E_sampling_vs_ASD"] = res


def part_F(out):
    """Wavenumber sensitivity of the estimator, and a prediction of the OU-noise errors from the in-band ASD.
    For stationary noise with one-sided spatial PSD S(k) the estimator variance is sum_k S(k) |K_hat(k)|^2.
    Prediction: error(L) ~ error(white) * sqrt(S_OU(k_c; L) / S_white) at the band center k_c."""
    res = {}
    for name in SURVEYS:
        sx = survey(name)
        G = prism_cell_operator(BND, sx)
        K = lsK(G)                                     # (3, 3n)
        n = len(sx)
        W = K.reshape(3, n, 3)
        P = np.zeros(n // 2 + 1)
        for i in range(3):
            for j in range(3):
                P += np.abs(np.fft.rfft(W[i, :, j])) ** 2
        k = np.fft.rfftfreq(n, d=SPACING)              # cycles per metre
        c = np.cumsum(P) / P.sum()
        q = {f"k{int(p*100)}": float(k[np.searchsorted(c, p)]) for p in (0.05, 0.25, 0.5, 0.75, 0.95)}
        kc = q["k50"]
        white = 2 * SIGMA_B ** 2 * SPACING
        pred = []
        for L in [0.005, 0.01, 0.02, 0.05, 0.1]:
            S = 2 * SIGMA_B ** 2 * SPACING * 0 + 4 * CN.SIGMA_C ** 2 * L / (1 + (2 * np.pi * kc * L) ** 2) + 2 * CN.SIGMA_W ** 2 * SPACING
            pred.append({"L_m": L, "asd_ratio": float(np.sqrt(S / white))})
        res[name] = {"quantiles_cyc_per_m": q, "white_psd_nT2_per_cyc_per_m": white,
                     "white_asd_nT_per_rtHz_at_1cm_s": float(np.sqrt(white / 0.01)), "ou_prediction": pred}
        print(f"[F {name}] estimator band (cyc/m) 5/25/50/75/95%: " + ", ".join(f"{v:.2f}" for v in q.values())
              + f"; white ASD at 1 cm/s = {np.sqrt(white/0.01):.2f} nT/rtHz; OU ASD ratios "
              + ", ".join(f"L={p['L_m']}: {p['asd_ratio']:.2f}" for p in pred), flush=True)
    out["F_estimator_band"] = res


def gradient_basis(sxyz):
    """Five station-major columns for a uniform field gradient (symmetric traceless 3x3 tensor T,
    field b = T r about the survey centroid)."""
    r = sxyz - sxyz.mean(0)
    E = []
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        T = np.zeros((3, 3)); T[i, j] = T[j, i] = 1.0 / np.sqrt(2); E.append(T)
    T = np.diag([1.0, -1.0, 0.0]) / np.sqrt(2); E.append(T)
    T = np.diag([1.0, 1.0, -2.0]) / np.sqrt(6); E.append(T)
    return np.column_stack([(r @ T.T).reshape(-1) for T in E]), E


DRIFT_SD_NT = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5]


def drift_rows(name, sxyz, d0, K_off, n, label, tag):
    """Drift requirement (parts G2 and J): white 0.5 nT plus an OU drift (L = 100 m) of per-axis std sd, inverted
    with rover- and local-frame offsets fitted (K_off), n inversions per point. The seed label is
    f"a06:{label}:{name}:{sd}"."""
    drows = []
    for sd in DRIFT_SD_NT:
        rng = np.random.default_rng(seed_for(f"a06:{label}:{name}:{sd}"))
        nz = rng.normal(0.0, SIGMA_B, size=(n, d0.size))
        if sd > 0:
            ou = CN.make_noise(sxyz, 100.0, n, rng)          # total std SIGMA_TOT, mostly OU
            nz = nz + ou * (sd / CN.SIGMA_TOT)
        M = (d0[None] + nz) @ K_off[:3].T
        r = {"drift_std_nT": sd, "offsets_fitted": metrics(M, M_TRUE)}
        # rms of the drift about its per-survey mean
        if sd > 0:
            v = (ou * (sd / CN.SIGMA_TOT)).reshape(n, -1, 3)
            r["rms_about_mean_nT"] = float(np.sqrt(np.mean((v - v.mean(1, keepdims=True)) ** 2)))
        drows.append(r)
        print(f"[{tag} {name}] drift std {sd}: {r['offsets_fitted']['dir_err_deg']['mean']:.1f} deg"
              + (f" (rms about mean {r['rms_about_mean_nT']:.3f} nT)" if sd > 0 else ""), flush=True)
    return drows


def part_G(out):
    """G1: a regional field gradient (uniform tensor, Frobenius norm g in nT/m) added to the data, inverted
    with offsets fitted, with and without the five gradient terms fitted.
    G2: slow drift about a constant (OU, L = 100 m, per-axis std s) added to the white 0.5 nT, inverted with
    offsets fitted; s that gives a mean direction error of 5 deg."""
    res = {}
    for name in SURVEYS:
        sxyz = survey(name)
        G = prism_cell_operator(BND, sxyz)
        d0 = G @ M_TRUE
        N0 = nuisance_columns(sxyz, "rover_world_dc")
        A_off, _ = augmented(G, N0)
        Bg, E = gradient_basis(sxyz)
        A_grad, kept_g = augmented(A_off, Bg)
        K_off, K_grad = lsK(A_off), lsK(A_grad)
        rows = []
        for g in [0.0, 0.003, 0.01, 0.03, 0.1, 0.3]:
            rng = np.random.default_rng(seed_for(f"a06:G1:{name}:{g}"))
            n = 500
            coef = rng.normal(size=(n, 5)); coef *= g / np.linalg.norm(coef, axis=1, keepdims=True)
            D = d0[None] + coef @ Bg.T + rng.normal(0.0, SIGMA_B, size=(n, d0.size))
            r = {"gradient_nT_per_m": g,
                 "offsets_fitted": metrics(D @ K_off[:3].T, M_TRUE),
                 "offsets_and_gradient_fitted": metrics(D @ K_grad[:3].T, M_TRUE)}
            rows.append(r)
            print(f"[G1 {name}] g={g}: offsets only {r['offsets_fitted']['dir_err_deg']['mean']:.1f} deg; "
                  f"+gradient {r['offsets_and_gradient_fitted']['dir_err_deg']['mean']:.1f} deg "
                  f"(|dM| {r['offsets_and_gradient_fitted']['dM_abs_pct']['median']:.1f}%)", flush=True)
        # G2: drift requirement
        drows = drift_rows(name, sxyz, d0, K_off, n=300, label="G2", tag="G2")
        res[name] = {"n_gradient_terms_kept": kept_g, "gradient": rows, "drift_L100m": drows}
    out["G_gradient_and_drift"] = res


def part_H(out):
    """H1: only local-frame (world) offsets fitted, scanned over 42 directions. The static site field cannot be
    removed using the lander and must be fitted; a rover-frame offset can be calibrated by rotating the rover.
    H2: tolerance to an unfitted constant offset (local frame, random direction) of amplitude c.
    H3: deeper buried extent (b up to 5 m), three magnetization directions, drive-around only."""
    res = {}
    for name in SURVEYS:
        sxyz = survey(name)
        G = prism_cell_operator(BND, sxyz)
        A, kept = augmented(G, nuisance_columns(sxyz, "world_dc"))
        r = dir_scan(G, seed=seed_for(f"a06:H1:{name}"), A=A); r["n_nuisance"] = kept
        vert = [x for x in r["rows"] if abs(x["inc_deg"]) == 90][0]
        print(f"[H1 {name}] local offsets fitted: {r['dir_err_mean_min']:.2f}-{r['dir_err_mean_max']:.2f} deg; "
              f"|dM| med {min(x['dM_abs_median'] for x in r['rows']):.1f}-{max(x['dM_abs_median'] for x in r['rows']):.1f}%; "
              f"vertical {vert['dir_err_mean']:.2f} deg / {vert['dM_abs_median']:.1f}%", flush=True)
        d0 = G @ M_TRUE; K = lsK(G)
        rows = []
        for c in [0.0, 0.003, 0.01, 0.03, 0.1, 0.3]:
            rng = np.random.default_rng(seed_for(f"a06:H2:{name}:{c}"))
            n = 1000
            b = rng.normal(size=(n, 3)); b *= c / np.linalg.norm(b, axis=1, keepdims=True)
            off = np.repeat(b[:, None, :], len(sxyz), axis=1).reshape(n, -1)
            D = d0[None] + off + rng.normal(0.0, SIGMA_B, size=(n, d0.size))
            m = metrics(D @ K.T, M_TRUE); m["offset_nT"] = c
            rows.append(m)
            print(f"[H2 {name}] unfitted constant offset {c} nT: {m['dir_err_deg']['mean']:.1f} deg, |dM| {m['dM_abs_pct']['median']:.1f}%", flush=True)
        res[name] = {"H1_world_offsets_fitted": r, "H2_unfitted_offset": rows}
    # H3: deeper burial (drive-around)
    sxyz = survey("drive_around")
    Gi = prism_cell_operator(BND, sxyz); K = lsK(Gi)
    dirs = {"-z": np.array([0, 0, -1.0]), "+x": np.array([1.0, 0, 0]), "inc45": np.array([0.7071, 0, -0.7071])}
    rows = []
    for b in [1.0, 2.0, 3.0, 5.0]:
        bt = (BND[0], BND[1], BND[2], BND[3], BND[4] - b, BND[5])
        Gt = prism_cell_operator(bt, sxyz)
        row = {"buried_depth_m": b}
        for dn, u in dirs.items():
            mt = M_NRM_MIN * u
            mnf = K @ (Gt @ mt)
            row[dn] = {"noise_free_dir_err_deg": dir_err_deg(mnf, mt), "noise_free_dM_pct": float(100 * (np.linalg.norm(mnf) / M_NRM_MIN - 1))}
        rows.append(row)
        print(f"[H3] b={b}: " + "  ".join(f"{dn} {row[dn]['noise_free_dir_err_deg']:.1f} deg / {row[dn]['noise_free_dM_pct']:+.0f}%" for dn in dirs), flush=True)
    res["H3_deep_burial_drive_around"] = rows
    out["H_offsets_and_burial"] = res


def part_I(out):
    """Survey designs that separate a constant site field (a local-frame offset, which must be fitted because the
    lander cannot remove the static crustal field at the rover) from the outcrop field. Drive-around loop at 1 m
    alone, with a 10 m straight approach leg, or with a second loop at 3 m standoff. Direction scan over 42
    directions with local-frame offsets fitted, and with offsets plus the five gradient terms fitted."""
    loop1 = survey("drive_around")
    loop3 = survey("drive_around", standoff=3.0)
    z0 = loop1[0, 2]
    leg = np.c_[np.linspace(-13.0, -3.0 - SPACING, 1000), np.full(1000, -3.0), np.full(1000, z0)]
    designs = {"loop_1m": loop1, "loop_1m+10m_leg": np.vstack([leg, loop1]),
               "loop_1m+loop_3m": np.vstack([loop1, loop3])}
    res = {}
    for dn, sx in designs.items():
        G = prism_cell_operator(BND, sx)
        A1, k1 = augmented(G, nuisance_columns(sx, "world_dc"))
        Bg, _ = gradient_basis(sx)
        A2, k2 = augmented(A1, Bg)
        r1 = dir_scan(G, seed=seed_for(f"a06:I:{dn}:dc"), A=A1)
        r2 = dir_scan(G, seed=seed_for(f"a06:I:{dn}:dcgrad"), A=A2)
        r0 = dir_scan(G, seed=seed_for(f"a06:I:{dn}:none"))
        rng_dm = lambda r: (min(x["dM_abs_median"] for x in r["rows"]), max(x["dM_abs_median"] for x in r["rows"]))
        res[dn] = {"n_stations": int(len(sx)), "no_offsets": r0, "local_offsets": r1, "local_offsets_and_gradient": r2,
                   "n_kept": [k1, k2]}
        print(f"[I {dn}] n={len(sx)} | none {r0['dir_err_mean_min']:.1f}-{r0['dir_err_mean_max']:.1f} | "
              f"offsets {r1['dir_err_mean_min']:.1f}-{r1['dir_err_mean_max']:.1f} deg, |dM| {rng_dm(r1)[0]:.1f}-{rng_dm(r1)[1]:.1f}% | "
              f"+gradient {r2['dir_err_mean_min']:.1f}-{r2['dir_err_mean_max']:.1f} deg, |dM| {rng_dm(r2)[0]:.1f}-{rng_dm(r2)[1]:.1f}%", flush=True)
    out["I_survey_designs_for_offsets"] = res


def part_J(out):
    """Figure 8 panel B at 500 inversions per point (drive-around, offsets fitted), the check quoted in the
    figure caption. Survey, operator and noise model as in G2 (n = 300, stored in G_gradient_and_drift), with
    n = 500 and either the G2 seed label (same_seed) or the independent label a06:G2N500 (independent_seed).
    The result is the mean direction error at each drift level, to compare with the stored n = 300 points."""
    name = "drive_around"
    sxyz = survey(name)
    G = prism_cell_operator(BND, sxyz)
    d0 = G @ M_TRUE
    A_off, _ = augmented(G, nuisance_columns(sxyz, "rover_world_dc"))
    K_off = lsK(A_off)
    if "G_gradient_and_drift" in out:
        n300 = out["G_gradient_and_drift"][name]["drift_L100m"]
    else:
        n300 = json.load(open(OUT))["G_gradient_and_drift"][name]["drift_L100m"]
    runs = {"same_seed": drift_rows(name, sxyz, d0, K_off, n=500, label="G2", tag="J same_seed"),
            "independent_seed": drift_rows(name, sxyz, d0, K_off, n=500, label="G2N500", tag="J independent_seed")}
    cmp = {"drift_std_nT": DRIFT_SD_NT,
           "rms_about_mean_nT": [r.get("rms_about_mean_nT", 0.0) for r in runs["independent_seed"]],
           "N300_mean_dir_err_deg": [r["offsets_fitted"]["dir_err_deg"]["mean"] for r in n300]}
    for k, rows in runs.items():
        cmp[f"N500_{k}_mean_dir_err_deg"] = [r["offsets_fitted"]["dir_err_deg"]["mean"] for r in rows]
        cmp[f"{k}_diff_deg"] = [a - b for a, b in zip(cmp[f"N500_{k}_mean_dir_err_deg"], cmp["N300_mean_dir_err_deg"])]
        cmp[f"{k}_max_abs_diff_deg"] = max(abs(v) for v in cmp[f"{k}_diff_deg"])
    out["J_fig8_n500_check"] = {"N": 500, "survey": name,
                                "seed_labels": {"same_seed": "a06:G2:<survey>:<sd>", "independent_seed": "a06:G2N500:<survey>:<sd>"},
                                "runs": runs, "comparison_with_N300": cmp}
    print(f"[J] max |mean direction error difference| vs N = 300: same seed {cmp['same_seed_max_abs_diff_deg']:.2f} deg, "
          f"independent seed {cmp['independent_seed_max_abs_diff_deg']:.2f} deg", flush=True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--parts", default=ALL_PARTS); a = ap.parse_args()
    t0 = time.time()
    # With --parts, keep the stored results of the other parts
    out = json.load(open(OUT)) if a.parts != ALL_PARTS and os.path.exists(OUT) else {}
    out["meta"] = dict(script="analysis/a06_isolated_robustness.py", sigma_B_nT=SIGMA_B, spacing_m=SPACING,
                       M_true_Am=M_TRUE.tolist(), geometry="as a02", L_grid=[("inf" if np.isinf(L) else L) for L in L_GRID])
    for p in a.parts:
        {"A": part_A, "B": part_B, "C": part_C, "D": part_D, "E": part_E, "F": part_F, "G": part_G, "H": part_H, "I": part_I, "J": part_J}[p](out)
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"part {p} done ({time.time()-t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
