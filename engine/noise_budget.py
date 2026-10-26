"""Table 1 of the paper (per-axis magnetic field uncertainty budget) as data, and the noise level used in the simulations.

Narrett, Zareck, Weiss & Hoffman, Table 1 (Section 3.2.4, "Noise Budget"): 1-sigma contributions per axis for
the rover and the lander, for 1 s averages (one measurement per centimeter at 1 cm/s). The rover total is the
noise sigma_B in every forward model and inversion; the correlated-noise tests split the same 0.50 nT into a
0.05 nT white instrument term and a correlated remainder.

The noise level is imported from here by:
    inversion_helpers.SIGMA_B                = ROVER_ADOPTED_NT      (0.50 nT)
    skatepark_1cm_dr_sweep.noise_floor       = ROVER_ADOPTED_NT
    correlated_noise_sensitivity.SIGMA_TOT   = ROVER_ADOPTED_NT
    correlated_noise_sensitivity.SIGMA_W     = ROVER_INSTRUMENT_NT   (0.05 nT)
and notebook/paper_numbers.py checks the printed table and its totals against ROWS.

Each row of ROWS gives the value in nT per axis for the rover and the lander (None if the term is not in that
column), whether it counts toward the total, what it is, and its basis: either a published reference or the
design allocation adopted for this study. Values are the table's rounded entries. validate(), which runs on
import, checks that the quadrature sum rounded to two decimals plus the margin equals the adopted total.
"""
import math

UNITS = "nT per Cartesian axis, 1-sigma, for 1 s (one-centimeter) averages"

# Instrument noise before averaging (Section 3.2.4)
SAMPLE_NOISE_PER_SAMPLE_NT = 0.35     # nT per axis per sample at 50 Hz
SAMPLES_PER_AVERAGE = 50              # 50 Hz at 1 cm/s

_ALLOC = "design allocation adopted for this study"

ROWS = [
    {
        "key": "instrument_noise",
        "source": "Instrument noise (1 s average)",
        "rover_nT": 0.05, "lander_nT": 0.05,
        "in_total": True,
        "what": "Magnetometer noise after averaging the 50 samples collected per centimeter.",
        "basis": ("Per-sample noise of 0.35 nT per axis divided by sqrt(50) gives about 0.05 nT per axis "
                  "(Section 3.2.4); 50 Hz sampling and 1 cm/s drive speed from Clark (2025)."),
        "reference": "Clark2025 (sampling rate and rover speed)",
    },
    {
        "key": "imf_raw",
        "source": "IMF (raw; not in total)",
        "rover_nT": 4.0, "lander_nT": 4.0,
        "in_total": False,
        "what": "Interplanetary magnetic field at 1 AU before lander subtraction; shown for scale, not part of the total.",
        "basis": ("Typical IMF amplitude at 1 AU of about 4 nT (1-sigma per axis), varying over seconds to hours "
                  "(Section 3.2.2)."),
        "reference": "ness_interaction_1971 (Ness 1971); TangSW (Tang et al. 2025)",
    },
    {
        "key": "imf_residual",
        "source": "IMF (after lander subtraction)",
        "rover_nT": 0.18, "lander_nT": None,
        "in_total": True,
        "what": "IMF left on the rover after subtracting the static lander magnetometer's measurement.",
        "basis": (_ALLOC.capitalize() + ": the IMF is reduced from about 4 nT to about 0.18 nT per axis, limited by "
                  "spatial uncertainties between the lander and rover positions over baselines of 1 km or less "
                  "(Section 3.2.2). Equal to the lander's adopted total, so lander error enters the rover budget once."),
        "reference": None,
    },
    {
        "key": "stray_field",
        "source": "Spacecraft stray field (AC + DC)",
        "rover_nT": 0.40, "lander_nT": 0.05,
        "in_total": True,
        "what": "Residual uncertainty in the rover's (and lander's) own magnetic field at the sensor.",
        "basis": ("Field knowledge of better than 0.4 nT per axis on the rover and 0.05 nT per axis on the lander, "
                  "consistent with heritage boom-mounted magnetometers (Section 3.2.3); values from the "
                  + _ALLOC + "."),
        "reference": "Weiss2023 (Weiss et al. 2023, Space Sci. Rev., Psyche magnetometry investigation); "
                     "Ream2021 (Ream et al. 2022, Meas. Sci. Technol.)",
    },
    {
        "key": "lander_reconstruction",
        "source": "Lander reconstruction",
        "rover_nT": None, "lander_nT": 0.16,
        "in_total": True,
        "what": "Error in transferring the lander's estimate of the external field to the rover's position and time.",
        "basis": (_ALLOC.capitalize() + ": 0.16 nT per axis total reconstructed lander uncertainty (Table 1 notes). "
                  "With the instrument and stray-field terms this gives 0.17 nT, adopted as 0.18 nT with margin."),
        "reference": None,
    },
]

# Rover total of the design allocation: sqrt(0.45^2 + 0.18^2) + 0.02 nT margin = 0.50 nT.
# Table 1 sums 0.05, 0.18 and 0.40 nT in quadrature (0.44 nT) and adds a 0.06 nT margin.
ALLOC_ROVER_RECONSTRUCTED_NT = 0.45
ALLOC_IMF_RESIDUAL_NT = 0.18
ALLOC_MARGIN_NT = 0.02
ALLOC_ROVER_TOTAL_NT = round(math.sqrt(ALLOC_ROVER_RECONSTRUCTED_NT ** 2 + ALLOC_IMF_RESIDUAL_NT ** 2)
                             + ALLOC_MARGIN_NT, 2)      # 0.50

# Margin added to the rounded quadrature sum, and the adopted totals (Table 1, last lines)
ROVER_MARGIN_NT = 0.06
LANDER_MARGIN_NT = 0.01
ROVER_ADOPTED_NT = 0.50       # noise sigma_B in all forward models and inversions
LANDER_ADOPTED_NT = 0.18

ROVER_INSTRUMENT_NT = ROWS[0]["rover_nT"]      # white instrument term for the correlated-noise tests


def rss(column):
    """Quadrature sum of the rows that count toward the total; column is 'rover_nT' or 'lander_nT'."""
    return math.sqrt(sum(r[column] ** 2 for r in ROWS if r["in_total"] and r[column] is not None))


def totals():
    """Totals of Table 1 computed from ROWS: {'rover': {...}, 'lander': {...}}, each with the quadrature sum
    (unrounded and rounded to two decimals), the margin and the adopted total."""
    out = {}
    for who, col, margin, adopted in (("rover", "rover_nT", ROVER_MARGIN_NT, ROVER_ADOPTED_NT),
                                      ("lander", "lander_nT", LANDER_MARGIN_NT, LANDER_ADOPTED_NT)):
        s = rss(col)
        out[who] = {"rss_nT": s, "rss_rounded_nT": round(s, 2), "margin_nT": margin, "adopted_nT": adopted}
    return out


def validate():
    """Check that the table is self-consistent; raises AssertionError otherwise."""
    t = totals()
    for who, v in t.items():
        assert math.isclose(v["rss_rounded_nT"] + v["margin_nT"], v["adopted_nT"], abs_tol=1e-9), \
            f"{who}: rounded quadrature sum {v['rss_rounded_nT']} + margin {v['margin_nT']} != adopted {v['adopted_nT']}"
    assert round(SAMPLE_NOISE_PER_SAMPLE_NT / math.sqrt(SAMPLES_PER_AVERAGE), 2) == ROWS[0]["rover_nT"]
    assert ALLOC_ROVER_TOTAL_NT == ROVER_ADOPTED_NT
    imf = next(r for r in ROWS if r["key"] == "imf_residual")
    assert imf["rover_nT"] == LANDER_ADOPTED_NT, "rover IMF residual should equal the lander's adopted total"
    return t


validate()
