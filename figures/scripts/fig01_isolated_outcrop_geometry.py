"""Figure 1 - Isolated outcrop geometry and survey paths.

The geometry comes from engine/inversion_helpers.py (PRISM_DIMS, SENSOR_H, STANDOFF, PAPER_SURVEYS, prism_bounds,
paper_survey); only the prism and the survey points are built. Prism 4 x 4 x 1 m on the ground (x, y in [-2, 2],
z in [0, 1]), M = 2.4e-3 A/m along -z; sensor at z = 0.5 m. Drive-around: closed loop 1 m outside the faces
(x, y = +-3 m) starting and ending at (-3, -3). Drive-by: straight pass along y = -3 m from x = -8 to +8 m.
(A) plan view (xy); (B) vertical section through y = 0 (xz), to scale, under (A) on the same x axis.
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import paperstyle as ps

sys.path.insert(0, os.path.join(ps.ROOT, "engine"))
import inversion_helpers as IH  # simulation definitions

# Prism on the ground (z in [0, sz]) and the two surveys
PX0, PX1, PY0, PY1, PZ0, PZ1 = IH.prism_bounds((0.0, 0.0, IH.PRISM_DIMS[2] / 2), IH.PRISM_DIMS)
HALF, H = (PX1 - PX0) / 2, PZ1 - PZ0
LOOP_XYZ = IH.paper_survey("drive_around")   # closed loop, starts and ends at the same corner
BY_XYZ = IH.paper_survey("drive_by")
XL, YA, ZB = (-8.8, 8.8), (-7.0, 7.0), (-1.1, 2.1)   # axis limits


def main():
    hA, hB = ps.h_true(XL, YA), ps.h_true(XL, ZB)
    T = ps.GRID_TOP
    bottom = T + hA + ps.SHARE_GAP + hB
    g = ps.GeoGrid(ps.fig_height(bottom)); fig = g.fig
    x0 = (ps.GEO_W - ps.COL_W) / 2          # single column, centered under the colorbar
    axA = g.axes(x0, T, ps.COL_W, hA)
    axB = g.axes(x0, T + hA + ps.SHARE_GAP, ps.COL_W, hB, sharex=axA)

    # (A) plan view
    axA.set_facecolor(ps.BG)
    axA.add_patch(Rectangle((PX0, PY0), PX1 - PX0, PY1 - PY0, fc=ps.SRC, ec="black", lw=2, zorder=3))
    for xm in (-HALF / 2, HALF / 2):
        for ym in (-HALF / 2, HALF / 2):
            ps.mark_into(axA, xm, ym, s=190)
    lp = np.vstack([LOOP_XYZ, LOOP_XYZ[:1]])   # close the loop back to its start
    ps.trav(axA, lp[:, 0], lp[:, 1], ps.LS_LOOP, "Drive-around")
    ps.trav(axA, BY_XYZ[:, 0], BY_XYZ[:, 1], ps.LS_BY, "Drive-by", zorder=5)
    ps.start_dot(axA, LOOP_XYZ[0, 0], LOOP_XYZ[0, 1])
    ps.start_legend_handle(axA)
    ps.triad(axA, -7.6, -5.6, 1.5, ("x", "y"), color="white", out_of_page="z")
    axA.set_xlim(*XL); axA.set_ylim(*YA); axA.set_aspect("equal")
    axA.set_yticks([-5, 0, 5]); axA.set_ylabel(r"$y$ (m)")
    plt.setp(axA.get_xticklabels(), visible=False)
    ps.panel_letter(axA, "A", dark=True, size=24)

    # (B) section through y = 0, to scale
    ps.ground(axB, -12, 12, 5)
    axB.add_patch(Rectangle((PX0, PZ0), PX1 - PX0, H, fc=ps.SRC, ec="black", lw=2, zorder=3))
    for xa in (-HALF / 2, HALF / 2):
        ps.arrow(axB, xa, PZ0 + 0.88 * H, 0, -0.76 * H, lw=1.5, ms=10)
    ps.guide(axB, [-12, 12], float(LOOP_XYZ[0, 2]))
    axB.set_ylim(*ZB); axB.set_aspect("equal")
    axB.set_xticks([-5, 0, 5]); axB.set_yticks([0, 1, 2])
    axB.set_xlabel(r"$x$ (m)"); axB.set_ylabel(r"$z$ (m)")
    ps.panel_letter(axB, "B", dark=False, size=24, y=0.97)

    hA_, lA_ = axA.get_legend_handles_labels()
    mx = axB.scatter([], [], marker="x", s=120, c="black", linewidths=2.0)
    ps.bottom_legend(g, hA_ + [mx], lA_ + [r"$\mathbf{M}$ into page"], bottom)
    ps.colorbar_mid(g)
    ps.save_fixed(fig, "fig01_isolated_outcrop_geometry")


if __name__ == "__main__":
    main()
