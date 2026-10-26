"""Figure 3 - Three-layer wall geometry, from engine/three_layer_wall.py (LAYER_Z, M_LAYERS, BASE_Y, BASE_Z,
ELEVATED_Z, traverse()) and engine/inversion_helpers.three_layer_mesh (wall extents). The wall:

  flows: x in [-50, 50] m, y in [0, 50] m (into the wall), z in [0,5], [5,10], [10,15] m
  M: flow 1 (bottom) -z, flow 2 +y (into the wall), flow 3 (top) +z, |M| = 2.4e-3 A/m
  ground traverse: y = -1 m (1 m standoff from the wall face), z = 0.5 m, x = -120..+120 m
  elevated traverse (Table 4): same line at z = 12.5 m; the height sweep is Figure 7

(A) plan view (xy), to scale: top of the wall (flow 3) and the traverse line in front of the face
(B) the wall face seen from the rover side (xz), vertical exaggeration 3
(C) vertical section through x = 0 (yz), to scale
"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import paperstyle as ps

sys.path.insert(0, os.path.join(ps.ROOT, "engine"))
import three_layer_wall as TL   # scenario definition

# Wall extents from the mesh (x in [-XW, XW], y in [0, YW], z cells = the three flows)
_mesh = TL.H.three_layer_mesh()
XW = float(_mesh.nodes_x[-1])
YW = float(_mesh.nodes_y[-1])
assert np.allclose(_mesh.nodes_x[[0, -1]], [-XW, XW]) and _mesh.nodes_y[0] == 0.0
assert np.allclose(_mesh.nodes_z, [0.0] + [zt for _, zt in TL.LAYER_Z]), "LAYER_Z disagrees with three_layer_mesh"
# Flows: (z bottom, z top, label, unit magnetization direction)
LAYERS = [(zb, zt, f"Flow {i + 1}", tuple(int(round(v)) for v in TL.M_LAYERS[i] / TL.M_NRM_MIN))
          for i, (zb, zt) in enumerate(TL.LAYER_Z)]
# Traverse stations (1 cm spacing)
TRAV_LO = TL.traverse(TL.BASE_Y, TL.BASE_Z)        # ground traverse, y = -1 m, z = 0.5 m
TRAV_HI = TL.traverse(TL.BASE_Y, TL.ELEVATED_Z)    # elevated traverse (Table 4), z = 12.5 m
VE = 3.0      # vertical exaggeration of the wall face (B); (A) and (C) are to scale
YA, ZL, YC = (-10.0, 60.0), (-2.5, 17.5), (-8.0, 62.0)   # axis limits
XA = 75.0     # half-width of the plan view (A); the traverse continues to +-120 m


def panel_plan(ax):
    ax.set_facecolor(ps.BG)
    ax.add_patch(Rectangle((-XW, 0), 2 * XW, YW, fc=ps.SRC, ec="black", lw=2, zorder=3))
    for xa in (-35, -12, 12, 35):
        for ya in (15, 35):
            ps.mark_out(ax, xa, ya, s=110, lw=1.8)
    ps.trav(ax, TRAV_LO[:, 0], TRAV_LO[:, 1], ps.LS_LOOP, zorder=6)
    ps.triad(ax, -68, 16, 11, ("x", "y"), color="white", out_of_page="z")
    ax.set_xlim(-XA, XA); ax.set_ylim(*YA); ax.set_aspect("equal")
    ax.set_xticks([-50, 0, 50]); ax.set_yticks([0, 25, 50])
    ax.set_ylabel(r"$y$ (m)"); plt.setp(ax.get_xticklabels(), visible=False)
    ps.panel_letter(ax, "A", dark=True, size=24)


def panel_face(ax):
    ps.ground(ax, -140, 140, 12)
    for z0, z1, name, m in LAYERS:
        ax.add_patch(Rectangle((-XW, z0), 2 * XW, z1 - z0, fc=ps.SRC, ec="black", lw=2, zorder=3))
        zc = 0.5 * (z0 + z1)
        for xa in (-35, -12, 12, 35):
            if m[1] == 0:
                ps.arrow(ax, xa, zc - m[2] * 1.6, 0, m[2] * 3.2)
            else:
                ps.mark_into(ax, xa, zc, s=80, lw=1.8)
    lo, hi = TRAV_LO[np.abs(TRAV_LO[:, 0]) <= XA + 1e-6], TRAV_HI[np.abs(TRAV_HI[:, 0]) <= XA + 1e-6]
    # Only stations inside the plotted window are drawn, so the dash pattern starts at the frame edge
    ps.trav(ax, lo[:, 0], lo[:, 2], ps.LS_LOOP, "Ground traverse", zorder=6)
    # Dashes as in ps.LS_ELEV, with a phase offset of 7 that places the pattern as if it started just outside the frame
    ps.trav(ax, hi[:, 0], hi[:, 2], (7.0, (5, 3)), zorder=6)
    ax.plot([], [], color=ps.TRAV_DARK, lw=ps.GEO_TRAV_LW, ls=ps.LS_ELEV, label="Elevated traverse",
            solid_capstyle="butt", dash_capstyle="butt")   # legend handle, default dash phase
    ax.set_ylim(*ZL); ax.set_aspect(VE)
    ax.set_xticks([-50, 0, 50]); ax.set_yticks([0, 5, 10, 15])
    ax.set_xlabel(r"$x$ (m, along strike)"); ax.set_ylabel(r"$z$ (m)")
    ps.panel_letter(ax, "B", dark=False, size=24, y=0.97)


def panel_section(ax):
    ps.ground(ax, -14, 92, 18)
    for z0, z1, name, m in LAYERS:
        ax.add_patch(Rectangle((0, z0), YW, z1 - z0, fc=ps.SRC, ec="black", lw=2, zorder=3))
        zc = 0.5 * (z0 + z1)
        for ya in (12, 25, 38):
            if m[1] == 0:
                ps.arrow(ax, ya, zc - m[2] * 1.7, 0, m[2] * 3.4, lw=1.5, ms=10)
            else:
                ps.arrow(ax, ya - 4, zc, 8, 0, lw=1.5, ms=10)
        ax.text(YW + 1.2, zc, name, ha="left", va="center", fontsize=11.5, zorder=8)
    ax.set_xlim(*YC); ax.set_ylim(*ZL); ax.set_aspect("equal")
    ax.set_xticks([0, 25, 50]); ax.set_yticks([0, 5, 10, 15])
    ax.set_xlabel(r"$y$ (m, into the wall)"); ax.set_ylabel(r"$z$ (m)")
    ps.panel_letter(ax, "C", dark=False, size=24, y=0.97)


def main():
    T = ps.GRID_TOP; XLIM = (-XA, XA)
    hA, hB, hC = ps.h_true(XLIM, YA), ps.h_true(XLIM, ZL, VE), ps.h_true(YC, ZL)
    bottom = T + hA + ps.SHARE_GAP + hB
    g = ps.GeoGrid(ps.fig_height(bottom)); fig = g.fig
    x0, x1 = ps.COL_X
    axA = g.axes(x0, T, ps.COL_W, hA)
    axB = g.axes(x0, T + hA + ps.SHARE_GAP, ps.COL_W, hB, sharex=axA)
    axC = g.axes(x1, T + (bottom - T - hC) / 2 - 0.2, ps.COL_W, hC)   # centered beside A and B (space for labels below)
    panel_plan(axA); panel_face(axB); panel_section(axC)
    hL, lL = axB.get_legend_handles_labels()
    mx = axC.scatter([], [], marker="x", s=120, c="black", linewidths=2.0)
    mo = axC.scatter([], [], marker="o", s=120, facecolors="none", edgecolors="black", linewidths=2.0)
    ps.bottom_legend(g, hL + [mx, mo], lL + [r"$\mathbf{M}$ into page", r"$\mathbf{M}$ out of page"], bottom)
    ps.colorbar_mid(g)
    ps.save_fixed(fig, "fig03_three_layer_geometry")


if __name__ == "__main__":
    main()
