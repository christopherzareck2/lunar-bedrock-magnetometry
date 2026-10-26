"""Figure 8 - How correlated noise and slow drift degrade the isolated-outcrop recovery (drive-around).

Data: results/analysis/a06_isolated_robustness.json
  (A) part A3: mean direction error vs OU correlation length L (total 0.50 nT per axis), inverted without
      nuisance terms and with the five rover- and local-frame offsets fitted. L = inf is a constant offset.
      Dotted: prediction for short L, white-noise error x sqrt(2L/dx) (dx = 1 cm).
  (B) part G2: mean direction error vs the rms variation of a slow drift (L = 100 m) about its survey mean,
      added to 0.50 nT white noise, offsets fitted.
"""
import os, json
import numpy as np
import matplotlib.pyplot as plt
import paperstyle as ps


def main():
    a6 = json.load(open(os.path.join(ps.RES_ANALYSIS, "a06_isolated_robustness.json")))
    A3 = a6["A_dc_offsets"]["drive_around"]["A3_correlated"]
    finite = [r for r in A3 if r["L_m"] != "inf"]
    L = np.array([r["L_m"] for r in finite]) * 100.0          # cm
    none = np.array([r["none"]["dir_err_deg"]["mean"] for r in finite])
    fit = np.array([r["rover_world_dc"]["dir_err_deg"]["mean"] for r in finite])
    inf = [r for r in A3 if r["L_m"] == "inf"][0]
    white = a6["A_dc_offsets"]["drive_around"]["A1_iid"]["none"]["dir_err_deg"]["mean"]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(15, 5.8), gridspec_kw={"width_ratios": [1.45, 1]})
    ax.axhline(5, color="#8C8C8C", ls="--", lw=1.8, zorder=1)
    ax.axhline(white, color="black", ls=":", lw=1.4, zorder=1, label="White noise only, 0.5 nT")
    ax.plot(L, none, "-o", color=ps.ORANGE, lw=2.4, ms=8, label="Offsets not fitted", zorder=4)
    ax.plot(L, fit, "-s", color=ps.BLUE, lw=2.4, ms=8, label="Rover and local offsets fitted", zorder=4)
    Ls = np.logspace(np.log10(0.5), np.log10(8), 50)
    ax.plot(Ls, white * np.sqrt(2 * Ls / 1.0), ls="-.", color="#6A3D9A", lw=2.6,
            label=r"$\theta_{\rm white}\sqrt{2L/\Delta x}$", zorder=6)
    xinf = L.max() * 3.2
    ax.set_xscale("log")
    ax.set_xlim(0.35, xinf * 1.8)
    ticks = [1, 10, 100, 1000, 10000]
    ax.set_xticks(ticks + [xinf])
    ax.set_xticklabels(["1 cm", "10 cm", "1 m", "10 m", "100 m", r"$\infty$"])
    ax.set_ylim(0, 92)
    ax.set_xlabel(r"Noise correlation length $L$ (1 cm = 1 s at 1 cm s$^{-1}$)")
    ax.set_ylabel(r"$\vec{M}_{\rm NRM}$ Direction Error ($\degree$)")
    ax.text(120, 6.0, r"$5\degree$ goal", fontsize=13, color="#6E6E6E")
    ax.legend(loc="upper center", bbox_to_anchor=(0.56, 0.995), fontsize=13, ncol=2)
    ps.panel_letter(ax, "A", dark=False)

    G2 = a6["G_gradient_and_drift"]["drive_around"]["drift_L100m"]
    x = np.array([0.0] + [r["rms_about_mean_nT"] for r in G2[1:]])
    y = np.array([r["offsets_fitted"]["dir_err_deg"]["mean"] for r in G2])
    bx.axhline(5, color="#8C8C8C", ls="--", lw=1.8, zorder=1)
    bx.plot(x, y, "-s", color=ps.BLUE, lw=2.4, ms=8, zorder=4)
    bx.set_xlabel("Slow drift about survey mean (nT rms per axis)")
    bx.set_ylabel(r"$\vec{M}_{\rm NRM}$ Direction Error ($\degree$)")
    bx.set_xlim(-0.005, 0.155); bx.set_ylim(0, 16)
    bx.text(0.100, 5.3, r"$5\degree$ goal", fontsize=13, color="#6E6E6E")
    ps.panel_letter(bx, "B", dark=False)
    fig.tight_layout(w_pad=3)

    # Axis break: panel A is split into [finite L | gap | L = inf] within the area given by tight_layout.
    # The L = inf points are drawn in their own narrow axes (shared y) at x = xinf, apart from the curves,
    # which end at L = 100 m.
    pos = ax.get_position()
    full_hi = xinf * 1.8                                    # right limit of the unbroken axis
    ax_hi = L.max() * 1.8                                   # right limit of the finite-L axes
    frac = np.log10(ax_hi / 0.35) / np.log10(full_hi / 0.35)   # same cm-per-decade scale as the unbroken axis
    gap = 0.03
    ax.set_xlim(0.35, ax_hi)
    ax.set_xticks(ticks)
    ax.set_xticklabels(["1 cm", "10 cm", "1 m", "10 m", "100 m"])
    ax.set_position([pos.x0, pos.y0, pos.width * frac, pos.height])
    axi = fig.add_axes([pos.x0 + pos.width * (frac + gap), pos.y0, pos.width * (1 - frac - gap), pos.height], sharey=ax)
    axi.axhline(5, color="#8C8C8C", ls="--", lw=1.8, zorder=1)
    axi.axhline(white, color="black", ls=":", lw=1.4, zorder=1)
    axi.plot([xinf], [inf["none"]["dir_err_deg"]["mean"]], "o", color=ps.ORANGE, ms=9, zorder=5)
    axi.plot([xinf], [inf["rover_world_dc"]["dir_err_deg"]["mean"]], "s", color=ps.BLUE, ms=9, zorder=5)
    axi.set_xlim(xinf * 0.5, xinf * 1.5)
    axi.set_xticks([xinf])
    axi.set_xticklabels([r"$\infty$"])
    axi.tick_params(axis="y", left=False, labelleft=False)
    ax.spines["right"].set_visible(False)
    axi.spines["left"].set_visible(False)
    ax.set_zorder(axi.get_zorder() + 1)                     # keeps the legend (a child of ax) on top
    brk = dict(marker=[(-1, -1), (1, 1)], markersize=13, mec="black", mew=2.0, ls="none", color="black",
               clip_on=False, zorder=30, gid="axisbreak")
    for a, xe in ((ax, 1.0), (axi, 0.0)):                   # break marks where each axes meets the break
        for ye in (0.0, 1.0):
            a.plot([xe], [ye], transform=a.transAxes, **brk)[0].set_in_layout(False)   # keep the figure size unchanged
    # Legend and panel letter at fixed positions on the page
    ax.get_legend().set_bbox_to_anchor((pos.x0 + 0.56 * pos.width, pos.y0 + 0.995 * pos.height), transform=fig.transFigure)
    [t for t in ax.texts if t.get_text() == "A"][0].set_position((0.03 / frac, 0.95))
    ax.xaxis.label.set_x(0.5 / frac)                        # x label centered on the whole panel
    ps.save(fig, "fig08_noise_requirement")


if __name__ == "__main__":
    main()
