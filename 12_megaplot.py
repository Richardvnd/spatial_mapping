import numpy as np
import matplotlib.pyplot as plt
import qnmfits
import CCE
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

plt.style.use("stylesheet.mplstyle")

"""

This code generates a figure illustrating the angle averaged multimode QNM fit to NR data.

"""

sim = CCE.SXS_CCE("0013", lev="Lev5", radius="R2")

ell_max, n_max = 4, 1
show_fundamental = True
t0 = 10.0
modes = [
    (l, m, n, p)
    for l in range(2, ell_max + 1)
    for m in range(-l, l + 1)
    for n in range(n_max + 1)
    for p in [-1, 1]
]
spherical_modes = [(l, m) for l in range(2, ell_max + 1) for m in range(-l, l + 1)]

best_fit = qnmfits.multimode_ringdown_fit(
    sim.times,
    sim.h,
    modes,
    sim.Mf,
    sim.chif_mag,
    t0,
    T=100,
    spherical_modes=spherical_modes,
)

fig, axs = plt.subplots(
    ell_max - 1,
    2 * ell_max + 1,
    figsize=(0.85 * 3 * (2 * ell_max + 1), 0.83 * 4.4 * (ell_max - 2)),
)

for ell in np.arange(2, ell_max):
    num_del = ell_max - ell
    for count in range(num_del):
        axs[ell - 2, count].axis("off")
        axs[ell - 2, 2 * ell_max - count].axis("off")

colors = plt.cm.get_cmap("viridis", 4)
colors2 = plt.cm.get_cmap("YlOrBr", 12)
C0 = colors(0)
C1 = colors(2)
C2 = colors2(3)
C3 = colors(3)

C0patch = mpatches.Patch(color=C0, label="NR data")
C1patch = mpatches.Patch(color=C1, label="QNM model")
kpatch = mpatches.Patch(color="k", label="Residuals\n(shifted and\nx10 for clarity)")
rpatch = mpatches.Patch(color=C2, label="Fundamental QNM\n(shifted for clarity)")
if show_fundamental:
    axs[1, -1].legend(
        handles=[C0patch, C1patch, kpatch, rpatch],
        frameon=False,
        fontsize=12,
        loc="upper right",
    )
else:
    axs[1, -1].legend(
        handles=[C0patch, C1patch, kpatch],
        frameon=False,
        fontsize=12,
        loc="upper right",
    )

line_re = Line2D([0], [0], label="Real part", color=C0, lw=2, ls="-")
line_im = Line2D([0], [0], label="Imag part", color=C0, lw=2, ls=":")
axs[0, -1].legend(
    handles=[line_re, line_im], frameon=False, fontsize=12, loc="lower left"
)

for ell in np.arange(2, ell_max + 1):
    for m in np.arange(-ell, ell + 1):
        ax = axs[ell - 2, m + ell_max]
        y_max = {2: 0.25, 3: 0.05, 4: 0.02, 5: 0.01, 6: 0.004, 7: 0.002, 8: 0.001}[ell]

        ax.plot(sim.times, sim.h[ell, m].real, c=C0, ls="-", lw=2)
        ax.plot(sim.times, sim.h[ell, m].imag, c=C0, ls=":", lw=2)

        ax.plot(
            best_fit["model_times"], best_fit["model"][ell, m].real, c=C1, ls="-", lw=2
        )
        ax.plot(
            best_fit["model_times"], best_fit["model"][ell, m].imag, c=C1, ls=":", lw=2
        )

        ax.plot(
            best_fit["model_times"],
            -y_max
            + 10 * (best_fit["model"][ell, m].real - best_fit["data"][ell, m].real),
            c="k",
            ls="-",
            lw=1,
        )
        ax.plot(
            best_fit["model_times"],
            -y_max
            + 10 * (best_fit["model"][ell, m].imag - best_fit["data"][ell, m].imag),
            c="k",
            ls=":",
            lw=1,
        )

        if m == 2 and show_fundamental:
            idx = [label == f"(2, 2, 0, 1)" for label in best_fit["mode_labels"]]
            # C = best_fit['C'][idx] * qnmfits.qnm.mu(ell, m, 2, 2, 0, 1, sims[2].chif_mag)
            # om = qnmfits.qnm.omega(2, 2, 0, 1, sims[2].chif_mag)
            C = best_fit["C"][idx] * qnmfits.qnm.mu(ell, m, 2, 2, 0, 1, sim.chif_mag)
            om = qnmfits.qnm.omega(2, 2, 0, 1, sim.chif_mag)
            t = best_fit["model_times"]
            ax.plot(
                t,
                y_max + np.real(C * np.exp(-1j * om * (t - t0))),
                c=C2,
                ls="-",
                lw=1,
                label=r"$C_{220+}\mu^{" + str(ell) + str(m) + "}_{220+}$",
            )
            ax.plot(
                t, y_max + np.imag(C * np.exp(-1j * om * (t - t0))), c=C2, ls=":", lw=1
            )
            ax.legend(
                frameon=False,
                fontsize=9,
                loc="lower right",
                bbox_to_anchor=(1.05, 0.67),
            )
            ax.set_facecolor((1, 0, 0, 0.07))

        ax.set_xlim(-30, 90)
        ax.set_xticks([-30, 0, 30, 60, 90])

        if show_fundamental:
            ax.set_ylim(-1.22 * y_max, 1.4 * y_max)
        else:
            ax.set_ylim(-1.22 * y_max, 1.02 * y_max)
        ax.set_yticks([-y_max, 0, y_max])

        if m == 0:
            ax.set_title(r"$\ell\!=\!{},m\!=\!{:}$".format(ell, m))
        else:
            ax.set_title(r"$\ell\!=\!{},m\!=\!{:+}$".format(ell, m))

plt.tight_layout()
# plt.show()
fig.savefig("figs/mega_plot.pdf")
