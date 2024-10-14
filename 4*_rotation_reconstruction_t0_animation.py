import matplotlib.pyplot as plt
import numpy as np
import qnmfits
import CCE
from matplotlib.animation import FuncAnimation
from scipy.interpolate import UnivariateSpline
from qnmfits.spatial_mapping_functions import *

plt.style.use("stylesheet.mplstyle")

"""
This code generates an animated reconstruction across a range of t0 values.

"""

id = "0001"
sim = CCE.SXS_CCE(id, lev="Lev5", radius="R2")

l_max = 8
n_max = 7

tstart = -100
tend = 101
tstep = 1
times_list = []

lon = np.linspace(-np.pi, np.pi, 200)
lat = np.linspace(-np.pi / 2, np.pi / 2, 200)
Lon, Lat = np.meshgrid(lon, lat)

QNMs = [(lam, 2, n, p) for lam in np.arange(2, l_max + 1) for n in np.arange(0, n_max + 1) for p in (-1, +1)]

spherical_modes = [(l, 2) for l in np.arange(2, l_max + 1)]

sm_list = []
arg_list = []

mapping = [(2, 2, 0, 1)]
map = mapping[0]

fig = plt.figure(figsize=(6, 5.5), dpi=100)
gs = fig.add_gridspec(3, 3, height_ratios=[0.5, 1, 1])

axs_small = fig.add_subplot(gs[0, 1], projection="mollweide")
axs_large = [fig.add_subplot(gs[i + 1, :]) for i in range(2)]

axs_small.set_xticklabels([])
axs_small.set_yticklabels([])

axs_large[0].tick_params(axis="x", which="both", labelbottom=False)
axs_large[0].set_yscale("log")
axs_large[0].set_ylabel("$\mathcal{M}_s$")

f_actual = qnmfits.qnm.omega_list(mapping, sim.chif_mag, sim.Mf)
axs_large[1].axhline(np.real(f_actual), color="b", ls="--", lw=1)

axs_large[1].set_xlabel("$t_0 \ [M]$")
axs_large[1].set_ylabel("$\mathrm{dArg}(\mathcal{M}_s)/\mathrm{d}t_0$")

axs_large[0].set_xlim(tstart, tend - 1)
axs_large[1].set_xlim(tstart, tend - 1)
axs_large[0].set_ylim(1e-7, 1)
axs_large[1].set_ylim(0, 0.6)


def update(step):

    if step in times_list:
        return fig
    else:
        times_list.append(step)

        best_fit = qnmfits.mapping_multimode_ringdown_fit(
            sim.times,
            sim.h,
            modes=QNMs,
            Mf=sim.Mf,
            chif=sim.chif_mag,
            t0=step,
            mapping_modes=mapping,
            spherical_modes=spherical_modes,
        )

        F = spatial_reconstruction(np.pi / 2 - Lat, Lon, best_fit, map, l_max)
        sm, arg, _ = spatial_mismatch_linear(best_fit, map, sim.chif_mag, l_max)
        sm_list.append(sm)
        arg_list.append(arg)

        axs_small.pcolormesh(Lon, Lat, np.real(F), cmap=plt.cm.jet)
        axs_large[0].plot(times_list, sm_list, color="black")

        if len(times_list) > 4:
            uvsarg = UnivariateSpline(times_list, -np.unwrap(arg_list), k=4, s=0)
            axs_large[1].plot(times_list, uvsarg.derivative()(times_list), color="black")

        return fig


ani = FuncAnimation(fig, update, frames=range(tstart, tend, tstep), interval=100)
ani.save(f"mapping_animation_{map}_{id}.mp4", writer="ffmpeg", dpi=500)
