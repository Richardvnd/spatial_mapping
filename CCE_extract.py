import scri
import numpy as np
import pickle
import json

from scri.asymptotic_bondi_data.map_to_superrest_frame import MT_to_WM

"""

This code takes the CCE data and transforms it into the superrest frame.

"""

filenums = [
    "0001",
    "0002",
    "0003",
    "0004",
    "0005",
    "0006",
    "0007",
    "0008",
    "0009",
    "0010",
    "0011",
    "0012",
    "0013",
]
r1 = [
    "0070",
    "0064",
    "0063",
    "0061",
    "0066",
    "0067",
    "0067",
    "0073",
    "0067",
    "0070",
    "0063",
    "0065",
    "0070",
]
r2 = [
    "0292",
    "0261",
    "0250",
    "0236",
    "0274",
    "0273",
    "0270",
    "0305",
    "0270",
    "0235",
    "0222",
    "0223",
    "0237",
]
r3 = [
    "0513",
    "0458",
    "0438",
    "0410",
    "0482",
    "0479",
    "0472",
    "0538",
    "0472",
    "0400",
    "0381",
    "0382",
    "0403",
]
r4 = [
    "0735",
    "0655",
    "0625",
    "0585",
    "0690",
    "0685",
    "0675",
    "0770",
    "0675",
    "0565",
    "0540",
    "0540",
    "0570",
]
rss = [r1, r2, r3, r4]
rssnames = ["R1", "R2", "R3", "R4"]
levs = ["Lev3", "Lev4", "Lev5"]

for lev in levs:
    for rssname, radii in zip(rssnames, rss):
        for fileindex, filenum in enumerate(filenums):

            radius = radii[fileindex]

            print(
                "Level: ",
                lev,
                "Radius name: ",
                rssname,
                "Radius: ",
                radius,
                "File: ",
                filenum,
            )

            # FILEPATHS GO HERE

            filepath = f""
            filepathraw = f""

            abd = scri.SpEC.file_io.create_abd_from_h5(
                h=f"{filepathraw}/{lev}:rhOverM_BondiCce_R{radius}.h5",
                Psi4=f"{filepathraw}/{lev}:rMPsi4_BondiCce_R{radius}.h5",
                Psi3=f"{filepathraw}/{lev}:r2Psi3_BondiCce_R{radius}.h5",
                Psi2=f"{filepathraw}/{lev}:r3Psi2OverM_BondiCce_R{radius}.h5",
                Psi1=f"{filepathraw}/{lev}:r4Psi1OverM2_BondiCce_R{radius}.h5",
                Psi0=f"{filepathraw}/{lev}:r5Psi0OverM3_BondiCce_R{radius}.h5",
                file_format="RPDMB",
            )

            # Define t = 0 using peak luminosity
            abd.t -= abd.t[np.argmax(MT_to_WM(2.0 * abd.sigma.bar.dot).norm())]

            # Interpolate to the merger/ringdown regime (to speed things up)
            abd_ringdown = abd.interpolate(np.arange(-100, abd.t[-1], 0.1))

            # Map to remnant superrest frame
            abd_ringdown_superrest, _, _ = abd_ringdown.map_to_superrest_frame(
                t_0=300, padding_time=20
            )

            # Get mass and spin from abd objects
            M_f = abd_ringdown_superrest.bondi_rest_mass()[-1]
            chi_f = abd_ringdown_superrest.bondi_dimensionless_spin()[-1]
            chi_f = chi_f.tolist()

            metadata = {"M_f": M_f, "chi_f": chi_f}

            # Compute the strain waveform
            h = MT_to_WM(2.0 * abd_ringdown_superrest.sigma.bar)
            h_dict = {"times": h.t}

            for ell in range(2, h.ell_max + 1):
                for m in range(-ell, ell + 1):
                    h_dict[ell, m] = np.array(h.data[:, h.index(ell, m)])

            # Save files

            with open(
                f"{filepath}/superrest_data/SXS:BBH_ExtCCE_superrest:{filenum}/SXS:BBH_ExtCCE_superrest:{filenum}_{lev}_{rssname}.pickle",
                "wb",
            ) as f:
                pickle.dump(h_dict, f)

            with open(
                f"{filepath}/superrest_data/SXS:BBH_ExtCCE_superrest:{filenum}/SXS:BBH_ExtCCE_superrest:{filenum}_{lev}_{rssname}_metadata.json",
                "w",
            ) as f:
                json.dump(metadata, f)
