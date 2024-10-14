import pickle
import json
import qnmfits

# FILEPATH GOES HERE

filepath = "/data/rvnd2-2/CCE_data/superrest_data"


def SXS_CCE(ID, zero_time=(2, 2), lev="Lev5", radius="R2"):
    """
    This function reads in the transformed CCE data and returns a qnmfits.Custom object.

    Parameters
    ----------
    ID : str
        The simulation ID.
    zero_time : tuple, optional
        The (ell,m) mode to use as the zero time. The default is (2,2).
    lev : str, optional
        The simulation level. The default is 'Lev5'.
    radius : str, optional
        The simulation radius. The default is 'R2'.

    Returns
    -------
    sim : qnmfits.Custom
        The qnmfits.Custom object containing the CCE data.

    """

    with open(
        f"{filepath}/SXS:BBH_ExtCCE_superrest:{ID}/SXS:BBH_ExtCCE_superrest:{ID}_{lev}_{radius}.pickle",
        "rb",
    ) as f:
        h_prime_dict = pickle.load(f)
    with open(
        f"{filepath}/SXS:BBH_ExtCCE_superrest:{ID}/SXS:BBH_ExtCCE_superrest:{ID}_{lev}_{radius}_metadata.json",
        "r",
    ) as f:
        metadata = json.load(f)

    times = h_prime_dict.pop("times")

    sim = qnmfits.Custom(
        times,
        h_prime_dict,
        metadata={
            "remnant_mass": metadata["M_f"],
            "remnant_dimensionless_spin": metadata["chi_f"],
        },
        zero_time=zero_time,
    )

    return sim
