import re
import numpy as np
from matplotlib.axes import Axes
from karyopyploter.constants import CHR_PATT, COLOUR_LOOKUP, STATIC_PATH
from karyopyploter.constants import GENOME
from pathlib import Path
import pandas as pd

def set_xmargin(ax: Axes, left: float = 0.1, right: float = 0.0) -> None:
    """
    Adjust the x-axis margin of a given Axes object.

    :param ax: The Axes object to modify.
    :param left: Factor by which to expand the left x-axis limit.
    :param right: Factor by which to expand the right x-axis limit.
    
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.set_xlim(0, 10)
    (0.0, 10.0)
    >>> set_xmargin(ax)
    >>> ax.get_xlim()
    (np.float64(-1.0), np.float64(10.0)
    """

    # https://stackoverflow.com/a/49382894
    ax.set_xmargin(0)
    ax.autoscale_view()
    lim = ax.get_xlim()
    delta = lim[1] - lim[0]
    left = lim[0] - delta * left
    right = lim[1] + delta * right
    ax.set_xlim(left, right)


def set_ymargin(ax: Axes, top: float = 0.1, bottom: float = 0.0) -> None:
    """
    Adjust the y-axis margin of a given Axes object.

    :param ax: The Axes object to modify.
    :param top: Factor by which to expand the top y-axis limit.
    :param bottom: Factor by which to expand the bottom y-axis limit.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax.set_ylim(0, 10)
    (0.0, 10.0)
    >>> set_ymargin(ax)
    >>> ax.get_ylim()
    (np.float64(0.0), np.float64(11.0))
    """

    # https://stackoverflow.com/a/49382894
    ax.set_ymargin(0)
    ax.autoscale_view()
    lim = ax.get_ylim()
    delta = lim[1] - lim[0]
    bottom = lim[0] - delta * bottom
    top = lim[1] + delta * top
    ax.set_ylim(bottom, top)



def get_cytobands(genome: GENOME) -> Path:
    """
    Return the cytobands file for the given genome.

    :param genome: The genome variant to get the cytobands file for.
    :return: The path to the cytobands file associated with the provided genome variant.
    :raises ValueError: If the provided genome variant is not recognized.

    >>> get_cytobands(GENOME.HG38) # doctest: +ELLIPSIS
    PosixPath('/.../src/karyopyploter/static/cytobands_hg38.bed')
    >>> get_cytobands(GENOME.CHM13) # doctest: +ELLIPSIS
    PosixPath('/.../src/karyopyploter/static/cytobands_chm13.bed')
    >>> get_cytobands(GENOME.HS1) # doctest: +ELLIPSIS
    PosixPath('/.../src/karyopyploter/static/cytobands_chm13.bed')
    >>> get_cytobands("invalid_genome")
    Traceback (most recent call last):
    ...
    ValueError: Unknown genome: invalid_genome
    """
    match genome:
        case GENOME.HG19:
            return STATIC_PATH / "cytobands_hg19.bed"
        case GENOME.HG38:
            return STATIC_PATH / "cytobands_hg38.bed"
        case GENOME.CHM13:
            return STATIC_PATH / "cytobands_chm13.bed"
        case GENOME.HS1:
            return STATIC_PATH / "cytobands_chm13.bed"
        case _:
            msg = f"Unknown genome: {genome}"
            raise ValueError(msg)
        
def chr_to_ord(x: str):
    out = re.match(CHR_PATT, x.lower())
    if out[3]:
        return (out[3], 0, "")
    chr = out[1].lower()
    if chr == "x":
        chr = 23
    elif chr == "y":
        chr = 24
    elif chr == "m":
        chr = 25
    if out[2]:
        return ("chr", int(chr), out[2])
    return (".", int(chr), "")  # make sure the canonical part stays on top

def get_cytoband_df(genome: GENOME, relative: bool = True) -> pd.DataFrame:
    """
    Convert the cytogram file for the given genome into a dataframe.
    :param genome: The genome to plot the ideogram for.
    :return: A DataFrame containing chromosome cytoband DETAILs.

    >>> dummy_genome = GENOME.HG38  # replace with a test genome path or identifier
    >>> result_df = get_cytoband_df(dummy_genome)
    >>> result_df["chrom"].tolist()[:2]
    ['chr1', 'chr1']
    >>> result_df["chromStart"].tolist()[:2]
    [0, 2300000]
    >>> result_df["arm"].tolist()[:2]
    ['p', 'p']
    """
    cytobands = pd.read_csv(
        get_cytobands(genome), sep="\t", names=["chrom", "chromStart", "chromEnd", "name", "gieStain"]
    )
    cytobands["arm"] = cytobands["name"].str[0]
    cytobands["colour"] = cytobands["gieStain"].map(COLOUR_LOOKUP)
    cytobands["width"] = cytobands["chromEnd"] - cytobands["chromStart"]
    if relative:
        return cytobands
    # Sort the chromosomes in a canonical order
    cytobands = cytobands.sort_values(
        by="chrom",
        key=lambda x: x.map(chr_to_ord),
    )
    # Add a column for the cumulated sum of the chromosome length..
    cumsumlengths = (cytobands.groupby("chrom",sort=False)["chromEnd"].max() - cytobands.groupby("chrom",sort=False)["chromStart"].min()).cumsum()
    # .. to produce the absolute offset of each chromosome
    offset = pd.Series(np.hstack([[0], cumsumlengths[:-1]]), cumsumlengths.index, name="offset")
    cytobands = cytobands.merge(
        offset,
        left_on="chrom",
        right_index=True,
    )
    cytobands["chromStart"] += cytobands["offset"]
    cytobands["chromEnd"] += cytobands["offset"]
    cytobands = cytobands.drop(columns=["offset"])
    return cytobands
