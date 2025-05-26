import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.axes._secondary_axes import SecondaryAxis
from matplotlib.patches import PathPatch, Rectangle, Patch
from matplotlib.path import Path as MplPath
from matplotlib.typing import ColorType
from typeguard import check_type
from typing import Literal
from typing import Callable, Union

from karyopyploter.utils import get_cytoband_df
from karyopyploter.constants import (
    GENOME,
    DETAIL,
    ORIENTATION,
)


def human_readable_coordinates(coordinate: float) -> str:
    """
    Convert a coordinate to a human-readable format.
    :param coordinate: Coordinate to convert.
    :return: Human-readable coordinate as a string.
    """
    if coordinate >= 1e9:
        return f"{coordinate / 1e9:.2f} Gb"
    elif coordinate >= 1e6:
        return f"{coordinate / 1e6:.2f} Mb"
    elif coordinate >= 1e3:
        return f"{coordinate / 1e3:.2f} Kb"
    else:
        return f"{coordinate:.0f} bp"


def add_ideogram_coordinates(
    ax: Axes,
    coordinates_number: int = 8,
    coordinate_format: Union[str, Callable[[float], str]] = human_readable_coordinates,
    orientation: ORIENTATION = ORIENTATION.HORIZONTAL,
    **tick_params,
):
    if isinstance(coordinate_format, str):
        coordinate_format = lambda x: f"{x:{coordinate_format}}"
    elif not callable(coordinate_format):
        raise ValueError("coordinate_format must be a string or a callable function")
    assert coordinates_number > 0, "coordinates_number must be greater than 0"
    extents = [child.get_extents() for child in ax.get_children() if isinstance(child, PathPatch)]
    if orientation == ORIENTATION.HORIZONTAL:
        lims = min([extent.x0 for extent in extents]), max([extent.x1 for extent in extents])
        lims = ax.transData.inverted().transform([[lims[0], 0], [lims[1], 0]])[:, 0]
        lims = max(lims[0], ax.get_xlim()[0]), min(lims[1], ax.get_xlim()[1])
    else:
        lims = min([extent.y0 for extent in extents]), max([extent.y1 for extent in extents])
        lims = ax.transData.inverted().transform([[0, lims[0]], [0, lims[1]]])[:, 1]
        lims = max(lims[0], ax.get_ylim()[0]), min(lims[1], ax.get_ylim()[1])

    ticks = [lims[0] + (lims[1] - lims[0]) * i / (coordinates_number - 1) for i in range(coordinates_number)]
    labels = [coordinate_format(tick) for tick in ticks]
    if orientation == ORIENTATION.HORIZONTAL:
        ax.xaxis.set_visible(True)
        ax.set_xticks(ticks, labels=labels, **tick_params)
    else:
        ax.yaxis.set_visible(True)
        ax.set_yticks(ticks, labels=labels, **tick_params)


def annotate_ideogram(
    ax: Axes,
    regions: list[tuple[int, int, ColorType]] | None = None,
    lower_anchor: int = 0,
    height: int = 1,
    orientation: ORIENTATION = ORIENTATION.HORIZONTAL,
    **kwargs,
):
    """
    Annotate ideogram regions with rectangles.
    :param ax: The axis to annotate.
    :param regions: List of regions to annotate.
    :param lower_anchor: Lower anchor point for the ideogram, for outline.
    :param height: Height of the ideogram.
    :param orientation: Orientation of the ideogram (horizontal or vertical).
    :param kwargs: Additional keyword arguments for the rectangle patches.
    """
    assert isinstance(regions, list), "Regions must be a list of tuples"
    for region in regions:
        assert len(region) == 3, "Each region must be a tuple of (start, stop, colour)"
        assert isinstance(region[0], int), "Start must be an integer"
        assert isinstance(region[1], int), "Stop must be an integer"
        assert region[0] < region[1], "Start must be less than stop"
        assert check_type(region[2], ColorType), "Third element must be a colour"
    for region in regions:
        start, stop, colour = region
        ax.add_patch(
            Rectangle(
                (start, 0),
                stop - start,
                1,
                color=colour,
                alpha=0.5,
            )
        )
    for r_start, r_stop, r_colour in regions:
        x0, rwidth = r_start, r_stop - r_start
        y0 = lower_anchor + 0.02
        # print(f"x0 {x0}, width {width}, height: he")
        rheight = height
        if orientation == ORIENTATION.VERTICAL:
            x0 = lower_anchor + 0.03
            y0 = r_start
            rheight = rwidth
            rwidth = 0.94

        r = Rectangle(
            (x0, y0),  # +0.01 should shift us off outline of chromosome
            width=rwidth,
            height=rheight,
            fill=kwargs.get("fill", True),
            color=r_colour,
            joinstyle="round",
            zorder=3,
            alpha=kwargs.get("alpha", 0.5),
            lw=kwargs.get("lw", 1),
        )
        ax.add_patch(r)


def zoom(ax: Axes, start: int, stop: int, orientation: ORIENTATION = ORIENTATION.HORIZONTAL):
    """
    Zoom in on the specified region of the ideogram.
    :param ax: The axis to zoom in on.
    :param start: The start position of the region to zoom in on.
    :param stop: The stop position of the region to zoom in on.
    :param orientation: The orientation of the ideogram (horizontal or orientation).
    """
    # Zoom in on the specified region
    if orientation == ORIENTATION.HORIZONTAL:
        if start is None:
            start = ax.get_xlim()[0]
        if stop is None:
            stop = ax.get_xlim()[1]
        ax.set_xlim(start, stop)
    else:
        if start is None:
            start = ax.get_ylim()[0]
        if stop is None:
            stop = ax.get_ylim()[1]
        ax.set_ylim(start, stop)


def plot_ideogram(
    ax: Axes,
    target: str,
    genome: GENOME = GENOME.HG38,
    start: int | None = None,
    stop: int | None = None,
    lower_anchor: int = 0,
    height: int = 1,
    curve: float = 0.02,
    label: str | None = None,
    label_placement: Literal["height", "length"] = "height",
    label_axis_offset: float = None,
    label_params: dict = None,
    show_coordinates: bool = False,
    coordinates_params: dict = None,
    orientation: ORIENTATION = ORIENTATION.HORIZONTAL,
    regions: list[tuple[int, int, ColorType]] | None = None,
    regions_annotation_params: dict = None,
    cytobands_df: pd.DataFrame = None,
    cytobands: DETAIL = DETAIL.CYTOBAND,
    _arrange_absolute_ax_lims: bool = True,
):
    """
    Plot a chromosome ideogram with cytobands and optionally highlight a specific region.

    :param ax: Matplotlib axis object where the ideogram will be plotted.
    :param target: target chromosome to filter and plot.
    :param start: Starting base pair position for the region of interest (optional).
    :param stop: Ending base pair position for the region of interest (optional).
    :param lower_anchor: Lower anchor point for the ideogram, for outline.
    :param height: Height of the ideogram.
    :param curve: Curve factor for the ideogram edges.
    :param label: Label for the ideogram, displayed at the top or side.
    :param label_placement: Placement of the label, either "height" (at the top) or "length" (at the side).
    :param label_axis_offset: Offset for the label axis. If to show coordinates, the default is 0.3, otherwise 0.
    :param label_params: Additional parameters for the label, such as font size or rotation.
    :param show_coordinates: Whether to show coordinates on the ideogram.
    :param coordinates_params: Parameters for the coordinates.
    :param target_region_extent: Extent of the target region highlight.
    :param orientation: orientation of ideogram.
    :param regions: List of regions to colour in on the karyotype. Respects orientation kwarg - a region should
    be a tuple of format (start, stop, colour)
    :param cytobands_df: DataFrame containing cytoband data with columns "chrom", "chromStart",
      "chromEnd", "gieStain", and "colour".
    :param cytobands: Whether to render cytobands

    :return: Updated axis object with the plotted ideogram.

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> ax = plot_ideogram(ax, "chr1", start=50, stop=250, label="Chromosome 1")
    >>> ax.get_xlim()  # To test if the ideogram was plotted (not a direct measure but gives an idea)
    (np.float64(50.0), np.float64(250.0))

    # Test behaviour with a non-existent chromosome
    >>> ax = plot_ideogram(ax, "chr_1", start=50, stop=250, label="Chromosome 1")# doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ValueError: Chromosome chr_1 not found in cytoband data. Should be one of ...

    """
    if coordinates_params is None:
        coordinates_params = dict()
    if label_params is None:
        label_params = dict()
    if regions_annotation_params is None:
        regions_annotation_params = dict()
    if label_axis_offset is None:
        label_axis_offset = 0.3 if show_coordinates else 0
    # some checks for input before we start
    if label is not None:
        assert label_placement in ["height", "length"], "label_placement must be either 'height' or 'length'"
    if start is not None and stop is not None:
        assert start < stop, "Start must be less than stop"

    if cytobands_df is None:
        df = get_cytoband_df(genome)
    else:
        df = cytobands_df
    chr_names = df["chrom"].unique()

    df = df[df["chrom"].eq(target)]
    if df.empty:
        msg = f"Chromosome {target} not found in cytoband data. Should be one of {chr_names}"
        raise ValueError(msg)

    # Beginning with plotting
    yrange = (lower_anchor, height)
    xrange = df[["chromStart", "width"]].values
    chr_start = df["chromStart"].min()
    chr_end = df["chromEnd"].max()
    chr_len = chr_end - chr_start
    ymid = (max(yrange) - min(yrange)) / 2
    if cytobands == DETAIL.CYTOBAND:
        if orientation == ORIENTATION.VERTICAL:
            yranges = df[["chromStart", "width"]].values
            x_range = (lower_anchor, height)
            face_colours = iter(df["colour"])
            for yrange in yranges:
                ax.broken_barh([(lower_anchor, height - 0.01)], yrange, facecolors=next(face_colours), zorder=1)

            (max(x_range) - min(x_range)) / 2

        else:
            ax.broken_barh(xrange, yrange, facecolors=df["colour"], alpha=0.6)

    # Define and draw the centromere using the rows marked as 'cen' in the 'gieStain' column
    cen_df = df[df["gieStain"].str.contains("cen")]
    cen_start = cen_df["chromStart"].min()
    cen_end = cen_df["chromEnd"].max()

    cen_outline = [
        (MplPath.MOVETO, (cen_start, lower_anchor)),
        (MplPath.LINETO, (cen_start, height)),
        (MplPath.LINETO, ((cen_start + cen_end) / 2, ymid)),
        (MplPath.LINETO, (cen_end, height)),
        (MplPath.LINETO, (cen_end, lower_anchor)),
        (MplPath.LINETO, ((cen_start + cen_end) / 2, ymid)),
        (MplPath.CLOSEPOLY, (cen_start, lower_anchor)),
    ]
    chr_end_without_curve = chr_end - chr_len * curve
    chr_start_without_curve = chr_start + chr_len * curve
    # Define and draw the chromosome outline, taking into account the shape around the centromere
    outline = [
        (MplPath.MOVETO, (chr_start_without_curve, height)),
        # Top part
        (MplPath.LINETO, (cen_start, height)),
        (MplPath.LINETO, ((cen_start + cen_end) / 2, ymid)),
        (MplPath.LINETO, (cen_end, height)),
        (MplPath.LINETO, (chr_end_without_curve, height)),
        (MplPath.CURVE3, (chr_end, height)),
        (MplPath.CURVE3, (chr_end, ymid)),
        # Bottom part
        (MplPath.CURVE3, (chr_end, lower_anchor)),
        (MplPath.CURVE3, (chr_end_without_curve, lower_anchor)),
        (MplPath.LINETO, (chr_end_without_curve, lower_anchor)),
        (MplPath.LINETO, (cen_end, lower_anchor)),
        (MplPath.LINETO, ((cen_start + cen_end) / 2, ymid)),
        (MplPath.LINETO, (cen_start, lower_anchor)),
        (MplPath.LINETO, (chr_start_without_curve, lower_anchor)),
        (MplPath.CURVE3, (chr_start, lower_anchor)),
        (MplPath.CURVE3, (chr_start, ymid)),
        (MplPath.CURVE3, (chr_start, height)),
        (MplPath.CURVE3, (chr_start_without_curve, height)),
        (MplPath.MOVETO, (chr_start_without_curve, height)),
    ]

    def invert_with_curve(outline):
        outline = outline[::-1]
        new_outline = outline.copy()
        i = 0
        while i < len(outline):
            if outline[i][0] == MplPath.CURVE3:
                j = i + 1
                while j < len(outline) and outline[j][0] == MplPath.CURVE3:
                    j += 1
                new_outline[i:j] = outline[i + 1 : j] + [(MplPath.CURVE3, outline[j][1])]
                i = j + 1
            else:
                i += 1
        return new_outline

    outside_outline = [
        (MplPath.MOVETO, (chr_start, height)),
        (MplPath.LINETO, (chr_end, height)),
        (MplPath.LINETO, (chr_end, lower_anchor)),
        (MplPath.LINETO, (chr_start, lower_anchor)),
        (MplPath.CLOSEPOLY, (chr_start, lower_anchor)),
    ] + invert_with_curve(outline)
    if orientation == ORIENTATION.VERTICAL:
        outline = [(command, coords[::-1]) for command, coords in outline]
        cen_outline = [(command, coords[::-1]) for command, coords in cen_outline]
        outside_outline = [(command, coords[::-1]) for command, coords in outside_outline]
    cen_move, cen_poly = zip(
        *cen_outline,
        strict=True,
    )
    cen_patch = PathPatch(MplPath(cen_poly, cen_move), facecolor=(0.8, 0.4, 0.4), lw=0, alpha=1, zorder=2)

    ax.add_patch(cen_patch)

    chr_move, chr_poly = zip(
        *outline,
        strict=True,
    )
    mask_move, mask_poly = zip(
        *outside_outline,
        strict=True,
    )
    mask_patch = PathPatch(
        MplPath(mask_poly, mask_move), facecolor=(1.0, 1.0, 1.0), alpha=1, edgecolor=(1.0, 1.0, 1.0), zorder=2
    )
    ax.add_patch(mask_patch)
    chr_patch = PathPatch(MplPath(chr_poly, chr_move), fill=None, joinstyle="round", alpha=1, zorder=2)
    ax.add_patch(chr_patch)
    # If start and stop positions are provided, draw a rectangle to highlight this region
    if start is not None or stop is not None:
        zoom(ax, start, stop, orientation=orientation)
    else:
        if _arrange_absolute_ax_lims:
            if orientation == ORIENTATION.HORIZONTAL:
                ax.set_xlim(min(ax.get_xlim()[0], chr_start), max(ax.get_xlim()[1], chr_end))
            else:
                ax.set_ylim(min(ax.get_ylim()[0], chr_start), max(ax.get_ylim()[1], chr_end))

    if regions:
        annotate_ideogram(
            ax, regions, height=height, lower_anchor=lower_anchor, orientation=orientation, **regions_annotation_params
        )

    if orientation == ORIENTATION.VERTICAL:
        ax.set_xlim(lower_anchor - 0.05, height + 0.05)
    else:
        ax.set_ylim(lower_anchor - 0.05, height + 0.05)

    # Remove axis spines and ticks for a cleaner look
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    if show_coordinates:
        add_ideogram_coordinates(
            ax,
            **coordinates_params,
            orientation=orientation,
        )

    def get_secondary_axis(ax, which: str, label_axis_offset=0):
        for x in ax.get_children():
            if isinstance(x, SecondaryAxis):
                if which == "x" and x._loc == "bottom":
                    return x
                if which == "y" and x._loc == "left":
                    return x
        if which == "x":
            return ax.secondary_xaxis(-label_axis_offset if label_axis_offset else "bottom")
        else:
            return ax.secondary_yaxis(-label_axis_offset if label_axis_offset else "left")

    # Add chromosome name to the plot
    if label is not None:
        if label_placement == "height":
            to_place = height / 2
            if orientation == ORIENTATION.VERTICAL:
                sec = get_secondary_axis(ax, "x", label_axis_offset=label_axis_offset)
                labs = sec.get_xticklabels()
                locs = sec.get_xticks()
            else:
                sec = get_secondary_axis(ax, "y", label_axis_offset=label_axis_offset)
                labs = sec.get_yticklabels()
                locs = sec.get_yticks()
        elif label_placement == "length":
            to_place = (chr_start + chr_end) / 2
            if orientation == ORIENTATION.VERTICAL:
                sec = get_secondary_axis(ax, "y", label_axis_offset=label_axis_offset)
                labs = sec.get_yticklabels()
                locs = sec.get_yticks()
            else:
                sec = get_secondary_axis(ax, "x", label_axis_offset=label_axis_offset)
                labs = sec.get_xticklabels()
                locs = sec.get_xticks()

        def is_number(s):
            try:
                float(s)
                return True
            except ValueError:
                return False

        tk = [
            i
            for i, (l, x) in enumerate(zip(labs, locs))
            if not is_number(l.get_text()) or round(float(x), 2) != round(float(l.get_text()), 2)
        ]
        labs = [labs[i] for i in tk]
        locs = [locs[i] for i in tk]

        x = [i for i, (l, u) in enumerate(zip(locs[:-1], locs[1:])) if to_place > l and to_place <= u]
        if x:
            pos = x[0]
        else:
            if locs and to_place > locs[-1]:
                pos = len(locs)
            else:
                pos = 0
        locs.insert(pos, to_place)
        labs.insert(pos, label)

        if label_placement == "height":
            if orientation == ORIENTATION.VERTICAL:
                sec.set_xticks(locs, labs)
                if label_params:
                    plt.setp(sec.get_xticklabels()[pos], **label_params)
                sec.spines["bottom"].set_visible(False)
            else:
                sec.set_yticks(locs, labs)
                if label_params:
                    plt.setp(sec.get_yticklabels()[pos], **label_params)
                sec.spines["left"].set_visible(False)
        else:
            if orientation == ORIENTATION.VERTICAL:
                sec.set_yticks(locs, labs)
                if label_params:
                    plt.setp(sec.get_yticklabels()[pos], **label_params)
                sec.spines["left"].set_visible(False)
            else:
                sec.set_xticks(locs, labs)
                if label_params:
                    plt.setp(sec.get_xticklabels()[pos], **label_params)
                sec.spines["bottom"].set_visible(False)
        sec.tick_params(axis=u'both', which=u'both', length=0)

    return ax


if __name__ == "__main__":
    fig, axes = plt.subplots(
        ncols=1,
        nrows=22,
        figsize=(11, 11),
        facecolor="white",
    )
    genome = GENOME.CHM13
    for ax, contig_name in zip(axes, range(1, 23), strict=False):
        chromosome = f"chr{contig_name}"
        plot_ideogram(ax, target=chromosome, genome=genome)
    fig.savefig("ideogram.png", dpi=300)
