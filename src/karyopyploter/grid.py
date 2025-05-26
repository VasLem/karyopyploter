"""Ideogram plotting executables."""

from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.gridspec import SubplotSpec
from matplotlib.gridspec import GridSpec as gs
from matplotlib.gridspec import GridSpecFromSubplotSpec as gsFromSubplotSpec
from math import ceil
from typing import Callable, Dict, Union, List, Optional
from karyopyploter.ideogram import plot_ideogram
from karyopyploter.constants import GENOME
from karyopyploter.utils import chr_to_ord, get_cytoband_df


def _make_target_grid(
    target: Union[str, List[str]],
    target_stop: str = None,
    genome: GENOME = GENOME.HG38,
    start: int | None = None,
    stop: int | None = None,
    num_subplots=1,
    subplot_width=6,
    height_ratio=0.5,
    ideogram_factor: float = 0.1,
    fig: Optional[plt.Figure] = None,
    subplot_spec: Optional[SubplotSpec] = None,
    ideogram_params: Optional[dict] = None,
) -> tuple[plt.Figure, list[Axes], Axes]:
    """
    Create a grid of subplots, with an ideogram at the bottom. Meant to plot multiple features on the same chromosome.
    :param target: (Starting) target chromosome to filter and plot.
    :param target_stop: Ending target chromosome to filter and plot (optional).
    :param genome: Genome variant to use.
    :param start: Starting base pair position for the region of interest (optional). If start is None, stop must also be None.
    :param stop: Ending base pair position for the region of interest (optional). If stop is None, start must also be None.
    :param num_subplots: Number of subplots to create. If 0, no subplots are created, only the ideogram axis is made.
    :param subplot_width: Width of each subplot.
    :param height_ratio: Height ratio for the subplots.
    :param ideogram_factor: Height factor for the ideogram.
    :param fig: Figure object to use for plotting. If None, a new figure is created.
    :param ideogram_params: Additional keyword arguments for the ideogram plotting function.
    :return: A tuple containing the figure and a list of axes for the subplots.
    """
    if target_stop is None:
        target_stop = target
    if ideogram_params is None:
        ideogram_params = dict()
    relative = ideogram_params.get("relative", None)
    if relative is None:
        relative = target == target_stop
    cytobands_df = None
    chr_start = None
    chr_end = None
    if target != target_stop:
        cytobands_df = get_cytoband_df(genome, relative=False)
        chr_names = cytobands_df["chrom"].unique()
        chr_names = sorted(chr_names, key=chr_to_ord)
        targets = chr_names[chr_names.index(target) : chr_names.index(target_stop) + 1]
        cytobands_df = cytobands_df[cytobands_df['chrom'].isin(targets)]
        chr_start = cytobands_df['chromStart'].min()
        chr_end = cytobands_df['chromEnd'].max()
        if relative:
            cytobands_df['chromEnd'] = cytobands_df['chromEnd'] - cytobands_df['chromStart'].min()
            cytobands_df['chromStart'] = cytobands_df['chromStart'] - cytobands_df['chromStart'].min()
    else:
        targets = [target]
    pfactor = int(1 / ideogram_factor)
    axes = []

    if fig is None:
        fig = plt.figure(
            figsize=(
                subplot_width,
                subplot_width * ((height_ratio * num_subplots) if num_subplots > 0 else height_ratio),
            ),
            facecolor="white",
        )
    if num_subplots > 0:
        if subplot_spec is None:
            gspec = gs(pfactor * num_subplots + 2 * num_subplots - 1, 1, hspace=0.05, wspace=0.05)
        else:
            gspec = gsFromSubplotSpec(
                pfactor * num_subplots + 2 * num_subplots - 1, 1, subplot_spec=subplot_spec, hspace=0.05, wspace=0.05
            )
        for i in range(num_subplots):
            ax = fig.add_subplot(gspec[pfactor * i + i : pfactor * (i + 1) + i, 0])
            axes.append(ax)
        for ax in axes[1:]:
            ax.sharex(axes[0])
            ax.set_xlabel("")
        ideogram_ax = fig.add_subplot(gspec[-num_subplots:, 0], sharex=axes[0])
    else:
        if subplot_spec is None:
            ideogram_ax = fig.add_subplot()
        else:
            ideogram_ax = fig.add_subplot(subplot_spec)
    ideogram_ax.set_xticks([])
    ideogram_ax.set_xticklabels([])
    ideogram_ax.set_xlabel("")
    for cnt, target in enumerate(targets):
        ideogram_params.update(
            {
                "target": target,
                "genome": genome,
                "label": target,
                "label_placement": ideogram_params.get("label_placement", "height" if len(targets) == 1 else "length"),
                "start": None,
                "stop": None,
                "relative": ideogram_params.get("relative", False),
            }
        )
        ideogram_ax = plot_ideogram(
            ideogram_ax, cytobands_df=cytobands_df, _arrange_absolute_ax_lims=False, **ideogram_params
        )
    if start is None:
        start = chr_start
    if stop is None:
        stop = chr_end
    # for obj in ideogram_ax.get_children():
    #     if hasattr(obj, "set_clip_on"):
    #         obj.set_clip_on(False)
    if num_subplots > 0:
        for ax in axes:
            ax.set_xlim(ideogram_ax.get_xlim())
        axes[-1].spines["bottom"].set_visible(False)
        for ax in axes:
            plt.setp(ax.get_xticklabels(), visible=False)
        axes[-1].tick_params(axis=u'both', which=u'both', length=0)
    return fig, axes, ideogram_ax


def make_ideogram_grid(
    target: Union[str, List[str]],
    genome: GENOME = GENOME.HG38,
    start: Union[str, Dict[str, int]] | None = None,
    stop: Union[str, Dict[str, int]] | None = None,
    num_subplots: int = 1,
    subplot_width: float = 10,
    height_ratio: Optional[float] = None,
    fig: Optional[plt.Figure] = None,
    ideogram_factor: float = 0.1,
    grid_params: Optional[dict] = None,
    ideogram_params: Optional[dict] = None,
) -> tuple[plt.Figure, Dict[str, list[Axes]], Dict[str, Axes]]:
    """
    Create a grid of subplots, with an ideogram at the bottom. Meant to plot multiple features on the same chromosome.
    :param target: Target chromosome(s) to plot.
    :param genome: Genome variant to use.
    :param start: Starting base pair position for the region of interest per target(optional). It must be a dictionary if multiple targets are provided. If start is not given for a target, stop must also not be given.
    :param stop: Ending base pair position for the region of interest per target(optional). It must be a dictionary if multiple targets are provided. If stop is not given for a target, start must also not be given.
    :param num_subplots: Number of subplots to create. Can be 0 for no subplots and just ideograms.
    :param show_coordinates: Whether to show coordinates on the ideogram.
    :param subplot_width: Width of each subplot, in inches.
    :param coordinates_number: Number of coordinates to show on the ideogram.
    :param coordinate_format: Format of the coordinates to show on the ideogram.
    :param height_ratio: Height ratio for the subplots, defaults to 0.5 when number of subplots is larger than 1, otherwise 0.1.
    :param ideogram_factor: Height factor for the ideogram.
    :param grid_params: Dictionary of grid parameters for the GridSpec.
    :param ideogram_params: Dictionary of ideogram parameters for the ideogram plotting function.
    :return: A tuple containing the figure, a dictionary of list of subplot axes for each target, and a dictionary of ideogram axes for each target.
    """
    targets = target if isinstance(target, list) else [target]
    if grid_params is None:
        grid_params = dict()
    grid_params.update(
        {
            "hspace": grid_params.get('hspace', 0.05),
            "top": grid_params.get('top', 0.95),
            "bottom": grid_params.get('bottom', 0.05),
            "left": grid_params.get('left', 0.1),
            "right": grid_params.get('right', 1 - 0.1),
        }
    )
    if height_ratio is None:
        if num_subplots > 1:
            height_ratio = 0.5
        else:
            height_ratio = 0.1
    if len(targets) > 1:
        if isinstance(start, int):
            raise ValueError("If multiple targets are provided, start must be a dictionary")
        if isinstance(stop, int):
            raise ValueError("If multiple targets are provided, stop must be a dictionary")
    else:
        start = {targets[0]: start} if isinstance(start, int) else start
        stop = {targets[0]: stop} if isinstance(stop, int) else stop
    if start is None:
        start = dict()
    if stop is None:
        stop = dict()

    start = {t: start.get(t, None) for t in targets}
    stop = {t: stop.get(t, None) for t in targets}
    nrows = ncols = None
    if 'nrows' in grid_params:
        nrows = grid_params["nrows"]
        del grid_params["nrows"]
    if 'ncols' in grid_params:
        ncols = grid_params["ncols"]
        del grid_params["ncols"]
    if ncols is None and nrows is not None:
        ncols = ceil(len(targets) / nrows)
    if nrows is None and ncols is not None:
        nrows = ceil(len(targets) / ncols)
    if ncols is None and nrows is None:
        nrows = len(targets)
        ncols = 1
    if fig is None:
        fig = plt.figure(
            figsize=(
                subplot_width * ncols,
                subplot_width * (height_ratio * num_subplots if num_subplots else height_ratio) * nrows,
            ),
            facecolor="white",
        )

    gs0 = gs(nrows, ncols, figure=fig, **grid_params)
    gs0 = [gs0[i, j] for i in range(nrows) for j in range(ncols)]

    value_axes = {}
    ideogram_axes = {}
    for i, target in enumerate(targets):
        _, a0, a1 = _make_target_grid(
            target,
            genome=genome,
            start=start[target],
            stop=stop[target],
            num_subplots=num_subplots,
            subplot_width=subplot_width,
            height_ratio=height_ratio,
            ideogram_factor=ideogram_factor,
            subplot_spec=gs0[i],
            fig=fig,
            ideogram_params=ideogram_params,
        )
        value_axes[target] = a0
        ideogram_axes[target] = a1
    return fig, value_axes, ideogram_axes


def make_genome_grid(
    target_start: str,
    target_stop: str,
    genome: GENOME = GENOME.HG38,
    num_subplots=1,
    subplot_width: float = 10,
    height_ratio: float = 0.5,
    ideogram_factor: float = 0.1,
    ideogram_params: Dict = None,
) -> tuple[plt.Figure, list[Axes], Axes]:
    """
    Create a grid of subplots for a specific genome with a specified start and stop target.
    :param target_start: Starting target, e.g., chr1.
    :param target_stop: Ending target, e.g., chr5.
    :param genome: Genome variant to use.
    :param num_subplots: Number of subplots to create.
    :param subplot_width: Width of each subplot, in inches.
    :param height_ratio: Height to width ratio for the subplots.
    :param ideogram_factor: Height factor for the ideogram.
    :param ideogram_params: Additional keyword arguments for the ideogram plotting function.
    :return: A tuple containing the figure and a list of axes for the subplots.
    """
    fig, axes, genome_ax = _make_target_grid(
        target=target_start,
        target_stop=target_stop,
        genome=genome,
        num_subplots=num_subplots,
        subplot_width=subplot_width,
        height_ratio=height_ratio,
        ideogram_factor=ideogram_factor,
        ideogram_params=ideogram_params,
    )
    return fig, axes, genome_ax
