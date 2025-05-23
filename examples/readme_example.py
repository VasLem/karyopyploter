from karyopyploter import GENOME, plot_ideogram, make_ideogram_grid, annotate_ideogram, zoom
from matplotlib import pyplot as plt
from itertools import chain
from pathlib import Path

OUT_DIR = Path(__file__).parent.parent / "example_outputs" / "readme_example"
OUT_DIR.mkdir(parents=True, exist_ok=True)
genome = GENOME.CHM13

fig, axes = plt.subplots(
    ncols=1,
    nrows=22,
    figsize=(11, 11),
    facecolor="white",
)
for ax, contig_name in zip(axes, [f"chr{i}" for i in chain(range(1, 23), "XY")]):
    chromosome = contig_name
    plot_ideogram(ax, target=chromosome, genome=genome, label=contig_name)
fig.savefig(OUT_DIR / "ideogram1.png", dpi=300)
# similar to:
fig = plt.figure(figsize=(11, 11), facecolor="white")
fig, _, ideogram_axes = make_ideogram_grid(target=[f"chr{contig_name}" for contig_name in chain(range(1, 23), "XY")], num_subplots=0, genome=genome, fig=fig)
fig.savefig(OUT_DIR / "ideogram2.png", dpi=300)
# and with a subplots grid
fig, ax, ideogram_axes = make_ideogram_grid(subplot_width=15,ideogram_factor=0.2, target=[f"chr{contig_name}" for contig_name in chain(range(1, 23), "XY")], num_subplots=1, genome=genome)
fig.savefig(OUT_DIR / "ideogram3.png", dpi=300)
# and with some regions annotated
regions = {'chr1':[(1000000,2000000, "red")], 'chr2':[(3000000, 4000000, 'blue')], 'chr3':[(5000000,6000000, (0,1,0)), (7000000,8000000, (1,0,0))]}
for chr in regions:
    annotate_ideogram(ideogram_axes[chr], regions=regions[chr], genome=genome)
fig.savefig(OUT_DIR / "ideogram4.png", dpi=300)
# maybe we want to zoom in on specific regions
zoom_regions = {'chr1':(500000,2500000), 'chr10':(3000000, 4000000)}
for chr in zoom_regions:
    zoom(ideogram_axes[chr], start=zoom_regions[chr][0], stop=zoom_regions[chr][1])
fig.savefig(OUT_DIR / "ideogram5.png", dpi=300)