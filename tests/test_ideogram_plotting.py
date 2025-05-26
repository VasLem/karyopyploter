from itertools import chain
from pathlib import Path

from matplotlib import pyplot as plt
from karyopyploter import plot_ideogram
from karyopyploter.ideogram import GENOME, DETAIL, ORIENTATION

OUT_DIR = Path(__file__).parent.parent / "example_outputs" / "ideogram"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def test_simple_vertical_chr3():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(3, 25),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr3", left_margin=0, label="", orientation=ORIENTATION.VERTICAL)

    fig.savefig(OUT_DIR / "testing_vert_chr3.png", bbox_inches="tight")
def test_simple_horizontal_chr3():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(25, 3),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr3", left_margin=0, label="", orientation=ORIENTATION.HORIZONTAL)

    fig.savefig(OUT_DIR / "testing_horz_chr3.png", bbox_inches="tight")
    
def test_simple_horizontal_chr1():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(25, 3),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr1", left_margin=0, label="", orientation=ORIENTATION.HORIZONTAL)

    fig.savefig(OUT_DIR / "testing_horz_chr1.png", bbox_inches="tight")
    

def test_simple_vertical_chr1_start_stop():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(3, 25),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr1", left_margin=0, label="", orientation=ORIENTATION.VERTICAL, start=150000, stop=50000000)

    fig.savefig(OUT_DIR / "testing_vert_chr1_start_stop.png", bbox_inches="tight")
    
    
def test_simple_horizontal_chr1_start_stop():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(25, 3),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr1", left_margin=0, label="", orientation=ORIENTATION.HORIZONTAL, start=150000, stop=50000000)

    fig.savefig(OUT_DIR / "testing_horz_start_stop.png", bbox_inches="tight")
    

def test_23_vertical_chm13():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=24, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax, target=f"chr{i}", label=f"Chr. {i}", 
            orientation=ORIENTATION.VERTICAL, genome=genome,
            grid_params=dict(left_margin=0),
            label_params={"rotation": 90}
        )

    fig.savefig(OUT_DIR / "testing_vert_23.png", bbox_inches="tight")


def test_23_vertical_hg38():
    genome = GENOME.HG38
    fig, axes = plt.subplots(ncols=22, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23)), strict=True):
        _ax = plot_ideogram(
            ax, target=f"chr{i}", label=f"Chr. {i}", left_margin=0, 
            label_params={"rotation": 90},
            orientation=ORIENTATION.VERTICAL, 
            genome=genome, relative=True
        )

    fig.savefig(OUT_DIR / "testing_vert_23_hg38.png", bbox_inches="tight")


def test_23_vertical_hg19():
    genome = GENOME.HG19
    fig, axes = plt.subplots(ncols=22, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23)), strict=True):
        _ax = plot_ideogram(
            ax, target=f"chr{i}", label=f"Chr. {i}", 
            left_margin=0, orientation=ORIENTATION.VERTICAL, 
            label_params={"rotation": 90},
            genome=genome
        )

    fig.savefig(OUT_DIR / "testing_vert_23_hg19.png", bbox_inches="tight")


def test_23_vertical_chm13_bare():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=24, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            left_margin=0,
            orientation=ORIENTATION.VERTICAL,
            label_params={"rotation": 90},
            cytobands=DETAIL.BARE,
            genome=genome,
        )

    fig.savefig(OUT_DIR / "testing_vert_bare_23.png", bbox_inches="tight")


def test_23_horizontal_chm13_bare():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=1, nrows=24, figsize=(25, 15), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            left_margin=0,
            orientation=ORIENTATION.HORIZONTAL,
            cytobands=DETAIL.BARE,
            genome=genome,
        )

    fig.savefig(OUT_DIR / "testing_horz_bare_23.png", bbox_inches="tight")


def test_23_vertical_chm13_regions():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=24, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    per_chr_regions = {"chr1": [(0, 1000000, "black"), (20_000_000, 35_000_000, "red")]}
    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        regions = per_chr_regions.get(f"chr{i}")
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            left_margin=0,
            orientation=ORIENTATION.VERTICAL,
            label_params={"rotation": 90},
            genome=genome,
            regions=regions,
        )

    fig.savefig(OUT_DIR / "testing_vert_23_regions_chm13.png", bbox_inches="tight")


def test_23_horz_chm13_regions():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=1, nrows=24, figsize=(15, 25), facecolor="white", sharey=True)

    per_chr_regions = {"chr1": [(0, 1000000, "black"), (20_000_000, 25_000_000, "red")]}
    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        regions = per_chr_regions.get(f"chr{i}", None)
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            y_label=f"Chr. {i}",
            left_margin=0,
            height=0.99,
            orientation=ORIENTATION.HORIZONTAL,
            genome=genome,
            regions=regions,
        )

    fig.savefig(OUT_DIR / "testing_horz_23_regions.png", bbox_inches="tight")