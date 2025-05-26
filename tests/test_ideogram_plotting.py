from itertools import chain
from pathlib import Path

from matplotlib import pyplot as plt
from karyopyploter import plot_ideogram
from karyopyploter.ideogram import GENOME, DETAIL, ORIENTATION
from filecmp import cmp as filecmp

TEST_DIR = Path(__file__).parent.parent / "test_data" / "ideogram"
TEST_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR = Path(__file__).parent.parent / "example_outputs" / "ideogram"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def test_simple_vertical_chr3():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(3, 25),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr3", label="", orientation=ORIENTATION.VERTICAL)

    fig.savefig(TEST_DIR / "vert_chr3.png")
    assert filecmp(TEST_DIR / "vert_chr3.png", OUT_DIR / "vert_chr3.png")


def test_simple_horizontal_chr3():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(25, 3),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr3", label="", orientation=ORIENTATION.HORIZONTAL)

    fig.savefig(TEST_DIR / "horz_chr3.png")
    assert filecmp(TEST_DIR / "horz_chr3.png", OUT_DIR / "horz_chr3.png")


def test_simple_horizontal_chr1():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(25, 3),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr1", label="", orientation=ORIENTATION.HORIZONTAL)

    fig.savefig(TEST_DIR / "horz_chr1.png")
    assert filecmp(TEST_DIR / "horz_chr1.png", OUT_DIR / "horz_chr1.png")


def test_simple_horizontal_chr1_with_coordinates():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(25, 3),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr1", label="", orientation=ORIENTATION.HORIZONTAL, show_coordinates=True)

    fig.savefig(TEST_DIR / "horz_chr1_with_coordinates.png")
    assert filecmp(TEST_DIR / "horz_chr1_with_coordinates.png", OUT_DIR / "horz_chr1_with_coordinates.png")


def test_simple_vertical_chr1_start_stop():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(3, 25),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr1", label="", orientation=ORIENTATION.VERTICAL, start=150000, stop=50000000)

    fig.savefig(TEST_DIR / "vert_chr1_start_stop.png")
    assert filecmp(
        TEST_DIR / "vert_chr1_start_stop.png",
        OUT_DIR / "vert_chr1_start_stop.png",
    )


def test_simple_vertical_chr1_start_stop_with_coordinates():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(3, 25),
        facecolor="white",
    )

    plot_ideogram(
        ax,
        target="chr1",
        label="",
        orientation=ORIENTATION.VERTICAL,
        start=150000,
        stop=50000000,
        show_coordinates=True,
    )

    fig.savefig(TEST_DIR / "vert_chr1_start_stop_with_coordinates.png")
    assert filecmp(
        TEST_DIR / "vert_chr1_start_stop_with_coordinates.png",
        OUT_DIR / "vert_chr1_start_stop_with_coordinates.png",
    )


def test_simple_horizontal_chr1_start_stop():
    fig, ax = plt.subplots(
        ncols=1,
        nrows=1,
        figsize=(25, 3),
        facecolor="white",
    )

    plot_ideogram(ax, target="chr1", label="", orientation=ORIENTATION.HORIZONTAL, start=150000, stop=50000000)

    fig.savefig(TEST_DIR / "horz_start_stop.png")
    assert filecmp(
        TEST_DIR / "horz_start_stop.png",
        OUT_DIR / "horz_start_stop.png",
    )


def test_23_vertical_chm13():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=24, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            orientation=ORIENTATION.VERTICAL,
            genome=genome,
            label_params={"rotation": 90},
        )

    fig.savefig(TEST_DIR / "vert_23.png")
    assert filecmp(
        TEST_DIR / "vert_23.png",
        OUT_DIR / "vert_23.png",
    )


def test_23_vertical_hg38():
    genome = GENOME.HG38
    fig, axes = plt.subplots(ncols=22, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23)), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            label_params={"rotation": 90},
            orientation=ORIENTATION.VERTICAL,
            genome=genome,
            relative=True,
        )

    fig.savefig(TEST_DIR / "vert_23_hg38.png")
    assert filecmp(
        TEST_DIR / "vert_23_hg38.png",
        OUT_DIR / "vert_23_hg38.png",
    )


def test_23_vertical_hg19():
    genome = GENOME.HG19
    fig, axes = plt.subplots(ncols=22, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23)), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            orientation=ORIENTATION.VERTICAL,
            label_params={"rotation": 90},
            genome=genome,
        )

    fig.savefig(TEST_DIR / "vert_23_hg19.png")
    assert filecmp(
        TEST_DIR / "vert_23_hg19.png",
        OUT_DIR / "vert_23_hg19.png",
    )


def test_23_vertical_chm13_bare():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=24, nrows=1, figsize=(15, 25), facecolor="white", sharey=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            orientation=ORIENTATION.VERTICAL,
            label_params={"rotation": 90},
            cytobands=DETAIL.BARE,
            genome=genome,
        )

    fig.savefig(TEST_DIR / "vert_bare_23.png")
    assert filecmp(
        TEST_DIR / "vert_bare_23.png",
        OUT_DIR / "vert_bare_23.png",
    )


def test_23_horizontal_chm13_bare():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=1, nrows=24, figsize=(25, 15), facecolor="white", sharex=True)

    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            orientation=ORIENTATION.HORIZONTAL,
            cytobands=DETAIL.BARE,
            genome=genome,
        )

    fig.savefig(TEST_DIR / "horz_bare_23.png")
    assert filecmp(
        TEST_DIR / "horz_bare_23.png",
        OUT_DIR / "horz_bare_23.png",
    )


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
            orientation=ORIENTATION.VERTICAL,
            label_params={"rotation": 90},
            genome=genome,
            regions=regions,
        )

    fig.savefig(TEST_DIR / "vert_23_regions_chm13.png")
    assert filecmp(
        TEST_DIR / "vert_23_regions_chm13.png",
        OUT_DIR / "vert_23_regions_chm13.png",
    )


def test_23_horz_chm13_regions():
    genome = GENOME.CHM13
    fig, axes = plt.subplots(ncols=1, nrows=24, figsize=(15, 25), facecolor="white", sharey=True)

    per_chr_regions = {"chr1": [(0, 1000000, "black"), (20_000_000, 25_000_000, "red")]}
    for ax, i in zip(axes, chain(range(1, 23), iter("XY")), strict=True):
        regions = per_chr_regions.get(f"chr{i}", None)
        _ax = plot_ideogram(
            ax,
            target=f"chr{i}",
            label=f"Chr. {i}",
            height=0.99,
            orientation=ORIENTATION.HORIZONTAL,
            genome=genome,
            regions=regions,
        )

    fig.savefig(TEST_DIR / "horz_23_regions.png")
    assert filecmp(
        TEST_DIR / "horz_23_regions.png",
        OUT_DIR / "horz_23_regions.png",
    )
