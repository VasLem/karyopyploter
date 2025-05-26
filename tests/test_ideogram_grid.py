from pathlib import Path

from karyopyploter import make_ideogram_grid
from filecmp import cmp as filecmp

TEST_DIR = Path(__file__).parent.parent / "test_data" / "ideogram_grid"
TEST_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR = Path(__file__).parent.parent / "example_outputs" / "ideogram_grid"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def test_ideogram_grid_generation():
    fig, axes, ideogram_ax = make_ideogram_grid(
        target="chr1",
        num_subplots=2,
    )
    fig.savefig(TEST_DIR / "chr1.png", bbox_inches="tight")
    assert filecmp(
        TEST_DIR / "chr1.png",
        OUT_DIR / "chr1.png",
    ), "Generated ideogram grid does not match expected output."


def test_ideogram_grid_generation_with_coordinates():
    fig, axes, ideogram_ax = make_ideogram_grid(
        target="chr1", num_subplots=2, ideogram_params=dict(show_coordinates=True)
    )
    fig.savefig(TEST_DIR / "chr1_with_coordinates.png", dpi=300)
    assert filecmp(
        TEST_DIR / "chr1_with_coordinates.png",
        OUT_DIR / "chr1_with_coordinates.png",
    ), "Generated ideogram grid with coordinates does not match expected output."


def test_ideogram_grid_generation_three_targets():
    fig, axes, ideogram_ax = make_ideogram_grid(
        target=["chr1", "chr2", "chr22"],
        num_subplots=2,
    )
    fig.savefig(TEST_DIR / "three_targets.png", bbox_inches="tight")
    assert filecmp(
        TEST_DIR / "three_targets.png",
        OUT_DIR / "three_targets.png",
    ), "Generated ideogram grid with three targets does not match expected output."


def test_ideogram_grid_generation_four_targets_two_columns():
    fig, axes, ideogram_ax = make_ideogram_grid(
        target=["chr1", "chr2", "chr5"], num_subplots=2, grid_params=dict(ncols=2)
    )
    fig.savefig(TEST_DIR / "three_targets_two_cols.png", bbox_inches="tight")
    assert filecmp(
        TEST_DIR / "three_targets_two_cols.png",
        OUT_DIR / "three_targets_two_cols.png",
    ), "Generated ideogram grid with three targets in two columns does not match expected output."


def test_ideogram_grid_generation_three_targets_with_start_stop():
    # should raise
    try:
        fig, axes, ideogram_ax = make_ideogram_grid(
            target=["chr1", "chr2", "chr22"],
            start=0,
            stop=1,
            num_subplots=2,
        )
        raise ValueError("Should have raised an error")
    except ValueError as e:
        pass
    fig, axes, ideogram_ax = make_ideogram_grid(
        target=["chr1", "chr2", "chr22"],
        start={"chr1": 0, "chr2": 0, "chr22": 0},
        stop={"chr1": 50000000, "chr2": 25000000, "chr22": 1500000},
        num_subplots=2,
    )
    fig.savefig(TEST_DIR / "three_targets_with_start_stop.png", bbox_inches="tight")
    assert filecmp(
        TEST_DIR / "three_targets_with_start_stop.png",
        OUT_DIR / "three_targets_with_start_stop.png",
    ), "Generated ideogram grid with three targets and start/stop does not match expected output."


def test_ideogram_grid_generation_three_targets_no_subplots():
    fig, axes, ideogram_ax = make_ideogram_grid(
        target=["chr1", "chr2", "chr22"],
        num_subplots=0,
    )
    fig.savefig(TEST_DIR / "three_targets_no_subplots.png", bbox_inches="tight")
    assert filecmp(
        TEST_DIR / "three_targets_no_subplots.png",
        OUT_DIR / "three_targets_no_subplots.png",
    ), "Generated ideogram grid with three targets and no subplots does not match expected output."
