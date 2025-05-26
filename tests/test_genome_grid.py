from pathlib import Path

from karyopyploter import make_genome_grid
from filecmp import cmp as filecmp

TEST_DIR = Path(__file__).parent.parent / "test_data" / "genome_grid"
TEST_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR = Path(__file__).parent.parent / "example_outputs" / "genome_grid"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def test_genome_grid_chr1to5():
    fig, axes, genome_ax = make_genome_grid(
        target_start="chr1",
        target_stop="chr5",
        num_subplots=2,
    )
    fig.savefig(TEST_DIR / "chr1_to_chr5.png", bbox_inches="tight")
    assert filecmp(
        TEST_DIR / "chr1_to_chr5.png",
        OUT_DIR / "chr1_to_chr5.png",
    ), "Generated genome grid does not match expected output."


def test_genome_grid_with_coordinates():
    fig, axes, genome_ax = make_genome_grid(
        target_start="chr1",
        target_stop="chr5",
        num_subplots=2,
        ideogram_params=dict(show_coordinates=True, label_axis_offset=1.2, coordinates_params=dict(rotation=90)),
    )
    fig.savefig(TEST_DIR / "chr1_to_chr5_with_coordinates.png", dpi=300)
    assert filecmp(
        TEST_DIR / "chr1_to_chr5_with_coordinates.png",
        OUT_DIR / "chr1_to_chr5_with_coordinates.png",
    ), "Generated genome grid with coordinates does not match expected output."
