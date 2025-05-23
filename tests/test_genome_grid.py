from pathlib import Path

from karyopyploter import make_genome_grid

OUT_DIR = Path(__file__).parent.parent / "example_outputs" / "genome_grid"
OUT_DIR.mkdir(parents=True, exist_ok=True)

def test_genome_grid():
    fig, axes, genome_ax = make_genome_grid(
        target_start="chr1",
        target_stop="chr5",
        num_subplots=2,
    )
    fig.savefig(OUT_DIR / "chr1_to_chr5.png", bbox_inches="tight")