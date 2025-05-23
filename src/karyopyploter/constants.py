from enum import Enum
from pathlib import Path
import re

class GENOME(Enum):
    HG19 = "hg19"
    HG38 = "hg38"
    CHM13 = "chm13"
    HS1 = "hs1"


class ORIENTATION(Enum):
    VERTICAL = "Vertical"
    HORIZONTAL = "Horizontal"


class DETAIL(Enum):
    CYTOBAND = "Cytoband"
    BARE = "Bare"
    


COLOUR_LOOKUP = {
    "gneg": (1.0, 1.0, 1.0),
    "gpos25": (0.6, 0.6, 0.6),
    "gpos50": (0.4, 0.4, 0.4),
    "gpos75": (0.2, 0.2, 0.2),
    "gpos100": (0.0, 0.0, 0.0),
    # 'acen': (.8, .4, .4),
    # Set acen to be white as we use a
    #   polygon to add it in later
    "acen": (1.0, 1.0, 1.0),
    "gvar": (0.8, 0.8, 0.8),
    "stalk": (0.9, 0.9, 0.9),
}
STATIC_PATH = Path(__file__).parent / "static"

CHR_PATT = re.compile(r"^(?:chr([0-9]+|x|y|m)(.*)|(.*))$")
