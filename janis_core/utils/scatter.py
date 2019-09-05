from dataclasses import dataclass
from enum import Enum
from typing import List

ScatterMethod = str


class ScatterMethods(Enum):
    dot: ScatterMethod = "dot"
    cross: ScatterMethod = "cross"


@dataclass(unsafe_hash=True)
class ScatterDescription:
    """
    Class for keeping track of scatter information
    """

    fields: List[str]
    method: ScatterMethod = None
