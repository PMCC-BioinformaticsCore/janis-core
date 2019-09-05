from enum import Enum
from typing import List

ScatterMethod = str


class ScatterMethods(Enum):
    dot: ScatterMethod = "dot"
    cross: ScatterMethod = "cross"


class ScatterDescription:
    """
    Class for keeping track of scatter information
    """

    def __init__(self, fields: List[str], method: ScatterMethod = None):
        self.fields = fields
        self.method = method
        if len(fields) > 1 and method is None:
            raise Exception(
                "When there is more than one field, a ScatterMethod must be selected"
            )
