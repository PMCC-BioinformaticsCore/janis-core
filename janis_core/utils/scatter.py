from enum import Enum
from typing import List

ScatterMethod = str


class ScatterMethods(Enum):
    dot: ScatterMethod = "dot"
    cross: ScatterMethod = "cross"

    def cwl(self):
        if self == ScatterMethods.dot:
            return "dotproduct"
        elif self == ScatterMethods.cross:
            return "flat_crossproduct"  # "nested_crossproduct"

        raise Exception(f"Unrecognised scatter method: '{self.value}'")


class ScatterDescription:
    """
    Class for keeping track of scatter information
    """

    def __init__(self, fields: List[str], method: ScatterMethods = None):
        self.fields = fields
        self.method: ScatterMethods = method
        if len(fields) > 1 and method is None:
            raise Exception(
                "When there is more than one field, a ScatterMethod must be selected"
            )
