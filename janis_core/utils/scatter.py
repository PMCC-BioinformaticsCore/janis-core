from enum import Enum
from typing import List

ScatterMethod = str


class ScatterMethods(Enum):
    """
    The scatter methods that Janis supports:
        - Dot: Inner product of two vectors of the same length, eg: Dot [A, B, C] · [1, 2, 3] => [A1, B2, C3]
        - Cross: Cartesian product of two vectors, producing every combination of the scattered inputs [A, B] x [1, 2] => [A1, A2, B1, B2]
    """

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
        """

        :param fields: The fields of the the tool that should be scattered on.
        :param method: The method that should be used to scatter the two arrays
        :type method: ScatterMethods
        """
        self.fields = fields
        self.method: ScatterMethods = method
        if len(fields) > 1 and method is None:
            raise Exception(
                "When there is more than one field, a ScatterMethod must be selected"
            )
