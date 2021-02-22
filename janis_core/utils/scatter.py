from enum import Enum
from typing import List, Union

from janis_core.operators.selectors import (
    Selector,
    InputNodeSelector,
    StepOutputSelector,
    InputSelector,
)
from janis_core.operators.operator import Operator


class ScatterMethod(Enum):
    """
    The scatter methods that Janis supports:
        - Dot: Inner product of two vectors of the same length, eg: Dot [A, B, C] · [1, 2, 3] => [A1, B2, C3]
        - Cross: Cartesian product of two vectors, producing every combination of the scattered inputs [A, B] x [1, 2] => [A1, A2, B1, B2]
    """

    dot = "dot"
    cross = "cross"

    def cwl(self):
        if self == ScatterMethod.dot:
            return "dotproduct"
        elif self == ScatterMethod.cross:
            return "flat_crossproduct"  # "nested_crossproduct"

        raise Exception(f"Unrecognised scatter method: '{self.value}'")


ScatterMethods = ScatterMethod


class ScatterDescription:
    """
    Class for keeping track of scatter information
    """

    def __init__(
        self,
        fields: List[str],
        method: ScatterMethod = None,
        labels: Union[Selector, List[str]] = None,
    ):
        """

        :param fields: The fields of the the tool that should be scattered on.
        :param method: The method that should be used to scatter the two arrays
        :param labels: (JANIS ONLY) -
        :type method: ScatterMethod
        """
        self.fields = fields
        self.method: ScatterMethod = method

        self.labels = None
        if labels is not None:
            if isinstance(labels, list):
                self.labels = map(str, labels)
            elif isinstance(labels, InputNodeSelector):
                labels = InputSelector(labels.id())
            elif isinstance(labels, StepOutputSelector):
                raise Exception(
                    f"Forbidden: Unable to use StepOutputSelector '{str(labels)}' for scatter label."
                )
            elif isinstance(labels, Operator):
                if any(isinstance(l, StepOutputSelector) for l in labels.get_leaves()):
                    raise Exception(
                        f"Forbidden: There was a StepOutputSelector as a parameter to the scatter "
                        f"label operator '{str(labels)}' which is not allowed."
                    )

            self.labels = labels

        if len(fields) > 1 and method is None:
            raise Exception(
                "When there is more than one field, a ScatterMethod must be selected"
            )
