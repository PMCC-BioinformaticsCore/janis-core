from enum import Enum
from typing import Optional, Union, List, Dict


class DocumentationMeta:
    """
    Base class to provide doc tags with more information tags could be a
    """

    def __init__(self, doc: Optional[str]):
        self.doc = doc

    def __str__(self):
        return self.doc or ""


class InputQualityType(Enum):
    user = "user"
    static = "static"
    configuration = "configuration"


class InputDocumentation(DocumentationMeta):
    def __init__(
        self,
        doc: Optional[str],
        quality: InputQualityType = InputQualityType.user,
        example: Optional[Union[str, List[str]]] = None,
        source: Optional[
            Union[str, List[str], Dict[str, Union[str, List[str]]]]
        ] = None,
    ):
        """
        Extended documentation for inputs
        :param doc:
        :param quality:
        :param example:
        :param source:
        """
        super().__init__(doc)

        self.quality = quality
        self.example = example
        self.source = source

    @staticmethod
    def try_parse_from(doc: Union[str, Dict[str, str], any]):
        if doc is None or isinstance(doc, str):
            return InputDocumentation(doc=doc)
        elif isinstance(doc, InputDocumentation):
            return doc
        elif isinstance(doc, dict):
            return InputDocumentation(**doc)
        else:
            raise TypeError(
                f"Unexpected type when parsing InputDocumentation, expected "
                f"'Union[str, Dict[str, any], InputDocumentation]', received '{type(doc)}'."
            )


class OutputDocumentation(DocumentationMeta):
    def __init__(self, doc: Optional[str]):
        super().__init__(doc)
