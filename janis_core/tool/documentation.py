from enum import Enum
from typing import Optional


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
        quality: InputQualityType = None,
        example: Optional[str] = None,
    ):
        super().__init__(doc)

        self.quality = quality
        self.example = example


class OutputDocumentation(DocumentationMeta):
    def __init__(self, doc: Optional[str]):
        super().__init__(doc)
