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
        quality: Union[InputQualityType, str] = InputQualityType.user,
        example: Optional[Union[str, List[str]]] = None,
        source: Optional[
            Union[str, List[str], Dict[str, Union[str, List[str]]]]
        ] = None,
        skip_sourcing_secondary_files=False,
    ):
        """
        Extended documentation for inputs
        :param doc: Documentation string
        :type doc: str
        :param quality: quality of input, whether the inputs are best classified by user (data), static (references), configuration (like constants, but tweakable)
        :type quality: InputQualityType | "user" | "static" | "configuration"
        :param example: An example of the filename, displayed in the generated example input.yaml
        :type example: str | List[str]
        :param source: A URI of this input, that Janis could localise if it's not provided. For example, you might want to specify a gs://<path>
        :type source: str | List[str] | Dict[str, str | List[str]]
        :param skip_sourcing_secondary_files: Skip localising the secondary files from the source. You might want to do this if the secondary files depend on the version of the tool (eg: BWA)
        :type skip_sourcing_secondary_files: bool
        """
        super().__init__(doc)

        if quality is not None and not isinstance(quality, InputQualityType):
            quality = InputQualityType(quality)

        self.quality = quality
        self.example = example
        self.source = source
        self.skip_sourcing_secondary_files = skip_sourcing_secondary_files

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
