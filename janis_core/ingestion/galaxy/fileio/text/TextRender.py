
from abc import ABC, abstractmethod
from enum import Enum, auto
from typing import Tuple

from janis_core.ingestion.galaxy.messages import Message


class CurationState(Enum):
    RAW     = auto()
    EDITED  = auto()
    CURATED = auto()

class Confidence(Enum):
    LOW     = auto()
    MEDIUM  = auto()
    HIGH    = auto()


class TextRender(ABC):
    """
    A TextRender is a text representation of an entity (tool input, workflow step etc).
    Similar to __str__ or __repr__ dunder methods but not coupled to the class. 
    Includes other information like imports and messages related to that entity. 
    Used to create the output definitions users will receive. 
    \n
    Accounts for all the different ways the entity could be rendered. 
    Example 1: whether or not imports will be rendered here
    Example 2: how messages are 
    displayed to the user etc. 
    """

    def __init__(self):
        self.messages: list[Message] = []
        self.curation_state: CurationState = CurationState.RAW
        self.confidence: Confidence = Confidence.LOW

    @property
    @abstractmethod
    def imports(self) -> list[Tuple[str, str]]:
        """
        collects all imports this entity will need.
        each import should be returned as [path(module), object]. 
        """
        ...

    @abstractmethod
    def render(self) -> str:
        """renders an entity to text"""
        ...
    




