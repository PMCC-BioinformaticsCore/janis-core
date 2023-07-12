

from __future__ import annotations
from typing import Optional, Protocol
from janis_core.ingestion.galaxy.gxtool.model import XMLParam
from janis_core.ingestion.galaxy import datatypes
from enum import Enum, auto

class ComponentConfidence(Enum):
    HIGH    = 3
    MEDIUM  = 2
    LOW     = 1


# defines what a CommandComponent should look like
class CommandComponent(Protocol):
    uuid: str
    gxparam: Optional[XMLParam]
    forced_optionality: Optional[bool]

    @property
    def datatype(self) -> datatypes.JanisDatatype:
        ...

    @property
    def name(self) -> str:
        ...

    @property
    def optional(self) -> bool:
        ...

    @property
    def array(self) -> bool:
        ...
    
    @property
    def docstring(self) -> Optional[str]:
        ...

    def set_confidence(self, level: str) -> None:
        """set the confidence level ('high'| 'medium' | 'low') of this component."""
        ...

