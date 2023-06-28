

from __future__ import annotations
from typing import Optional, Protocol
from janis_core.ingestion.galaxy.gx.gxtool.param import Param


# defines what a CommandComponent should look like
class CommandComponent(Protocol):
    uuid: str
    gxparam: Optional[Param]
    forced_optionality: Optional[bool]

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


