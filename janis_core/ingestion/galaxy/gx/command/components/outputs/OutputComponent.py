



from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional

from janis_core.ingestion.galaxy import tags
from janis_core.ingestion.galaxy import datatypes

from uuid import uuid4
from janis_core.ingestion.galaxy.gx.gxtool.param import Param


class OutputComponent(ABC):
    def __init__(self):
        self.uuid: str = str(uuid4())
        self.gxparam: Optional[Param] = None
        self.forced_optionality: Optional[bool] = None
        self.forced_array: Optional[bool] = None
 
    @property
    @abstractmethod
    def name(self) -> str:
        """
        returns a name for this component. created depending on what
        information is available to the component
        """
        ...
    
    @property
    def tag(self) -> str:
        return tags.get(self.uuid)
    
    @property
    def datatype(self) -> datatypes.JanisDatatype:
        return datatypes.get(self)

    @property
    @abstractmethod
    def optional(self) -> bool:
        """
        returns whether the component is optional or not.
        uses galaxy param information if available, otherwise uses the presence array. flags are always optional
        """
        ...

    @property
    def array(self) -> bool:
        if self.forced_array is not None:
            return self.forced_array
        elif self.gxparam:
            return self.gxparam.array
        return False

    @property
    @abstractmethod
    def docstring(self) -> Optional[str]:
        """
        gets helptext for the component. uses galaxy param if available,
        otherwise usually just presents the witnessed values as examples. 
        """
        ...



