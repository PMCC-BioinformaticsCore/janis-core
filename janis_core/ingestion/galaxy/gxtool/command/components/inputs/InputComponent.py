
from abc import ABC, abstractmethod
from typing import Any, Optional
from janis_core.ingestion.galaxy.gxtool.model import XMLParam
from uuid import uuid4

from janis_core.ingestion.galaxy import tags
from janis_core.ingestion.galaxy import datatypes
from ...components import ComponentConfidence



class InputComponent(ABC):
    def __init__(self):
        self.uuid: str = str(uuid4())
        self.gxparam: Optional[XMLParam] = None
        self.cmd_pos: int = -1
        self.forced_optionality: Optional[bool] = None
        self.forced_array: Optional[bool] = None
        self.forced_default: Optional[bool] = None
        self.confidence: ComponentConfidence = ComponentConfidence.LOW
 
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
    def default_value(self) -> Any:
        """
        gets the default value of the component.
        uses the components observed values (positionals / options) and
        galaxy param information if available.
        """
        ...

    @property
    @abstractmethod
    def optional(self) -> bool:
        """
        returns whether the component is optional or not.
        uses galaxy param information if available, otherwise uses the presence array. flags are always optional
        """
        ...

    @property
    @abstractmethod
    def array(self) -> bool:
        """
        returns whether the component is an array or not
        uses galaxy param information if available.
        flags components are never arrays.
        """
        ...

    @property
    @abstractmethod
    def docstring(self) -> Optional[str]:
        """
        gets helptext for the component. uses galaxy param if available,
        otherwise usually just presents the witnessed values as examples. 
        """
        ...

    @abstractmethod
    def update(self, incoming: Any) -> None:
        """
        updates this component with information from another similar component. 
        also sets the component as being present in the cmdstr being parsed.
        example: 
        str1: abricate $input.fastq --db card > out.txt
        str2: abricate $input.fastq --db resfinder > out.txt
        the --db option should be updated so we know its possible values
        are (at least) card, resfinder
        """
        ...

    def set_confidence(self, level: str) -> None:
        """set the confidence level ('high'| 'medium' | 'low') of this component."""
        if level == 'high':
            self.confidence = ComponentConfidence.HIGH
        elif level == 'medium':
            self.confidence = ComponentConfidence.MEDIUM
        elif level == 'low':
            self.confidence = ComponentConfidence.LOW
        else:
            raise ValueError(f'invalid confidence level: {level}')


