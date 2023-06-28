

from dataclasses import dataclass
from typing import Any, Optional
from uuid import uuid4

from janis_core.ingestion.galaxy.datatypes import JanisDatatype
from janis_core.ingestion.galaxy import tags


@dataclass
class WorkflowInput:
    _name: str
    array: bool
    optional: bool
    is_runtime: bool
    datatype: JanisDatatype
    value: Any = None

    def __post_init__(self):
        self.uuid: str = str(uuid4())

    @property
    def name(self) -> str:
        if not self.is_runtime:
            if not self._name.startswith('in'):
                return f'in_{self._name}'
        return self._name
    
    @property
    def tag(self) -> str:
        return tags.get(self.uuid)
    
    @property
    def docstring(self) -> Optional[str]:
        return 'None yet!'
        

