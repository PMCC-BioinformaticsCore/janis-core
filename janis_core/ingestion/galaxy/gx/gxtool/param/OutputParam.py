
from typing import Any, Optional
from .Param import Param


class OutputParam(Param):
    def __init__(self, name: str):
        self.name: str = name
        self.label: str = ''
        self.formats: list[str] = []
        self.discover_pattern: Optional[str] = None

    @property
    def default(self) -> Any:
        return None

    @property
    def docstring(self) -> str:
        if self.label:
            return str(self.label)
        return ''
    
    @property
    def optional(self) -> bool:
        return False
    
    @property
    def array(self) -> bool:
        return False


class DataOutputParam(OutputParam):
    def __init__(self, name: str, from_work_dir: Optional[str]=None, discover_pattern: Optional[str]=None):
        super().__init__(name)
        self.from_work_dir = from_work_dir
        self.discover_pattern = discover_pattern

    @property
    def array(self) -> bool:
        if self.discover_pattern and '*' in self.discover_pattern:
            return True
        return False


class CollectionOutputParam(OutputParam):
    def __init__(self, name: str, discover_pattern: Optional[str]=None):
        super().__init__(name)
        self.discover_pattern = discover_pattern
        self.collection_type: str = 'list'

    @property
    def array(self) -> bool:
        return True

