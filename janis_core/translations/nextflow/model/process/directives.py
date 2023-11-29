

from abc import ABC, abstractmethod
from dataclasses import dataclass

from janis_core import settings

from ...casefmt import to_case
from ...params import Param
from ...scope import Scope


@dataclass
class NFProcessDirective(ABC):

    @abstractmethod
    def get_string(self) -> str:
        ...

@dataclass
class NFPublishDirDirective(NFProcessDirective):
    name: str

    def get_string(self) -> str:
        path = to_case(self.name, settings.translate.nextflow.NF_OUTDIR_CASE)
        if path == '':
            return f"publishDir \"$params.outdir\""
        else:
            return f"publishDir \"${{params.outdir}}/{path}\""

        #return f"publishDir \"$params.outdir/$task.process\""
        #return f"publishDir \"$launchDir/{self.process_name}\""
    
@dataclass
class NFCacheDirective(NFProcessDirective):
    enabled: bool

    def get_string(self) -> str:
        return f"cache {str(self.enabled).lower()}"

@dataclass
class NFContainerDirective(NFProcessDirective):
    container: str

    def get_string(self) -> str:
        return f'container "{self.container}"'

@dataclass
class NFDebugDirective(NFProcessDirective):
    debug: str

    def get_string(self) -> str:
        return f''
        # return f'debug {self.debug}'

@dataclass
class NFCpusDirective(NFProcessDirective):
    value: Param | str

    def get_string(self) -> str:
        if isinstance(self.value, Param):
            return f'cpus "${{params.{self.value.name}}}"'
        return f'cpus "{self.value}"'

@dataclass
class NFDiskDirective(NFProcessDirective):
    value: Param | str
    
    def get_string(self) -> str:
        if isinstance(self.value, Param):
            return f'disk "${{params.{self.value.name}}}"'
        return f'disk "{self.value}"'

@dataclass
class NFMemoryDirective(NFProcessDirective):
    value: Param | str
    
    def get_string(self) -> str:
        if isinstance(self.value, Param):
            return f'memory "${{params.{self.value.name}}}"'
        return f'memory "{self.value}"'

@dataclass
class NFTimeDirective(NFProcessDirective):
    value: Param | str
    
    def get_string(self) -> str:
        if isinstance(self.value, Param):
            return f'time "${{params.{self.value.name}}}"'
        return f'time "{self.value}"'