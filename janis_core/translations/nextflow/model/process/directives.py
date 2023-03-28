

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
    scope: Scope

    def get_string(self) -> str:
        scope = self.scope.labels[1:]  # remove 'settings.NF_MAIN_NAME' from start of the scope
        subpath = '/'.join(scope)
        subpath = to_case(subpath, settings.translate.nextflow.NF_OUTDIR_CASE)
        if subpath == '':
            return f"publishDir \"$params.outdir\""
        else:
            return f"publishDir \"${{params.outdir}}/{subpath}\""

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
        return f'debug {self.debug}'

@dataclass
class NFCpusDirective(NFProcessDirective):
    param: Param

    def get_string(self) -> str:
        return f'cpus \"${{params.{self.param.name}}}\"'

@dataclass
class NFDiskDirective(NFProcessDirective):
    param: Param
    
    def get_string(self) -> str:
        return f'disk \"${{params.{self.param.name}}}\"'

@dataclass
class NFMemoryDirective(NFProcessDirective):
    param: Param
    
    def get_string(self) -> str:
        return f'memory \"${{params.{self.param.name}}}\"'

@dataclass
class NFTimeDirective(NFProcessDirective):
    param: Param
    
    def get_string(self) -> str:
        return f'time \"${{params.{self.param.name}}}\"'