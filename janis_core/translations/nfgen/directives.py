
from abc import ABC, abstractmethod

from janis_core.translations.nfgen.common import NFBase


# TODO: Create enums for relevant directives: https://www.nextflow.io/docs/latest/process.html#directives
# why? the module acts as an enum. currently can access directives via `directives.ProcessDirective`  etc

# implemented
# class ProcessDirective(NFBase, ABC):
#     def __init__(self, name: str, value: Any):
#         self.name = name
#         self.value = value

#     def get_string(self) -> str:
#         return f'{self.name} "{self.value}"'


class ProcessDirective(NFBase, ABC):
    @abstractmethod
    def get_string(self) -> str:
        ...

class CacheDirective(ProcessDirective):
    def __init__(self, enabled: bool):
        self.enabled = enabled

    def get_string(self) -> str:
        return f"cache {str(self.enabled).lower()}"


class CpusDirective(ProcessDirective):
    def __init__(self, varname: str):
        self.var = varname

    def get_string_old(self) -> str:
        return f"cpus \"${{{self.var} ? {self.var} : ''}}\""
    
    def get_string(self) -> str:
        return f"cpus \"${{{self.var}}}\""


class ContainerDirective(ProcessDirective):
    def __init__(self, container: str):
        self.container = container

    def get_string(self) -> str:
        return f'container "{self.container}"'


class DiskDirective(ProcessDirective):
    def __init__(self, varname: str):
        self.var = varname
    
    def get_string_old1(self) -> str:
        return f"disk \"${{{self.var} ? {self.var} + 'GB': ''}}\""
    
    def get_string_old2(self) -> str:
        return f"disk {self.var} ? \"{{{self.var} + 'GB'}}\" : null"
    
    def get_string(self) -> str:
        return f"disk \"${{{self.var}}}\""


class MemoryDirective(ProcessDirective):
    def __init__(self, varname: str):
        self.var = varname
    
    def get_string_old(self) -> str:
        return f"memory \"${{{self.var} ? {self.var} + 'GB': ''}}\""
    
    def get_string(self) -> str:
        return f"memory \"${{{self.var}}}\""


class PublishDirDirective(ProcessDirective):
    def __init__(self, process_name: str):
        self.pname = process_name

    def get_string(self) -> str:
        return f"publishDir \"$launchDir/{self.pname}\""


class TimeDirective(ProcessDirective):
    def __init__(self, varname: str):
        self.var = varname
    
    def get_string_old(self) -> str:
        return f"time \"${{{self.var} + 's'}}\""
    
    def get_string(self) -> str:
        return f"time \"${{{self.var}}}\""



# not implemented

class AcceleratorDirective(ProcessDirective):
    def __init__(self, accelerator):
        super().__init__("accelerator", accelerator)

class AfterScriptDirective(ProcessDirective):
    def __init__(self, afterScript):
        super().__init__("afterScript", afterScript)

class BeforeScriptDirective(ProcessDirective):
    def __init__(self, beforeScript):
        super().__init__("beforeScript", beforeScript)

class CondaDirective(ProcessDirective):
    def __init__(self, conda):
        super().__init__("conda", conda)

class ContainerOptionsDirective(ProcessDirective):
    def __init__(self, containerOptions):
        super().__init__("containerOptions", containerOptions)

class ClusterOptionsDirective(ProcessDirective):
    def __init__(self, clusterOptions):
        super().__init__("clusterOptions", clusterOptions)

class EchoDirective(ProcessDirective):
    def __init__(self, echo):
        super().__init__("echo", echo)

class ErrorStrategyDirective(ProcessDirective):
    def __init__(self, errorStrategy):
        super().__init__("errorStrategy", errorStrategy)

class ExecutorDirective(ProcessDirective):
    def __init__(self, executor):
        super().__init__("executor", executor)

class ExtDirective(ProcessDirective):
    def __init__(self, ext):
        super().__init__("ext", ext)

class LabelDirective(ProcessDirective):
    def __init__(self, label):
        super().__init__("label", label)

class MachineTypeDirective(ProcessDirective):
    def __init__(self, machineType):
        super().__init__("machineType", machineType)

class MaxErrorsDirective(ProcessDirective):
    def __init__(self, maxErrors):
        super().__init__("maxErrors", maxErrors)

class MaxForksDirective(ProcessDirective):
    def __init__(self, maxForks):
        super().__init__("maxForks", maxForks)

class MaxRetriesDirective(ProcessDirective):
    def __init__(self, maxRetries):
        super().__init__("maxRetries", maxRetries)

class ModuleDirective(ProcessDirective):
    def __init__(self, module):
        super().__init__("module", module)

class PenvDirective(ProcessDirective):
    def __init__(self, penv):
        super().__init__("penv", penv)

class PodDirective(ProcessDirective):
    def __init__(self, pod):
        super().__init__("pod", pod)

class QueueDirective(ProcessDirective):
    def __init__(self, queue):
        super().__init__("queue", queue)

class ScratchDirective(ProcessDirective):
    def __init__(self, scratch):
        super().__init__("scratch", scratch)

class StageInModeDirective(ProcessDirective):
    def __init__(self, stageInMode):
        super().__init__("stageInMode", stageInMode)

class StageOutModeDirective(ProcessDirective):
    def __init__(self, stageOutMode):
        super().__init__("stageOutMode", stageOutMode)

class StoreDirDirective(ProcessDirective):
    def __init__(self, storeDir):
        super().__init__("storeDir", storeDir)

class TagDirective(ProcessDirective):
    def __init__(self, tag):
        super().__init__("tag", tag)

class ValidExitStatusDirective(ProcessDirective):
    def __init__(self, validExitStatus):
        super().__init__("validExitStatus", validExitStatus)
