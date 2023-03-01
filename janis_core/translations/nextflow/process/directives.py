

from abc import ABC, abstractmethod
from typing import Any
from dataclasses import dataclass

from janis_core import CommandTool, PythonTool, Selector
from janis_core.types import Int

from ..casefmt import to_case
from ..params import Param
from ..scope import Scope
from janis_core import settings
from .. import params

# TODO: Create enums for relevant directives: https://www.nextflow.io/docs/latest/process.html#directives
# why? the module acts as an enum. currently can access directives via `directives.ProcessDirective`  etc - GH Dec 2022


@dataclass
class ProcessDirective(ABC):

    @abstractmethod
    def get_string(self) -> str:
        ...

@dataclass
class PublishDirDirective(ProcessDirective):
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
class CacheDirective(ProcessDirective):
    enabled: bool

    def get_string(self) -> str:
        return f"cache {str(self.enabled).lower()}"

@dataclass
class ContainerDirective(ProcessDirective):
    container: str

    def get_string(self) -> str:
        return f'container "{self.container}"'

@dataclass
class DebugDirective(ProcessDirective):
    debug: str

    def get_string(self) -> str:
        return f'debug {self.debug}'

@dataclass
class CpusDirective(ProcessDirective):
    param: Param

    def get_string(self) -> str:
        return f'cpus \"${{params.{self.param.name}}}\"'

@dataclass
class DiskDirective(ProcessDirective):
    param: Param
    
    def get_string(self) -> str:
        return f'disk \"${{params.{self.param.name}}}\"'

@dataclass
class MemoryDirective(ProcessDirective):
    param: Param
    
    def get_string(self) -> str:
        return f'memory \"${{params.{self.param.name}}}\"'

@dataclass
class TimeDirective(ProcessDirective):
    param: Param
    
    def get_string(self) -> str:
        return f'time \"${{params.{self.param.name}}}\"'



def gen_directives_for_process(tool: CommandTool | PythonTool, resources: dict[str, Any], scope: Scope) -> list[ProcessDirective]:
    
    # TODO REFACTOR
    nf_directives: dict[str, ProcessDirective] = {}
    nf_directives['publishDir'] = PublishDirDirective(scope)
    nf_directives['debug'] = DebugDirective(debug='true')

    # Add directives from input resources
    for res, val in resources.items():
        if res.endswith("runtime_cpu"):
            param = params.add(janis_tag='cpus', scope=scope, default=val, janis_dtype=Int())
            nf_directives['cpus'] = CpusDirective(param)
        
        elif res.endswith("runtime_memory"):
            param = params.add(janis_tag='memory', scope=scope, default=val, janis_dtype=Int())
            nf_directives['memory'] = MemoryDirective(param)
        
        elif res.endswith("runtime_seconds"):
            param = params.add(janis_tag='time', scope=scope, default=val, janis_dtype=Int())
            nf_directives['time'] = TimeDirective(param)
        
        elif res.endswith("runtime_disk"):
            param = params.add(janis_tag='disk', scope=scope, default=val, janis_dtype=Int())
            nf_directives['disk'] = DiskDirective(param)
    
    # Add directives from tool resources
    if 'cpus' not in nf_directives and tool.cpus({}) is not None:    
        param = params.add(janis_tag='cpus', scope=scope, default=tool.cpus({}), janis_dtype=Int())
        nf_directives['cpus'] = CpusDirective(param)
    
    if 'memory' not in nf_directives and tool.memory({}) is not None:
        param = params.add(janis_tag='memory', scope=scope, default=tool.memory({}), janis_dtype=Int())
        nf_directives['memory'] = MemoryDirective(param)
    
    if 'disk' not in nf_directives and tool.disk({}) is not None:
        param = params.add(janis_tag='disk', scope=scope, default=tool.disk({}), janis_dtype=Int())
        nf_directives['disk'] = DiskDirective(param)
    
    if 'time' not in nf_directives and tool.time({}) is not None:
        param = params.add(janis_tag='time', scope=scope, default=tool.time({}), janis_dtype=Int())
        nf_directives['time'] = TimeDirective(param)
    
    final_directives: list[ProcessDirective] = []
    for direc in nf_directives.values():
        if hasattr(direc, 'default') and isinstance(direc.default, Selector):
            continue
        else:
            final_directives.append(direc)
    return final_directives




# not implemented

# class AcceleratorDirective(ProcessDirective):
#     def __init__(self, accelerator):
#         super().__init__("accelerator", accelerator)

# class AfterScriptDirective(ProcessDirective):
#     def __init__(self, afterScript):
#         super().__init__("afterScript", afterScript)

# class BeforeScriptDirective(ProcessDirective):
#     def __init__(self, beforeScript):
#         super().__init__("beforeScript", beforeScript)

# class CondaDirective(ProcessDirective):
#     def __init__(self, conda):
#         super().__init__("conda", conda)

# class ContainerOptionsDirective(ProcessDirective):
#     def __init__(self, containerOptions):
#         super().__init__("containerOptions", containerOptions)

# class ClusterOptionsDirective(ProcessDirective):
#     def __init__(self, clusterOptions):
#         super().__init__("clusterOptions", clusterOptions)

# class EchoDirective(ProcessDirective):
#     def __init__(self, echo):
#         super().__init__("echo", echo)

# class ErrorStrategyDirective(ProcessDirective):
#     def __init__(self, errorStrategy):
#         super().__init__("errorStrategy", errorStrategy)

# class ExecutorDirective(ProcessDirective):
#     def __init__(self, executor):
#         super().__init__("executor", executor)

# class ExtDirective(ProcessDirective):
#     def __init__(self, ext):
#         super().__init__("ext", ext)

# class LabelDirective(ProcessDirective):
#     def __init__(self, label):
#         super().__init__("label", label)

# class MachineTypeDirective(ProcessDirective):
#     def __init__(self, machineType):
#         super().__init__("machineType", machineType)

# class MaxErrorsDirective(ProcessDirective):
#     def __init__(self, maxErrors):
#         super().__init__("maxErrors", maxErrors)

# class MaxForksDirective(ProcessDirective):
#     def __init__(self, maxForks):
#         super().__init__("maxForks", maxForks)

# class MaxRetriesDirective(ProcessDirective):
#     def __init__(self, maxRetries):
#         super().__init__("maxRetries", maxRetries)

# class ModuleDirective(ProcessDirective):
#     def __init__(self, module):
#         super().__init__("module", module)

# class PenvDirective(ProcessDirective):
#     def __init__(self, penv):
#         super().__init__("penv", penv)

# class PodDirective(ProcessDirective):
#     def __init__(self, pod):
#         super().__init__("pod", pod)

# class QueueDirective(ProcessDirective):
#     def __init__(self, queue):
#         super().__init__("queue", queue)

# class ScratchDirective(ProcessDirective):
#     def __init__(self, scratch):
#         super().__init__("scratch", scratch)

# class StageInModeDirective(ProcessDirective):
#     def __init__(self, stageInMode):
#         super().__init__("stageInMode", stageInMode)

# class StageOutModeDirective(ProcessDirective):
#     def __init__(self, stageOutMode):
#         super().__init__("stageOutMode", stageOutMode)

# class StoreDirDirective(ProcessDirective):
#     def __init__(self, storeDir):
#         super().__init__("storeDir", storeDir)

# class TagDirective(ProcessDirective):
#     def __init__(self, tag):
#         super().__init__("tag", tag)

# class ValidExitStatusDirective(ProcessDirective):
#     def __init__(self, validExitStatus):
#         super().__init__("validExitStatus", validExitStatus)

