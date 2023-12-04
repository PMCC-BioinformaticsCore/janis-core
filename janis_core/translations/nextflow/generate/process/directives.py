

from typing import Any

from janis_core import CommandToolBuilder, PythonTool, Selector, StringFormatter
from janis_core.types import Int
from janis_core import settings
from janis_core.translations.nextflow.unwrap import unwrap_expression

from ... import params
from ...variables import VariableManager

from ...model.process.directives import (
    NFProcessDirective,
    NFContainerDirective,
    NFPublishDirDirective,
    NFDebugDirective,
    NFCpusDirective,
    NFMemoryDirective,
    NFTimeDirective,
    NFDiskDirective
)


# TODO: Create enums for relevant directives: https://www.nextflow.io/docs/latest/process.html#directives
# why? the module acts as an enum. currently can access directives via `directives.ProcessDirective`  etc - GH Dec 2022


def gen_nf_process_directives(tool: CommandToolBuilder | PythonTool, vmanager: VariableManager) -> list[NFProcessDirective]:
    nf_directives: dict[str, NFProcessDirective] = {}
    
    # publishDir
    if settings.translate.nextflow.ENTITY == 'workflow':
        nf_directives['publishDir'] = NFPublishDirDirective(tool.id())
    
    # debug
    nf_directives['debug'] = NFDebugDirective(debug='true')

    # container
    uri = None
    if settings.translate.WITH_CONTAINER and tool.container():
        uri = tool.container()
        nf_directives['container'] = NFContainerDirective(uri)

    if uri is None and not settings.translate.ALLOW_EMPTY_CONTAINER:
        raise Exception(
            f"The tool '{tool.id()}' did not have a container and no container override was specified. "
            f"Although not recommended, Janis can export empty docker containers with the parameter "
            f"'allow_empty_container=True' or --allow-empty-container"
        )

    # # directives from input resources
    # for res, val in resources.items():
    #     if res.endswith("runtime_cpu"):
    #         param = params.add(tinput_id='cpus', task_id=tool.id(), default=val, janis_dtype=Int(), subtype='sub_tool')
    #         nf_directives['cpus'] = NFCpusDirective(param)
        
    #     elif res.endswith("runtime_memory"):
    #         param = params.add(tinput_id='memory', task_id=tool.id(), default=val, janis_dtype=Int(), subtype='sub_tool')
    #         nf_directives['memory'] = NFMemoryDirective(param)
        
    #     elif res.endswith("runtime_seconds"):
    #         param = params.add(tinput_id='time', task_id=tool.id(), default=val, janis_dtype=Int(), subtype='sub_tool')
    #         nf_directives['time'] = NFTimeDirective(param)
        
    #     elif res.endswith("runtime_disk"):
    #         param = params.add(tinput_id='disk', task_id=tool.id(), default=val, janis_dtype=Int(), subtype='sub_tool')
    #         nf_directives['disk'] = NFDiskDirective(param)
    
    # directives from tool resources
    if settings.translate.nextflow.ENTITY == 'workflow':
        if 'cpus' not in nf_directives and tool.cpus({}) is not None: 
            cpus = tool.cpus({})
            if isinstance(cpus, Selector) and not settings.ingest.SOURCE == 'wdl':
                res = unwrap_expression(val=cpus, context='process_script', variable_manager=vmanager, tool=tool, apply_braces=False)
                if isinstance(cpus, StringFormatter):
                    res = f'{res}'
                else:
                    res = f'${{{res}}}'
                nf_directives['cpus'] = NFCpusDirective(res)
            else:
                param = params.add(tinput_id='cpus', task_id=tool.id(), default=tool.cpus({}), janis_dtype=Int(), subtype='sub_tool')
                nf_directives['cpus'] = NFCpusDirective(param)
        
        if 'memory' not in nf_directives and tool.memory({}) is not None:
            memory = tool.memory({})
            if isinstance(memory, Selector) and not settings.ingest.SOURCE == 'wdl':
                res = unwrap_expression(val=memory, context='process_script', variable_manager=vmanager, tool=tool, apply_braces=False)
                if isinstance(memory, StringFormatter):
                    res = f'{res} MB'
                else:
                    res = f'${{{res}}} MB'
                nf_directives['memory'] = NFMemoryDirective(res)
            else:
                param = params.add(tinput_id='memory', task_id=tool.id(), default=tool.memory({}), janis_dtype=Int(), subtype='sub_tool')
                nf_directives['memory'] = NFMemoryDirective(param)
        
        if 'disk' not in nf_directives and tool.disk({}) is not None:
            disk = tool.disk({})
            if isinstance(disk, Selector) and not settings.ingest.SOURCE == 'wdl':  # bandaid fix for wdl
                res = unwrap_expression(val=disk, context='process_script', variable_manager=vmanager, tool=tool, apply_braces=False)
                if isinstance(disk, StringFormatter):
                    res = f'{res} MB'
                else:
                    res = f'${{{res}}} MB'
                nf_directives['disk'] = NFDiskDirective(res)
            else:
                param = params.add(tinput_id='disk', task_id=tool.id(), default=tool.disk({}), janis_dtype=Int(), subtype='sub_tool')
                nf_directives['disk'] = NFDiskDirective(param)
        
        if 'time' not in nf_directives and tool.time({}) is not None:
            time = tool.time({})
            if isinstance(time, Selector) and not settings.ingest.SOURCE == 'wdl':
                res = unwrap_expression(val=time, context='process_script', variable_manager=vmanager, tool=tool, apply_braces=False)
                if isinstance(time, StringFormatter):
                    res = f'{res}'
                else:
                    res = f'${{{res}}}'
                nf_directives['time'] = NFTimeDirective(res)
            else:
                param = params.add(tinput_id='time', task_id=tool.id(), default=tool.time({}), janis_dtype=Int(), subtype='sub_tool')
                nf_directives['time'] = NFTimeDirective(param)
    
    final_directives: list[NFProcessDirective] = []
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

