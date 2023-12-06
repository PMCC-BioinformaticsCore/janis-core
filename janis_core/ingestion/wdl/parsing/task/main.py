
import WDL 
from typing import Optional, Any

from janis_core import settings
from janis_core import CommandToolBuilder, ToolInput, ToolOutput, Selector

from .requirements import (
    TaskContainerParser, 
    TaskCpusParser, 
    TaskMemoryParser,
    TaskDiskParser, 
    TaskEnvVarsParser, 
    TaskFilesToCreateParser, 
    TaskDirsToCreateParser, 
)
from .io import TaskInputParser, TaskOutputParser
from .graph import get_decl_refs
from ..command import NativeSimpleParser, NativeArgumentParser, ShellCommandParser

def parse_task(task: WDL.Tree.Task) -> CommandToolBuilder:
    mark_ignored_inputs(task)
    cmdtool = CommandToolBuilder(
        tool=task.name,
        version='DEV',
        container='ubuntu:latest',
        base_command=None,
        inputs=[],
        outputs=[]
    )

    # requirements
    cmdtool._container = parse_container_requirement(task, cmdtool)
    cmdtool._cpus = parse_cpus_requirement(task, cmdtool)
    cmdtool._memory = parse_memory_requirement(task, cmdtool)
    cmdtool._disk = parse_disk_requirement(task, cmdtool)
    # cmdtool._env_vars = parse_env_vars(task, cmdtool)
    # cmdtool._files_to_create = parse_files_to_create(task, cmdtool)
    # cmdtool._directories_to_create = parse_dirs_to_create(task, cmdtool)
    
    # inputs / outputs / command
    cmdtool._inputs = parse_inputs(task, cmdtool)
    cmdtool._outputs = parse_outputs(task, cmdtool)
    cmdtool = parse_command(task, cmdtool)
    return cmdtool

def mark_ignored_inputs(task: WDL.Tree.Task) -> None:
    if task.inputs is None:
        return 
    decls = set([inp.name for inp in task.inputs])
    decl_refs = get_decl_refs(task, decls)
    ignore_inputs = set()
    for inp in task.inputs:
        if inp.name not in decl_refs:
            raise RuntimeError
        if decl_refs[inp.name]['inputs'] > 0:
            continue
        elif decl_refs[inp.name]['scopedvars'] > 0:
            continue
        elif decl_refs[inp.name]['command'] > 0:
            continue
        elif decl_refs[inp.name]['outputs'] > 0:
            continue
        else:
            ignore_inputs.add(inp.name) 
    task.ignored_inputs = ignore_inputs

def parse_container_requirement(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> str:
    parser = TaskContainerParser(task, cmdtool)
    return parser.parse()

def parse_cpus_requirement(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> Optional[int]:
    parser = TaskCpusParser(task, cmdtool)
    return parser.parse()

def parse_memory_requirement(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> Optional[float]:
    parser = TaskMemoryParser(task, cmdtool)
    return parser.parse()

def parse_disk_requirement(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> Optional[float]:
    parser = TaskDiskParser(task, cmdtool)
    return parser.parse()

# def parse_env_vars(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> dict[str, Any]:
#     parser = TaskEnvVarsParser(task, cmdtool)
#     return parser.parse()

# def parse_files_to_create(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> dict[str, Any]:
#     parser = TaskFilesToCreateParser(task, cmdtool)
#     return parser.parse()

# def parse_dirs_to_create(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> list[str | Selector]:
#     parser = TaskDirsToCreateParser(task, cmdtool)
#     return parser.parse()

def parse_inputs(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> list[ToolInput]:
    valid_inputs = get_valid_inputs(task)
    return [parse_input(task, cmdtool, wdl_inp) for wdl_inp in valid_inputs]

def get_valid_inputs(task: WDL.Tree.Task) -> Any:
    # no inputs 
    if task.inputs is None:
        return []
    assert task.ignored_inputs is not None
    normal_inputs = [decl for decl in task.inputs if decl.name not in task.ignored_inputs]
    post_inputs = [decl for decl in task.postinputs if isinstance(decl, WDL.Tree.Decl)]
    return normal_inputs + post_inputs
    # decls = set([inp.name for inp in task.inputs])
    # decl_refs = get_decl_refs(task, decls)
    # valid_inps = []
    # for inp in task.inputs:
    #     if inp.name not in decl_refs:
    #         raise RuntimeError
    #     if decl_refs[inp.name]['inputs'] > 0:
    #         valid_inps.append(inp)
    #     elif decl_refs[inp.name]['scopedvars'] > 0:
    #         valid_inps.append(inp)
    #     elif decl_refs[inp.name]['command'] > 0:
    #         valid_inps.append(inp)
    #     elif decl_refs[inp.name]['outputs'] > 0:
    #         valid_inps.append(inp)
    #     else:
    #         continue 
    # return valid_inps

def parse_input(task: WDL.Tree.Task, cmdtool: CommandToolBuilder, wdl_inp: WDL.Tree.Decl) -> ToolInput:
    parser = TaskInputParser(task, cmdtool, wdl_inp)
    return parser.parse()
    
def parse_outputs(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> list[ToolOutput]:
    # no outputs 
    if task.outputs is None:
        return []
    return [parse_output(task, cmdtool, wdl_out) for wdl_out in task.outputs]
    
def parse_output(task: WDL.Tree.Task, cmdtool: CommandToolBuilder, wdl_out: WDL.Tree.Decl) -> ToolOutput:
    parser = TaskOutputParser(task, cmdtool, wdl_out)
    return parser.parse()

def parse_command(task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> CommandToolBuilder:
    if task.command is None:
        return cmdtool
    
    if settings.ingest.wdl.COMMAND_PARSER == 'native_simple':
        p_classes = [NativeSimpleParser]
    elif settings.ingest.wdl.COMMAND_PARSER == 'native_arguments':
        p_classes = [NativeArgumentParser]
    elif settings.ingest.wdl.COMMAND_PARSER == 'shell':
        p_classes = [ShellCommandParser]
    else:
        p_classes = [NativeSimpleParser, NativeArgumentParser, ShellCommandParser]
    
    #try native approach
    for p_class in p_classes:
        sfmode = settings.ingest.SAFE_MODE
        settings.ingest.SAFE_MODE = True
        parser = p_class(task, cmdtool)
        parser.parse()
        if parser.success:
            cmdtool._base_command = parser.base_command
            cmdtool._env_vars = parser.env_vars
            cmdtool._files_to_create = parser.files_to_create
            cmdtool._directories_to_create = parser.directories_to_create
            return cmdtool
        settings.ingest.SAFE_MODE = sfmode

    # TODO error handling here
    raise RuntimeError

    