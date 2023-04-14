

from typing import Optional

from janis_core import (
    DataType,
    Workflow,
    CommandTool,
    PythonTool,
)

from ..casefmt import to_case
from ..scope import Scope
from janis_core import settings


### GENERAL

def get_construct_name(tool: CommandTool | PythonTool | Workflow, scope: Scope) -> str:
    construct_type = ''
    depth = len(scope.items)
    if isinstance(tool, CommandTool) or isinstance(tool, PythonTool):
        construct_type = 'process'
    elif isinstance(tool, Workflow) and depth == 1:  # scope = ['main']  (the main workflow)
        construct_type = 'main_workflow'
    elif isinstance(tool, Workflow) and depth > 1: # scope = ['main', 'sub', ...] 
        construct_type = 'sub_workflow'
    else:
        raise NotImplementedError
    return construct_type

def gen_varname_workflow(basename: str) -> str:
    return to_case(basename, settings.translate.nextflow.NF_PROCESS_CASE)

def gen_varname_process(basename: str) -> str:
    return to_case(basename, settings.translate.nextflow.NF_PROCESS_CASE)

def gen_varname_channel(janis_tag: str, name_override: Optional[str]=None, dtype: Optional[DataType]=None) -> str:
    basename = name_override if name_override else janis_tag
    # basename = _handle_plurals(basename, dtype)
    name = to_case(basename, settings.translate.nextflow.NF_CHANNEL_CASE)
    name = f'ch_{name}'
    return name

def gen_varname_file(janis_tag: str, name_override: Optional[str]=None, dtype: Optional[DataType]=None) -> str:
    basename = name_override if name_override else janis_tag
    # basename = _handle_plurals(basename, dtype)
    name = to_case(basename, settings.translate.nextflow.NF_CHANNEL_CASE)
    return name

def gen_varname_param(
    task_id: str, 
    tinput_id: Optional[str]=None, 
    name_override: Optional[str]=None, 
    is_subtask_param: bool=False
    ) -> str:
    assert(tinput_id or name_override)
    basename = name_override if name_override else tinput_id
    basename = to_case(basename, settings.translate.nextflow.NF_PARAM_CASE)
    if task_id and is_subtask_param:
        task_id = to_case(task_id, settings.translate.nextflow.NF_PARAM_CASE)
        name = f'{task_id}.{basename}'
    else:
        name = basename
    return name
