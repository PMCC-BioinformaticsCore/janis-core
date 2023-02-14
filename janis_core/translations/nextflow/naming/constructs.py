

from typing import Optional

from janis_core import (
    DataType,
    Workflow,
    CommandTool,
    PythonTool,
)

from ..casefmt import to_case
from ..scope import Scope
from .. import settings


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
    return to_case(basename, settings.NF_PROCESS_CASE)

def gen_varname_process(basename: str) -> str:
    return to_case(basename, settings.NF_PROCESS_CASE)

def gen_varname_channel(janis_tag: str, name_override: Optional[str]=None, dtype: Optional[DataType]=None) -> str:
    basename = name_override if name_override else janis_tag
    # basename = _handle_plurals(basename, dtype)
    name = to_case(basename, settings.NF_CHANNEL_CASE)
    name = f'ch_{name}'
    return name

def gen_varname_param(janis_tag: str, scope: Scope, name_override: Optional[str]=None, dtype: Optional[DataType]=None) -> str:
    basename = name_override if name_override else janis_tag
    # basename = _handle_plurals(basename, dtype)
    name = to_case(basename, settings.NF_PARAM_CASE)
    depth = len(scope.items)
    if depth > 1:
        scope_labels = scope.labels[1:]
        scope_labels = [to_case(x, settings.NF_PARAM_CASE) for x in scope_labels]
        name = f"{'.'.join(scope_labels)}.{name}"
    return name

# # @unused
# def _handle_plurals(basename: str, dtype: Optional[DataType]) -> str:
#     # add 's' plural for array types
#     if _should_make_plural(dtype):
#         basename = _make_plural(basename)
#     return basename

# # @unused
# def _should_make_plural(dtype: Optional[DataType]) -> bool:
#     if dtype and dtype.is_array():
#         return True
#     return False

# # @unused
# def _make_plural(basename: str) -> str:
#     if not basename.endswith('s'):
#         basename = f'{basename}s'
#     return basename


### ENTITY LABELS

