

from typing import Optional, Any

from janis_core.types import (
    File,
    Array
)
from janis_core import (
    DataType,
    ToolInput,
    TInput,
    Workflow,
    CommandTool,
    PythonTool,
)

from .casefmt import to_case
from .scope import Scope
from . import nfgen_utils
from . import params
from . import settings


"""

Fixing these items 

- array naming 
    - markDuplicates.nf
        inputs:
        path bam -> path bams
    
- base_recalibrator.nf
    def knownSites = gz.collect{ "--known-sites " + it }.join(' ')
    -> 
    def knownSites_gzs = knownSites_gzs.collect{ "--known-sites " + it }.join(' ')

- param naming
    - params.step.pname format    (params.align_and_sort.sortsam_tmp_dir, params.unicycler.kmers etc)
    - // process fastqc   // subworkflow align_and_sort etc headings 
    - fastqs = []  -> 
      fastqs = [
          // list files here
      ]

"""


### GENERAL

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

# @unused
def _handle_plurals(basename: str, dtype: Optional[DataType]) -> str:
    # add 's' plural for array types
    if _should_make_plural(dtype):
        basename = _make_plural(basename)
    return basename

# @unused
def _should_make_plural(dtype: Optional[DataType]) -> bool:
    if dtype and dtype.is_array():
        return True
    return False

# @unused
def _make_plural(basename: str) -> str:
    if not basename.endswith('s'):
        basename = f'{basename}s'
    return basename


### PROCESS SPECIFIC

def gen_varname_toolinput(inp: ToolInput | TInput, process_inputs: set[str], param_inputs: set[str], sources: dict[str, Any]) -> Optional[str]:
    if inp.id() in process_inputs:
        return _process_input_varname(inp, process_inputs, param_inputs)
    elif inp.id() in param_inputs:
        return _param_input_varname(inp, sources)
    else:
        return None
        # internal inputs in the future
        raise NotImplementedError

def _param_input_varname(inp: ToolInput | TInput, sources: dict[str, Any]) -> Optional[str]: 
    # data fed via global param
    src = sources[inp.id()]
    sel = src.source_map[0].source
    param = params.get(sel.input_node.uuid)
    return f'params.{param.name}'

def _process_input_varname(inp: ToolInput | TInput, process_inputs: set[str], param_inputs: set[str]) -> Optional[str]:
    # data fed via process input
    name: Optional[str] = None
    if inp.id() in process_inputs:
        name = inp.id()
        
        # # secondary files (name mapped to ext of primary file)
        # dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype
        # basetype = nfgen_utils.get_base_type(dtype)
        # if isinstance(basetype, File) and basetype.has_secondary_files():
        #     names = _gen_varname_toolinput_secondaries(dtype)
        #     name = names[0]
        # # everything else
        # else:
        #     # name = _handle_plurals(name, dtype)

    return name

def gen_varname_toolinput_secondaries(dtype: DataType) -> list[str]:
    """returns name of each file for File types with secondaries"""
    basetype: File = nfgen_utils.get_base_type(dtype)
    exts = nfgen_utils.get_extensions(basetype)
    exts = _remove_symbols(exts)
    # # @unused    
    # # Array(Secondary())
    # # plurals handled here as special case (not _handle_plurals())
    # if isinstance(dtype, Array):
    #     exts = [f'{ext}s' for ext in exts]
    #     return exts
    return exts

    
    
# @unused
def _remove_symbols(exts: list[str]) -> list[str]:
    return [x.rsplit('.')[-1] for x in exts]
    


### ENTITY LABELS

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
