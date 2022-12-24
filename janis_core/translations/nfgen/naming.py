

from typing import Optional, Any

from janis_core.workflow.workflow import StepNode
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

def get_varname_workflow(basename: str) -> str:
    return to_case(basename, settings.NF_PROCESS_CASE)

def get_varname_process(basename: str) -> str:
    return to_case(basename, settings.NF_PROCESS_CASE)

def get_varname_channel(basename: str) -> str:
    name = to_case(basename, settings.NF_CHANNEL_CASE)
    name = f'ch_{name}'
    return name

def get_varname_param(basename: str, scope: list[str]) -> str:
    name = to_case(basename, settings.NF_PARAM_CASE)
    if len(scope) > 1:
        scope = scope[1:]
        scope = [to_case(x, settings.NF_PARAM_CASE) for x in scope]
        name = f"{'.'.join(scope)}.{name}"
    return name



### PROCESS SPECIFIC

def get_varname_secondaries(dtype: DataType) -> list[str]:
    """returns name of each file for File types with secondaries"""
    basetype: File = nfgen_utils.get_base_type(dtype)
    exts = nfgen_utils.get_extensions(basetype)
    exts = _remove_symbols(exts)
    
    # Array(Secondary())
    if isinstance(dtype, Array):
        exts = [f'{ext}s' for ext in exts]
        return exts
    
    # Secondary
    else:
        return exts
    
def _remove_symbols(exts: list[str]) -> list[str]:
    return [x.rsplit('.')[-1] for x in exts]

def get_varname_toolinput(inp: ToolInput | TInput, process_inputs: set[str], param_inputs: set[str], sources: dict[str, Any]) -> Optional[str]:
    if inp.id() in process_inputs:
        return _process_input_varname(inp, process_inputs, param_inputs)
    elif inp.id() in param_inputs:
        return _param_input_varname(inp, sources)
    else:
        # internal inputs in the future
        raise NotImplementedError

def _process_input_varname(inp: ToolInput | TInput, process_inputs: set[str], param_inputs: set[str]) -> Optional[str]:
    # data fed via process input
    name: Optional[str] = None
    if inp.id() in process_inputs:
        # secondary files (name mapped to ext of primary file)
        dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype
        basetype = nfgen_utils.get_base_type(dtype)
        if isinstance(basetype, File) and basetype.has_secondary_files():
            names = get_varname_secondaries(basetype)
            name = names[0]
        # everything else
        else:
            name = inp.id()
    return name
    
def _param_input_varname(inp: ToolInput | TInput, sources: dict[str, Any]) -> Optional[str]: 
    # data fed via global param
    src = sources[inp.id()]
    sel = src.source_map[0].source
    param = params.get(sel.input_node.uuid)
    return f'params.{param.name}'



### ENTITY LABELS

def get_construct_name(tool: CommandTool | PythonTool | Workflow, scope: list[str]) -> str:
    construct_type = ''
    if isinstance(tool, CommandTool) or isinstance(tool, PythonTool):
        construct_type = 'process'
    elif isinstance(tool, Workflow) and len(scope) == 1:  # scope = ['main']  (the main workflow)
        construct_type = 'main_workflow'
    elif isinstance(tool, Workflow) and len(scope) > 1: # scope = ['main', 'sub', ...] 
        construct_type = 'sub_workflow'
    else:
        raise NotImplementedError
    return construct_type
