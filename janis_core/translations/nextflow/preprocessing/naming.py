


from typing import Optional, Any

from janis_core import (
    DataType,
    TInput,
)

from janis_core import translation_utils as utils
from janis_core import settings

from .. import params
from ..scope import Scope
from ..casefmt import to_case


from copy import deepcopy
from typing import Any

import janis_core.translation_utils as utils
from janis_core.workflow.workflow import Workflow
from janis_core import CommandTool, PythonTool, TInput
from janis_core import settings

from ..scope import Scope
from .. import naming
from .. import params
from ..process.inputs import data_sources


def register_process_internal_names(wf: Workflow) -> None:
    """
    register param(s) for each workflow input. 
    channel(s) may also be registered if necessary.
    """
    scope = Scope()
    do_register_process_internal_names(scope, wf)

def do_register_process_internal_names(scope: Scope, wf: Workflow):
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        
        if isinstance(step.tool, CommandTool) or isinstance(step.tool, PythonTool):
            generator = ProcessInternalNameGenerator(current_scope, step.tool, step.sources)
            generator.generate()
            for tinput_id, variable_name in generator.names_map.items():
                naming.process_internal.update(current_scope, tinput_id, variable_name)
        
        elif isinstance(step.tool, Workflow):
            do_register_process_internal_names(current_scope, step.tool)
        else:
            raise NotImplementedError


class ProcessInternalNameGenerator:
    """
    for each CommandTool / PythonTool, determines the internal variable_name for each tinput.id()
    the variable_name is how the particular tinput_id will be referenced inside the process. 
    """
    def __init__(self, scope: Scope, tool: CommandTool | PythonTool, sources: dict[str, Any]) -> None:
        self.scope = scope
        self.tool = tool
        self.sources = sources

        self.process_inputs = data_sources.process_inputs(self.scope)
        self.param_inputs = data_sources.param_inputs(self.scope)
        self.internal_inputs = data_sources.internal_inputs(self.scope)
        
        # want to calculate this
        self.names_map: dict[str, Optional[str | list[str]]] = {}

    # helper properties
    @property
    def tinputs(self) -> list[TInput]:
        if isinstance(self.tool, CommandTool):
            return self.tool.tool_inputs()
        elif isinstance(self.tool, PythonTool):  # type: ignore
            return self.tool.inputs()
        else:
            raise NotImplementedError
    
    # helper methods
    def duplicate_datatype_exists(self, inp: TInput) -> bool:
        """
        check if another TInput has the same dtype as this TInput.
        only checking for secondary and array secondary datatype duplicates.
        """
        rtype = inp.intype  # type: ignore
        rbasetype = utils.get_base_type(rtype)  # type: ignore
        
        for tinput in self.tinputs:
            # dont check tinput against itself
            if tinput.id() == inp.id():
                continue
            
            # check if types match
            qtype = tinput.intype  # type: ignore
            qbasetype = utils.get_base_type(qtype)  # type: ignore

            if utils.is_array_secondary_type(rtype) and utils.is_array_secondary_type(qtype):  # type: ignore
                if type(rbasetype) == type(qbasetype):
                    return True
            
            elif utils.is_secondary_type(rtype) and utils.is_secondary_type(qtype):  # type: ignore
                if type(rtype) == type(qtype):  # type: ignore
                    return True
                
        return False
        
    def generate(self) -> None:
        for tinput in self.tinputs:
            name = self.gen_varname(tinput)
            self.names_map[tinput.id()] = name

    # this should be kept in ProcessVariableManager or something
    def gen_varname(self, inp: TInput) -> Optional[str | list[str]]:
        if inp.id() in self.process_inputs:
            return self.gen_varname_process_input(inp)
        elif inp.id() in self.param_inputs:
            return self.gen_varname_param_input(inp)
        elif inp.id() in self.internal_inputs:
            return None
        else:
            raise RuntimeError
        
    def gen_varname_process_input(self, inp: TInput) -> str | list[str]:
        # data fed via process input
        dtype: DataType = inp.intype  # type: ignore
        
        if utils.is_array_secondary_type(dtype):
            name = self.gen_varname_process_input_secondaries_array(inp)
        elif utils.is_secondary_type(dtype):
            name = self.gen_varname_process_input_secondaries(inp)  # first extension
        else:
            name = self.gen_varname_process_input_generic(inp)  
        
        return name

    def gen_varname_process_input_generic(self, inp: TInput) -> str:
        return to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)

    def gen_varname_process_input_secondaries(self, inp: TInput) -> list[str]:
        """returns name of each file for File types with secondaries"""
        # src = sources[inp.id()]
        # srctype: DataType = get_src_type(src)

        # # datatype mismatch! get type info from srctype
        # if srctype.name() != desttype.name():
        #     basetype: File = utils.get_base_type(srctype)  # type: ignore 
        
        # # datatype match. get type info from dest
        # else:
        #     basetype: File = utils.get_base_type(desttype)  # type: ignore 
        
        
        dtype: DataType = inp.intype  # type: ignore
        # basetype: File = utils.get_base_type(dtype) # type: ignore 
        exts = utils.get_extensions(dtype, remove_prefix_symbols=True)
        exts = [x.replace('.', '_') for x in exts]
        
        # ['bam', 'bai'] -> ['normal_bam', 'normal_bai'] when other tinput 
        # with same secondary datatype exists
        if self.duplicate_datatype_exists(inp):
            identifier = to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
            exts = [f'{identifier}_{x}' for x in exts]

        return exts

    def gen_varname_process_input_secondaries_array(self, inp: TInput) -> str:
        if self.duplicate_datatype_exists(inp):
            identifier = to_case(inp.id(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
            return f'{identifier}_flat'
        else:
            dtype: DataType = inp.intype  # type: ignore
            basetype: DataType = utils.get_base_type(dtype)  # type: ignore
            basetype_name = to_case(basetype.name(), case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
            return f'{basetype_name}_flat'
    
    def gen_varname_param_input(self, inp: TInput) -> str: 
        # data fed via global param
        src = self.sources[inp.id()]
        sel = src.source_map[0].source
        param = params.get(sel.input_node.uuid)
        pname = to_case(param.name, case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
        return f'params.{pname}'






# def process_input_secondaries_array_primary_files(inp: ToolInput | TInput) -> str:
#     """
#     example: ToolInput = Array(BamBai)
#     process input name: indexed_bam_array_flat
#     primary files name: bams
#     """
#     dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype  # type: ignore
#     basetype: File = utils.get_base_type(dtype)  # type: ignore 
#     exts = utils.get_extensions(basetype, remove_symbols=True)
#     primary_ext = exts[0]
#     primary_name = f'{primary_ext}s' # primary extension -> make plural
#     return primary_name

# def get_src_type(src: Any) -> DataType:
#     # the srctype corresponds to either a workflow input, or step output.
#     # scattering doesn't matter. 
#     dtype = trace.trace_source_datatype(src)
#     if isinstance(dtype, Stdout):
#         return dtype.subtype
#     else:
#         return dtype

