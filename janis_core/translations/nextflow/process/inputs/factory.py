

from typing import Optional, Any

from janis_core import (
    ToolInput, 
    TInput,
    CommandTool,
    PythonTool,
    File, 
    Directory, 
    DataType
)

from janis_core import translation_utils as utils
from ... import nfgen_utils
from ... import naming
from ... import ordering

from .common import get_process_inputs
from .model import (
    ProcessInput, 
    PathProcessInput,
    ValProcessInput,
    TupleProcessInput
)


def create_nextflow_process_inputs(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> list[ProcessInput]:
    process_inputs: list[ProcessInput] = []
    tinput_ids = get_process_inputs(sources)
    tinputs = nfgen_utils.items_with_id(tool.inputs(), tinput_ids)
    tinputs = ordering.order_janis_process_inputs(tinputs)
    for i in tinputs:
        process_inputs.append(create_input(i, sources))
    return process_inputs

def create_input(inp: ToolInput | TInput, sources: dict[str, Any]) -> ProcessInput:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    basetype = utils.get_base_type(dtype)
    basetype = utils.ensure_single_type(basetype)
    assert(basetype)
    
    # @secondariesarray
    # secondaries array
    if utils.is_array_secondary_type(dtype):
        return create_path_input_secondaries_array(inp)
    
    # secondaries
    if utils.is_secondary_type(dtype):
        return create_tuple_input_secondaries(inp, sources)
    
    # filepair array
    elif dtype.name() == 'Array' and is_filepair_type(dtype):
        return create_path_input(inp)
    
    # filepair
    elif is_filepair_type(dtype):
        return create_path_input(inp)
    
    # file array
    elif dtype.is_array() and isinstance(basetype, (File, Directory)):
        return create_path_input(inp)
    
    # file
    elif isinstance(basetype, (File, Directory)):
        return create_path_input(inp)
    
    # nonfile array
    elif dtype.is_array(): 
        return create_val_input(inp)

    # nonfile 
    else:
        return create_val_input(inp)


def is_filepair_type(dtype: DataType) -> bool:
    basetype = utils.get_base_type(dtype)
    basetype = utils.ensure_single_type(basetype)
    if basetype.name() in ['FastqPair', 'FastqGzPair']:
        return True
    return False

def create_path_input_secondaries_array(inp: ToolInput | TInput) -> ProcessInput:
    # TODO ignoring secondaries_presents_as for now!
    name = naming.process_input_secondaries_array(inp)
    new_input = PathProcessInput(name=name)
    return new_input

def create_tuple_input_secondaries(inp: ToolInput | TInput, sources: dict[str, Any]) -> TupleProcessInput:
    # tuple sub-element for each file
    subnames = naming.process_input_secondaries(inp, sources)
    qualifiers = ['path'] * len(subnames)
    
    new_input = TupleProcessInput(
        name=inp.id(), 
        qualifiers=qualifiers, 
        subnames=subnames
    )
    return new_input

def create_path_input(inp: ToolInput | TInput) -> PathProcessInput:
    name = naming.process_input_name(inp)
    dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype
    presents_as = None
    if isinstance(inp, ToolInput):
        presents_as = inp.presents_as
    new_input = PathProcessInput(name=name, dtype=dtype, presents_as=presents_as)
    return new_input

def create_val_input(inp: ToolInput | TInput) -> ValProcessInput:
    name = naming.process_input_name(inp)
    new_input = ValProcessInput(name=name)
    return new_input
