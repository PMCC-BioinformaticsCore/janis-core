

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

from ... import nfgen_utils
from ... import naming
from ... import ordering

from .janis import get_process_inputs
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
        process_inputs.append(create_input(i))
    return process_inputs

def create_input(inp: ToolInput | TInput) -> ProcessInput:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    basetype: Optional[DataType] = nfgen_utils.get_base_type(dtype)
    assert(basetype)
    
    # @secondariesarray
    # secondaries array
    if nfgen_utils.is_array_secondary_type(dtype):
        return create_path_input_secondaries_array(inp)
    
    # secondaries
    if nfgen_utils.is_secondary_type(dtype):
        return create_tuple_input_secondaries(inp)
    
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

def create_path_input_secondaries_array(inp: ToolInput | TInput) -> ProcessInput:
    # TODO ignoring secondaries_presents_as for now!
    name = naming.process_input_secondaries_array(inp)
    new_input = PathProcessInput(name=name)
    return new_input

def create_tuple_input_secondaries(inp: ToolInput | TInput) -> TupleProcessInput:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    assert(isinstance(dtype, File))

    # tuple sub-element for each file
    subnames = naming.process_input_secondaries(dtype)
    qualifiers = ['path'] * len(subnames)
    
    new_input = TupleProcessInput(
        name=inp.id(), 
        qualifiers=qualifiers, 
        subnames=subnames
    )
    return new_input

def create_path_input(inp: ToolInput | TInput) -> PathProcessInput:
    name = naming.process_input_generic(inp)
    new_input = PathProcessInput(name=name)
    new_input.presents_as = None
    if isinstance(inp, ToolInput):
        new_input.presents_as = inp.presents_as
    return new_input

def create_val_input(inp: ToolInput | TInput) -> ValProcessInput:
    name = naming.process_input_generic(inp)
    new_input = ValProcessInput(name=name)
    return new_input
