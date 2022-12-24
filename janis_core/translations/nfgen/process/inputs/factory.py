

from typing import Optional, Any

from janis_core import (
    ToolInput, 
    TInput,
    CommandTool,
    PythonTool,
    File, 
    Directory, 
    Array, 
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
        process_inputs += create_inputs(i)
    return process_inputs

def create_inputs(inp: ToolInput | TInput) -> list[ProcessInput]:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    datatype: DataType = dtype
    if isinstance(datatype, Array):
        return create_inputs_array(inp)
    else:
        return create_inputs_single(inp)

def create_inputs_array(inp: ToolInput | TInput) -> list[ProcessInput]:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    basetype: Optional[DataType] = nfgen_utils.get_base_type(dtype)
    assert(basetype)

    # secondaries array
    if isinstance(basetype, File) and basetype.has_secondary_files():
        # a path input per file type
        return create_path_input_secondaries(inp)

    # file array
    if isinstance(basetype, (File, Directory)):
        return [create_path_input(inp)]

    # nonfile array
    return [create_val_input(inp)]


def create_inputs_single(inp: ToolInput | TInput) -> list[ProcessInput]:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    basetype: Optional[DataType] = nfgen_utils.get_base_type(dtype)
    assert(basetype)
        
    # file secondaries
    if isinstance(basetype, File) and basetype.has_secondary_files():
        inputs = [create_tuple_input_secondaries(inp)]
        return inputs # type: ignore
    
    # file
    if isinstance(basetype, (File, Directory)):
        return [create_path_input(inp)]
    
    # nonfile
    return [create_val_input(inp)]


def create_path_input(inp: ToolInput | TInput) -> PathProcessInput:
    new_input = PathProcessInput(name=inp.id())
    new_input.presents_as = None
    if isinstance(inp, ToolInput):
        new_input.presents_as = inp.presents_as
    return new_input

def create_val_input(inp: ToolInput | TInput) -> ValProcessInput:
    new_input = ValProcessInput(name=inp.id())
    return new_input

def create_path_input_secondaries(inp: ToolInput | TInput) -> list[ProcessInput]:
    # TODO ignoring secondaries_presents_as for now!
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    inputs: list[ProcessInput] = []
    names = naming.get_varname_secondaries(dtype)
    for name in names:
        new_input = PathProcessInput(name=name)
        inputs.append(new_input)
    return inputs

def create_tuple_input_secondaries(inp: ToolInput | TInput) -> TupleProcessInput:
    dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    assert(isinstance(dtype, File))
    qualifiers: list[str] = []
    subnames: list[str] = []

    # tuple sub-element for each file
    names = naming.get_varname_secondaries(dtype)
    for name in names:
        qualifiers.append('path')
        subnames.append(name)
    
    new_input = TupleProcessInput(
        name=inp.id(), 
        qualifiers=qualifiers, 
        subnames=subnames
    )
    return new_input