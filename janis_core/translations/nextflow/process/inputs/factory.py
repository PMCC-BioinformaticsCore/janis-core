
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
from ...scope import Scope

from .model import (
    ProcessInput, 
    PathProcessInput,
    ValProcessInput,
    TupleProcessInput
)

from . import data_sources


def create_nextflow_process_inputs(scope: Scope, tool: CommandTool | PythonTool) -> list[ProcessInput]:
    generator = ProcessInputGenerator(scope, tool)
    return generator.generate()


class ProcessInputGenerator:
    def __init__(self, scope: Scope, tool: CommandTool | PythonTool):
        self.scope = scope
        self.tool = tool

    def generate(self) -> list[ProcessInput]:
        process_inputs: list[ProcessInput] = []

        tinput_ids = data_sources.process_inputs(self.scope)
        tinputs = nfgen_utils.items_with_id(self.tool.inputs(), tinput_ids)
        tinputs = ordering.order_janis_process_inputs(tinputs)
        for inp in tinputs:
            process_inputs.append(self.create_input(inp))
        return process_inputs

    def create_input(self, inp: ToolInput | TInput) -> ProcessInput:
        dtype: DataType = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
        basetype = utils.get_base_type(dtype)
        basetype = utils.ensure_single_type(basetype)
        assert(basetype)
        
        # @secondariesarray
        # secondaries array
        if utils.is_array_secondary_type(dtype):
            return self.create_path_input_secondaries_array(inp)
        
        # secondaries
        if utils.is_secondary_type(dtype):
            return self.create_tuple_input_secondaries(inp)
        
        # filepair array
        elif dtype.name() == 'Array' and self.is_filepair_type(dtype):
            return self.create_path_input(inp)
        
        # filepair
        elif self.is_filepair_type(dtype):
            return self.create_path_input(inp)
        
        # file array
        elif dtype.is_array() and isinstance(basetype, (File, Directory)):
            return self.create_path_input(inp)
        
        # file
        elif isinstance(basetype, (File, Directory)):
            return self.create_path_input(inp)
        
        # nonfile array
        elif dtype.is_array(): 
            return self.create_val_input(inp)

        # nonfile 
        else:
            return self.create_val_input(inp)

    def is_filepair_type(self, dtype: DataType) -> bool:
        basetype = utils.get_base_type(dtype)
        basetype = utils.ensure_single_type(basetype)
        if basetype.name() in ['FastqPair', 'FastqGzPair']:
            return True
        return False

    def create_path_input_secondaries_array(self, inp: ToolInput | TInput) -> ProcessInput:
        # TODO ignoring secondaries_presents_as for now!
        name = naming.process_internal.get(self.scope, inp)
        assert(isinstance(name, str))
        new_input = PathProcessInput(name=name)
        return new_input

    def create_tuple_input_secondaries(self, inp: ToolInput | TInput) -> TupleProcessInput:
        # tuple sub-element for each file
        subnames = naming.process_internal.get(self.scope, inp)
        assert(isinstance(subnames, list))
        qualifiers = ['path'] * len(subnames)
        
        new_input = TupleProcessInput(
            name=inp.id(), 
            qualifiers=qualifiers, 
            subnames=subnames
        )
        return new_input

    def create_path_input(self, inp: ToolInput | TInput) -> PathProcessInput:
        name = naming.process_internal.get(self.scope, inp)
        assert(isinstance(name, str))
        dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype
        presents_as = None
        if isinstance(inp, ToolInput):
            presents_as = inp.presents_as
        new_input = PathProcessInput(name=name, dtype=dtype, presents_as=presents_as)
        return new_input

    def create_val_input(self, inp: ToolInput | TInput) -> ValProcessInput:
        name = naming.process_internal.get(self.scope, inp)
        assert(isinstance(name, str))
        new_input = ValProcessInput(name=name)
        return new_input
