

from janis_core import (
    ToolInput, 
    TInput,
    CommandTool,
    CommandToolBuilder,
    PythonTool,
    File, 
    DataType
)
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType
from janis_core import settings

from ... import naming
from ... import nfgen_utils
from ... import task_inputs

from ...model.process.inputs import (
    NFProcessInput, 
    NFPythonToolProcessInput,
    NFScriptProcessInput,
    NFPathProcessInput,
    NFValProcessInput,
    NFTupleProcessInput,
)


def gen_nf_process_inputs(tool: CommandTool | PythonTool) -> list[NFProcessInput]:
    generator = ProcessInputGenerator(tool)
    return generator.generate()


class ProcessInputGenerator:
    def __init__(self, tool: CommandTool | PythonTool):
        self.tool = tool
        self.process_inputs: list[NFProcessInput] = []

    @property
    def dtype(self) -> DataType:
        return self.tinput.input_type if isinstance(self.tinput, ToolInput) else self.tinput.intype # type: ignore
    
    @property
    def basetype(self) -> DataType:
        basetype = utils.get_base_type(self.dtype)
        basetype = utils.ensure_single_type(basetype)
        assert(basetype)
        return basetype
    
    def generate(self) -> list[NFProcessInput]:
        # CodeTools get additional 'code_file' input
        if isinstance(self.tool, PythonTool):
            self.generate_code_file_input()
        
        # files_to_create inputs
        if isinstance(self.tool, CommandTool):
            self.generate_files_to_create_inputs()

        # regular inputs
        self.generate_regular_inputs()
        
        return self.process_inputs

    def generate_code_file_input(self) -> None:
        # pythontool gets extra code_file input before normal inputs
        new_input = NFPythonToolProcessInput(
            name=settings.translate.nextflow.PYTHON_CODE_FILE, 
            tinput_id=settings.translate.nextflow.PYTHON_CODE_FILE, 
            dtype=File())
        self.process_inputs.append(new_input)

    def generate_files_to_create_inputs(self) -> None:
        if isinstance(self.tool, CommandToolBuilder) and not settings.ingest.SOURCE == 'galaxy':
            if self.tool._files_to_create:
                assert isinstance(self.tool._files_to_create, dict)
                for filename, filecontents in self.tool._files_to_create.items():
                    
                    # ignoring shell script parsed from WDL
                    if self.tool.is_shell_script and filename == 'script.sh':
                        continue 

                    # generate a name for this input
                    if len(self.tool._files_to_create) == 1:
                        name = 'script'
                    else:
                        name = naming.process.files_to_create_script(filename)
                    
                    # create the nf process input 
                    new_input = NFScriptProcessInput(
                        name=name,
                        tinput_id=name,
                        dtype=File(), 
                        presents_as=filename
                    )
                    self.process_inputs.append(new_input)
        
    def generate_regular_inputs(self) -> None:
        tinput_ids = task_inputs.task_inputs(self.tool.id())
        tinputs = nfgen_utils.items_with_id(self.tool.inputs(), tinput_ids)
        
        for inp in tinputs:
            self.tinput = inp
            new_input = self.create_input()
            self.process_inputs.append(new_input)

    def create_input(self) -> NFProcessInput:
        """create the correct NFProcessInput for the given input"""
        dtt = utils.get_dtt(self.dtype)

        if dtt == DTypeType.SECONDARY_ARRAY:
            # return self.create_path_input_secondaries_array(self.tinput)
            return self.create_path_input(self.tinput)
        
        # secondaries
        elif dtt == DTypeType.SECONDARY:
            # return self.create_tuple_input_secondaries(self.tinput)
            return self.create_path_input(self.tinput)
        
        # filepair array
        elif dtt == DTypeType.FILE_PAIR_ARRAY:
            return self.create_path_input_file_pair_array(self.tinput)
        
        # filepair
        elif dtt == DTypeType.FILE_PAIR:
            return self.create_tuple_input_file_pair(self.tinput)
        
        # file array
        elif dtt == DTypeType.FILE_ARRAY:
            return self.create_path_input(self.tinput)
        
        # file
        elif dtt == DTypeType.FILE:
            return self.create_path_input(self.tinput)
        
        # filename
        elif dtt == DTypeType.FILENAME:
            return self.create_val_input(self.tinput)
        
        # nonfile array
        elif dtt == DTypeType.GENERIC_ARRAY:
            return self.create_val_input(self.tinput)

        # nonfile 
        elif dtt == DTypeType.GENERIC:
            return self.create_val_input(self.tinput)
        
        elif dtt == DTypeType.FLAG_ARRAY:
            return self.create_val_input(self.tinput)
        
        elif dtt == DTypeType.FLAG:
            return self.create_val_input(self.tinput)
        
        else:
            raise NotImplementedError(f"Unknown input type: {dtt}")

    # def create_path_input_secondaries_array(self, inp: ToolInput | TInput) -> NFPathProcessInput:
    #     # TODO ignoring secondaries_presents_as for now
    #     ti = task_inputs.get(self.tool.id(), inp)
    #     name = ti.value
    #     assert(isinstance(name, str))
    #     new_input = NFPathProcessInput(name=name, tinput_id=inp.id(), dtype=self.dtype)
    #     return new_input
    
    # def create_tuple_input_secondaries(self, inp: ToolInput | TInput) -> NFTupleProcessInput:
    #     # tuple sub-element for each file
    #     ti = task_inputs.get(self.tool.id(), inp)
    #     subnames = ti.value
    #     assert(isinstance(subnames, list))
        
    #     new_input = NFTupleProcessInput(
    #         name=inp.id(), 
    #         tinput_id=inp.id(),
    #         dtype=self.dtype,
    #         subnames=subnames
    #     )
    #     return new_input
    
    def create_path_input_file_pair_array(self, inp: ToolInput | TInput) -> NFPathProcessInput:
        ti = task_inputs.get(self.tool.id(), inp)
        name = ti.value
        assert(isinstance(name, str))
        new_input = NFPathProcessInput(name=name, tinput_id=inp.id(), dtype=self.dtype)
        return new_input
    
    def create_tuple_input_file_pair(self, inp: ToolInput | TInput) -> NFTupleProcessInput:
        # tuple sub-element for each file
        ti = task_inputs.get(self.tool.id(), inp)
        assert(ti.value)
        assert(len(ti.value) == 2)
        
        new_input = NFTupleProcessInput(
            name=inp.id(), 
            tinput_id=inp.id(),
            dtype=self.dtype,
            subnames=[ti.value[0], ti.value[1]]
        )
        return new_input

    def create_path_input(self, inp: ToolInput | TInput) -> NFPathProcessInput:
        ti = task_inputs.get(self.tool.id(), inp)
        name = ti.value
        assert(isinstance(name, str))
        presents_as = None
        if isinstance(inp, ToolInput):
            presents_as = inp.presents_as
        new_input = NFPathProcessInput(name=name, tinput_id=inp.id(), dtype=self.dtype, presents_as=presents_as)
        return new_input

    def create_val_input(self, inp: ToolInput | TInput) -> NFValProcessInput:
        if isinstance(self.dtype, File) and isinstance(inp, ToolInput):
            if inp.presents_as:
                raise NotImplementedError
                presents_as = inp.presents_as
        ti = task_inputs.get(self.tool.id(), inp)
        name = ti.value
        assert(isinstance(name, str))
        new_input = NFValProcessInput(name=name, tinput_id=inp.id(), dtype=self.dtype)
        return new_input
