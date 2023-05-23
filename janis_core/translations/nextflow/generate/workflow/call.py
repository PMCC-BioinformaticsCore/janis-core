


from typing import Any, Optional
from textwrap import indent

from janis_core import CommandTool, PythonTool, Workflow, ScatterDescription, TInput
from janis_core.workflow.workflow import StepNode
from janis_core.types import DataType, Stdout
from janis_core import settings

from ....common import trace
from ... import params

from ...model.process import NFProcess
from ...model.process import NFProcessInput
from ...model.process import NFPythonToolProcessInput
from ...model.process import NFScriptProcessInput
from ...model.workflow import NFWorkflow
from ...model.workflow import NFWorkflowTake

from ...unwrap import unwrap_expression
from ...variables import VariableManager

from .datatype_mismatch import is_datatype_mismatch
from .datatype_mismatch import gen_datatype_mismatch_plumbing

from .edge_cases import satisfies_edge_case
from .edge_cases import handle_edge_case
from ... import nulls

NF_INDENT = settings.translate.nextflow.NF_INDENT



def gen_task_call(alias: str, task: NFProcess | NFWorkflow, vmanager: VariableManager, step: StepNode) -> list[str]:
    generator = TaskCallGenerator(alias, task, vmanager, step)
    return generator.generate()


class TaskCallGenerator:
    def __init__(self, alias: str, task: NFProcess | NFWorkflow, vmanager: VariableManager, step: StepNode) -> None:
        self.alias = alias
        self.task = task
        self.step = step
        self.vmanager = vmanager
        
        self.tool: CommandTool | PythonTool | Workflow  = step.tool     
        self.sources: dict[str, Any]                    = step.sources
        
        # want to calculate these
        self.args: list[str] = []
        self.call: list[str] = []

    @property
    def ordered_task_inputs(self) -> list[NFProcessInput | NFWorkflowTake]:
        if isinstance(self.task, NFProcess):
            return [x for x in self.task.ordered_inputs]
        elif isinstance(self.task, NFWorkflow): # type: ignore
            return [x for x in self.task.ordered_take]
        else:
            raise RuntimeError

    def generate(self) -> list[str]:
        self.gen_call_arguments()
        self.format_task_call()
        return self.call
    
    def gen_call_arguments(self) -> None:
        call_args: list[str] = []

        for task_input in self.ordered_task_inputs:
            arg = self.get_call_arg(task_input)
            call_args.append(arg)
        
        self.args = call_args
    
    def get_call_arg(self, task_input: NFProcessInput | NFWorkflowTake) -> str:
        generator = TaskCallArgumentGenerator(
            task_input=task_input,
            vmanager=self.vmanager,
            # calling_scope=self.calling_scope,
            step=self.step
        )
        return generator.generate()

    # formatting task call text
    def format_task_call(self, ind: int=0) -> None:    
        if len(self.args) == 0:
            call_lines = self.call_fmt0()
        else:
            call_lines = self.call_fmt2()

        # indenting for workflow scope
        call_lines = [indent(ln, ind * NF_INDENT) for ln in call_lines]
        self.call = call_lines

    def call_fmt0(self) -> list[str]:
        return [f'{self.alias}()']

    def call_fmt2(self) -> list[str]:
        arg_lines = self.args

        # adding commas for each argument
        arg_lines = self.add_commas(arg_lines)

        # adding comments for each argument
        if settings.translate.RENDER_COMMENTS:
            arg_lines = self.add_comments(arg_lines)

        # indenting each argument
        arg_lines = self.indent_arg_lines(arg_lines)
        
        # formatting call lines
        call_lines: list[str] = []
        call_lines.append(f'{self.alias}(')
        call_lines += arg_lines
        call_lines.append(')')
        return call_lines
    
    def add_commas(self, arg_lines: list[str]) -> list[str]:
        out: list[str] = []
        for i, inp in enumerate(arg_lines):
            comma = ',' if i < len(arg_lines) - 1 else ''
            out.append(f'{inp}{comma}')
        return out
    
    def add_comments(self, arg_lines: list[str]) -> list[str]:
        assert(len(arg_lines) == len(self.ordered_task_inputs))
        max_width = max([len(x) for x in arg_lines])
        out: list[str] = []
        for arg, task_input in zip(arg_lines, self.ordered_task_inputs):
            arg = f'{arg:<{max_width + 2}}// {task_input.tinput_id}'
            out.append(arg)
        return out
    
    def indent_arg_lines(self, arg_lines: list[str]) -> list[str]:
        return [f'{NF_INDENT}{arg}' for arg in arg_lines]

 

class TaskCallArgumentGenerator:
    def __init__(self, task_input: NFProcessInput | NFWorkflowTake, vmanager: VariableManager, step: StepNode) -> None:
        self.task_input = task_input
        self.vmanager = vmanager

        self.tool: CommandTool | PythonTool | Workflow  = step.tool     
        self.scatter: Optional[ScatterDescription]      = step.scatter
        self.src: Optional[Any]                         = None 

        # update if has source
        if self.tinput_id in step.sources:
            self.src = step.sources[self.tinput_id]

    @property
    def tinput_id(self) -> str:
        return self.task_input.tinput_id

    @property
    def tinput(self) -> TInput:
        return [x for x in self.tool.tool_inputs() if x.id() == self.tinput_id][0]

    @property
    def srctype(self) -> Optional[DataType]:
        """the datatype of the data source"""
        if self.src:
            dtype = trace.trace_source_datatype(self.src)
            assert(dtype)
            if isinstance(dtype, Stdout):
                return dtype.subtype
            else:
                return dtype
        return None
        
    @property
    def desttype(self) -> DataType:
        """the datatype of the relevant ToolInput"""
        tinputs = self.tool.tool_inputs()
        tinp = [x for x in tinputs if x.id() == self.tinput_id][0]
        return tinp.intype  # type: ignore
    
    # @property
    # def src_scatter(self) -> bool:
    #     return trace.trace_source_scatter(self.src)

    @property
    def dest_scatter(self) -> bool:
        if self.scatter and self.tinput_id in self.scatter.fields:
            return True
        return False

    def generate(self) -> str:
        """calculate the arg which will feed this task input"""
        if isinstance(self.task_input, NFPythonToolProcessInput):
            return self.generate_code_file()
        elif isinstance(self.task_input, NFScriptProcessInput):
            return self.generate_script()
        else:
            return self.generate_regular()

    def generate_script(self) -> str:
        param = params.get(
            tinput_id=self.task_input.tinput_id,
            task_id=self.tool.id()
        )
        return f'params.{param.name}'
    
    def generate_code_file(self) -> str:
        param = params.get(
            tinput_id=self.task_input.tinput_id, 
            task_id=self.tool.id()
        )
        return f'params.{param.name}'

    def generate_regular(self) -> str:
        arg = unwrap_expression(
            val=self.src,
            context='workflow',
            variable_manager=self.vmanager,
            quote_strings=True
        )
            
        if isinstance(arg, list):
            raise NotImplementedError
            call_args += arg

        if arg is not None:
            # handle misc edge case (takes priority over datatype mismatches)
            if satisfies_edge_case(self.src):
                # TODO optionality checking should be inside handle_edge_case()
                # if not self.srctype.optional:
                suffix = handle_edge_case(self.src)
                arg = f'{arg}{suffix}'

            # handle datatype relationship
            elif self.srctype is not None and is_datatype_mismatch(self.srctype, self.desttype, self.dest_scatter):
                # TODO optionality checking should be inside gen_datatype_mismatch_plumbing()
                # if not self.srctype.optional:
                suffix = gen_datatype_mismatch_plumbing(self.srctype, self.desttype, self.dest_scatter)
                arg = f'{arg}{suffix}'
        
        if arg is None:
            arg = nulls.get_null_value(self.desttype, as_param=True, should_add_file_cast=True)
        
        return arg
    
    

    




