


from typing import Any, Optional
from textwrap import indent
from copy import deepcopy

from janis_core import CommandTool, PythonTool, Workflow, ScatterDescription
from janis_core.workflow.workflow import StepNode
from janis_core.types import DataType, Stdout
from janis_core import settings

from .. import data_sources
from .. import nfgen_utils
from .. import ordering

from ..unwrap import unwrap_expression
from ..scope import Scope

from .. import trace
from .datatype_mismatch import is_datatype_mismatch
from .datatype_mismatch import gen_datatype_mismatch_plumbing

from .edge_cases import satisfies_edge_case
from .edge_cases import handle_edge_case



# from .scatter import is_scatter_relationship
# from .scatter import handle_scatter_relationship

NF_INDENT = settings.translate.nextflow.NF_INDENT

 
def gen_task_call(step: StepNode, scope: Scope, entity_name: str) -> str:
    generator = TaskCallGenerator(step, scope, entity_name)
    return generator.generate()


# TODO HERE

class TaskCallGenerator:
    def __init__(self, step: StepNode, scope: Scope, entity_name: str) -> None:
        self.step = step
        self.scope = scope
        self.entity_name = entity_name
        
        self.tool: CommandTool | PythonTool | Workflow  = self.step.tool     
        self.sources: dict[str, Any]                    = self.step.sources
        
        # want to calculate these
        self.args: list[str]                            = []
        self.call: str = ''

    @property
    def workflow_scope(self) -> Scope:
        return self.scope
    
    @property
    def task_scope(self) -> Scope:
        task_scope = deepcopy(self.scope)
        task_scope.update(self.step)
        return task_scope

    @property
    def ordered_process_input_ids(self) -> list[str]:
        task_ids = data_sources.process_inputs(self.task_scope)
        task_inputs = nfgen_utils.items_with_id(self.tool.tool_inputs(), task_ids)
        if isinstance(self.tool, Workflow):
            task_inputs = ordering.order_workflow_inputs(task_inputs)
        else:
            task_inputs = ordering.order_janis_process_inputs(task_inputs)
        return [x.id() for x in task_inputs]

    def generate(self) -> str:
        self.args = self.get_call_arguments()
        self.call = self.format_process_call()
        return self.call
    
    def get_call_arguments(self) -> list[str]:
        call_args: list[str] = []

        for tinput_id in self.ordered_process_input_ids:
            arg = self.get_call_arg(tinput_id)
            call_args.append(arg)
        
        # add extra arg in case of python tool - the code file.
        # a param with the same name will have already been created. 
        if isinstance(self.tool, PythonTool):
            scope_joined = self.task_scope.to_string(ignore_base_item=True)
            call_args = [f'params.{scope_joined}.code_file'] + call_args
        
        return call_args
    
    def get_call_arg(self, tinput_id: str) -> str:
        generator = TaskCallArgumentGenerator(
            tinput_id=tinput_id,
            scope=self.workflow_scope,
            step=self.step
        )
        return generator.generate()

    # formatting process call text
    def format_process_call(self, ind: int=0) -> str:
        if len(self.args) == 0:
            call_str = self.call_fmt0()
        else:
            call_str = self.call_fmt2()
        # elif len(inputs) == 1:
        #     call_str = call_fmt1(name, inputs[0])
        # elif len(inputs) > 1:
            # call_str = call_fmt2(name, inputs)
        return indent(call_str, ind * NF_INDENT)

    def call_fmt0(self) -> str:
        return f'{self.entity_name}()\n'

    def call_fmt2(self) -> str:
        call_str = f'{self.entity_name}(\n'
        for i, inp in enumerate(self.args):
            comma = ',' if i < len(self.args) - 1 else ''
            call_str += f'{NF_INDENT}{inp}{comma}\n'
        call_str += ')\n'
        return call_str

    # def _call_fmt1(self) -> str:
    #     return f'{self.entity_name}( {self.args[0]} )\n'

    # # helpers:
    # # identifying which tool inputs we need to provide an arg for 
    # def get_input_ids(self) -> list[str]:
    #     # input ids which we need args for (in correct order)
    #     if isinstance(self.tool, Workflow):
    #         return self.get_input_ids_workflow()
    #     else:
    #         return self.get_input_ids_tool()

    # def get_input_ids_workflow(self) -> list[str]:
    #     # sub Workflow - order via workflow inputs
    #     # (ignore workflow inputs which don't appear in the original janis step call)
    #     subwf: Workflow = self.tool  # type: ignore
    #     subwf_ids = set(self.sources.keys())
    #     subwf_inputs = nfgen_utils.items_with_id(list(subwf.input_nodes.values()), subwf_ids)
    #     subwf_inputs = ordering.order_workflow_inputs(subwf_inputs)
    #     return [x.id() for x in subwf_inputs]

    # def get_input_ids_tool(self) -> list[str]:
    #     # CommandTool / PythonTool - order via process inputs 
    #     tool: CommandTool | PythonTool = self.tool  # type: ignore
    #     process_ids = data_sources.process_inputs(self.scope)
    #     process_inputs = nfgen_utils.items_with_id(tool.inputs(), process_ids)
    #     process_inputs = ordering.order_janis_process_inputs(process_inputs)
    #     return [x.id() for x in process_inputs]


class TaskCallArgumentGenerator:
    def __init__(self, tinput_id: str, scope: Scope, step: StepNode) -> None:
        self.tinput_id = tinput_id
        self.scope = scope

        self.tool: CommandTool | PythonTool | Workflow  = step.tool     
        self.scatter: Optional[ScatterDescription]      = step.scatter
        self.src: Any                                   = step.sources[self.tinput_id]

    @property
    def srctype(self) -> DataType:
        """the datatype of the data source"""
        dtype = trace.trace_source_datatype(self.src)
        assert(dtype)
        if isinstance(dtype, Stdout):
            return dtype.subtype
        else:
            return dtype
        
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
        arg = unwrap_expression(
            val=self.src,
            context='workflow',
            scope=self.scope, 
        )
        if isinstance(arg, list):
            raise NotImplementedError
            call_args += arg

        # handle misc edge case (takes priority over datatype mismatches)
        if satisfies_edge_case(self.src):
            suffix = handle_edge_case(self.src)
            arg = f'{arg}{suffix}'

        # handle datatype relationship
        elif is_datatype_mismatch(self.srctype, self.desttype, self.dest_scatter):
            suffix = gen_datatype_mismatch_plumbing(self.srctype, self.desttype, self.dest_scatter)
            arg = f'{arg}{suffix}'
        
        return arg

    




