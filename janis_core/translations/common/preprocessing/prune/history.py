

from typing import Any
from dataclasses import dataclass, field
from copy import deepcopy

from janis_core import Tool, WorkflowBuilder, PythonTool, TInput, CommandToolBuilder
from janis_core import translation_utils as utils
from janis_core.workflow.workflow import InputNode
from janis_core.operators.selectors import InputNodeSelector
from janis_core.operators.selectors import StepOutputSelector
from janis_core.types import DataType
from janis_core.translations.common import trace



@dataclass
class TaskInputHistory:
    tinput: TInput
    sources: list[InputNodeSelector | StepOutputSelector] = field(default_factory=list)

    @property
    def dtype(self) -> DataType:
        return self.tinput.intype

    @property
    def is_file(self) -> bool:
        return utils.is_file_type(self.dtype)
    
    @property
    def is_optional(self) -> bool:
        if self.dtype.optional == True:
            return True
        return False
    
    @property
    def connections(self) -> bool:
        for src in self.sources:
            if isinstance(src, StepOutputSelector):
                return True
        return False
    
    @property
    def mandatory_input_sources(self) -> set[InputNodeSelector]:
        out: set[InputNodeSelector] = set()
        for src in self.sources:
            if isinstance(src, InputNodeSelector):
                if src.input_node.datatype.optional == False:
                    out.add(src)
        return out
    
    def add_value(self, src: InputNodeSelector | StepOutputSelector) -> None:
        self.sources.append(src)


@dataclass 
class InputNodeReference:
    step_id: str 
    tinput_id: str

@dataclass 
class InputNodeReferenceStore:
    node_id: str
    references: list[InputNodeReference] = field(default_factory=list)

    @property 
    def seems_like_dummy(self) -> bool:
        # only 1 reference 
        if len(self.references) == 1:
            # the 1 reference implies dummy format
            ref = self.references[0]
            if self.node_id == f'{ref.step_id}_{ref.tinput_id}':
                return True
        return False


class TaskInputCollector:
    """
    for a given tool_id, searches the workflow for each step calling that tool.
    records the values provided to each TInput in that step call. 
    """
    def __init__(self, tool: Tool) -> None:
        self.tool = tool
        self.tool_id = tool.id()
        self.histories: dict[str, TaskInputHistory] = {}
        self.step_count: int = 0

    @property
    def base_inputs_dict(self) -> dict[str, Any]:
        if isinstance(self.tool, CommandToolBuilder) or isinstance(self.tool, PythonTool):
            return self.init_inputs_dict_tool()
        else:
            return self.init_inputs_dict_workflow()

    def get_input_node_references(self, entity: Any, tool: Tool) -> list[str]:
        referenced_vars = trace.trace_referenced_variables(entity, tool)
        all_input_keys = tool.inputs_map().keys()
        input_vars = [x for x in referenced_vars if x in all_input_keys]
        return input_vars
            
    def init_inputs_dict_tool(self) -> dict[str, Any]:
        return {tinput.id(): None for tinput in self.tool.tool_inputs()}

    def init_inputs_dict_workflow(self) -> dict[str, Any]:
        """
        we assume all workflow inputs are valid, except those which are not files 
        and look like they are dummy inputs derived when a task has a static input value. 

        a dummy input has the name formatted as <step_id>_<tinput_id>,
        where tinput id is the id of the input in the tool being called.

        to find these, for each workflow input, we assume it is valid. 
        we then search through the step inputs of each step in the workflow for that input.

        if we find only 1 reference and it has the format of a dummy input, we conclude it 
        is likely to be one of these and remove it from the inputs dict.
        """
        assert(isinstance(self.tool, WorkflowBuilder))
        valid_inputs: set[str] = set()

        for node in self.tool.input_nodes.values():
            
            # for non files check if it is a dummy input
            if not utils.is_file_type(node.datatype):
                reference_store = self.gen_reference_store(node)
                if reference_store.seems_like_dummy:
                    continue
            
            valid_inputs.add(node.id())
        
        return {x: None for x in valid_inputs}

    def gen_reference_store(self, node: InputNode) -> InputNodeReferenceStore:
        assert(isinstance(self.tool, WorkflowBuilder))
        reference_store = InputNodeReferenceStore(node.id())

        # iterate through each step_input for each step in workflow
        for step in self.tool.step_nodes.values():
            for tinput_id, src in step.sources.items():
                node_names = self.get_input_node_references(src, self.tool)
                for name in node_names:
                    if name == node.id():
                        new_ref = InputNodeReference(step.id(), tinput_id)
                        reference_store.references.append(new_ref)
        
        return reference_store

    def collect(self, wf: WorkflowBuilder) -> None:
        for step in wf.step_nodes.values():
            # the uuid is a unique identifier for this Tool. use as key. 
            identifier = step.tool.id()
            if identifier == self.tool_id:
                self.step_count += 1
                # initialise inputs dict
                inputs_dict = deepcopy(self.base_inputs_dict)
                
                # update task inputs for tinputs with item in sources
                inputs_dict = self.replace_sources(step.sources, inputs_dict)

                # trace each value in task_inputs, if is InputNode with default and is not file type, 
                # update the step inputs with the default.
                # this is an example of a static value (eg inStr='hello')
                # inputs_dict = self.replace_static_values(inputs_dict)
                
                # update the dict
                self.update_histories(inputs_dict, step.tool)

            # recursive for subworkflows
            if isinstance(step.tool, WorkflowBuilder):
                self.collect(step.tool)
    
    def replace_sources(self, sources: dict[str, Any], inputs_dict: dict[str, Any]) -> dict[str, InputNodeSelector | StepOutputSelector]:
        for tinput_id in inputs_dict.keys():
            if tinput_id in sources:
                selector = sources[tinput_id].source_map[0].source
                inputs_dict[tinput_id] = selector
        return inputs_dict
    
    # def replace_static_values(self, inputs_dict: dict[str, Any]) -> dict[str, Any]:
    #     for tid, node in inputs_dict.items():
    #         if isinstance(node, InputNode):
    #             if node.default is not None:  # type: ignore
    #                 inputs_dict[tid] = node.default  # type: ignore
    #     return inputs_dict
    
    def update_histories(self, inputs_dict: dict[str, Any], tool: Tool) -> None:
        for tinput_id, src in inputs_dict.items():

            # add a TaskInputHistory for TInput if not exists
            if tinput_id not in self.histories:
                if tinput_id == 'domfiles':
                    print()
                tinput = [x for x in tool.tool_inputs() if x.id() == tinput_id][0]
                history = TaskInputHistory(tinput)
                self.histories[tinput_id] = history
            
            # add a value to this TInput's TaskInputHistory
            if src is not None:
                self.histories[tinput_id].add_value(src)

