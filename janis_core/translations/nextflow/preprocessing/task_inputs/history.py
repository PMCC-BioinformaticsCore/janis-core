

from typing import Any, Optional
from dataclasses import dataclass, field
from copy import deepcopy

from janis_core import Tool, Workflow, CommandTool, PythonTool, TInput, TOutput
from janis_core import translation_utils as utils
from janis_core.workflow.workflow import InputNode
from janis_core.types import DataType

from ... import trace


@dataclass
class TaskInputHistory:
    tinput_id: str
    dtype: DataType
    values: list[Optional[str]] = field(default_factory=list)

    @property
    def is_file(self) -> bool:
        return utils.is_file_type(self.dtype)
    
    @property
    def is_optional(self) -> bool:
        if self.dtype.optional == True:
            return True
        return False
    
    @property
    def unique_values(self) -> set[Optional[str]]:
        return set([x for x in self.values])
    
    @property
    def non_null_unique_values(self) -> set[str]:
        return set([x for x in self.values if x is not None])
    
    def add_value(self, val: Any) -> None:
        if val is not None:
            val = str(val)
        self.values.append(val)



class TaskInputCollector:
    """
    for a given tool_id, searches the workflow for each step calling that tool.
    records the values provided to each TInput in that step call. 
    """
    def __init__(self, tool: Tool) -> None:
        self.tool = tool
        self.tool_id = tool.id()
        self.histories: dict[str, TaskInputHistory] = {}

    @property
    def base_inputs_dict(self) -> dict[str, Any]:
        if isinstance(self.tool, CommandTool) or isinstance(self.tool, PythonTool):
            return self.init_inputs_dict_tool()
        else:
            return self.init_inputs_dict_workflow()

    def trace_to_node(self, entity: Any, tool: Tool) -> Optional[TInput | TOutput]:
        referenced_vars = trace.trace_referenced_variables(entity, tool)
        for var in referenced_vars:
            if var in tool.all_input_keys():
                return tool.inputs_map()[var]
            else:
                print()
            
    def init_inputs_dict_tool(self) -> dict[str, Any]:
        return {tinput_id: None for tinput_id in self.tool.inputs_map().keys()}

    def init_inputs_dict_workflow(self) -> dict[str, Any]:
        assert(isinstance(self.tool, Workflow))
        inputs_dict: dict[str, Any] = {}
        
        for node in self.tool.input_nodes.values():
            is_valid_wf_input = True
            
            for step in self.tool.step_nodes.values():
                for tinput_id, src in step.sources.items():
                    src_node = utils.resolve_node(src)
                    # src_node = self.trace_to_node(src, self.tool)
                    if isinstance(src_node, InputNode) and src_node.id() == node.id():
                        derived_fmt = f'{step.id()}_{tinput_id}'
                        if node.id() == derived_fmt:
                            is_valid_wf_input = False

            if is_valid_wf_input:
                inputs_dict[node.id()] = None
        
        return inputs_dict

    def collect(self, wf: Workflow) -> None:
        for step in wf.step_nodes.values():
            # the uuid is a unique identifier for this Tool. use as key. 
            identifier = step.tool.id()
            if identifier == self.tool_id:
                # initialise inputs dict
                inputs_dict = deepcopy(self.base_inputs_dict)
                
                # update task inputs for tinputs with item in sources
                inputs_dict = self.replace_sources(step.sources, inputs_dict)

                # trace each value in task_inputs, if is InputNode with default and is not file type, 
                # update the step inputs with the default.
                # this is an example of a static value (eg inStr='hello')
                inputs_dict = self.replace_static_values(inputs_dict)
                
                # update the dict
                self.update_histories(inputs_dict, step.tool)

            # recursive for subworkflows
            if isinstance(step.tool, Workflow):
                self.collect(step.tool)
    
    def replace_sources(self, sources: dict[str, Any], inputs_dict: dict[str, Any]) -> dict[str, Any]:
        for tinput_id in inputs_dict.keys():
            if tinput_id in sources:
                src = sources[tinput_id]
                node = utils.resolve_node(src)
                # node = self.trace_to_node(src, self.tool)
                inputs_dict[tinput_id] = node
        return inputs_dict
    
    def replace_static_values(self, inputs_dict: dict[str, Any]) -> dict[str, Any]:
        for tid, src in inputs_dict.items():
            node = utils.resolve_node(src)
            # node = self.trace_to_node(src, self.tool)
            if isinstance(node, InputNode):
                if node.default is not None:  # type: ignore
                    inputs_dict[tid] = node.default  # type: ignore
        return inputs_dict
    
    def update_histories(self, inputs_dict: dict[str, Any], tool: Tool) -> None:
        for tinput_id, val in inputs_dict.items():

            # add a TaskInputHistory for TInput if not exists
            if tinput_id not in self.histories:
                tinput = [x for x in tool.tool_inputs() if x.id() == tinput_id][0]
                history = TaskInputHistory(tinput_id, tinput.intype)
                self.histories[tinput_id] = history
            
            # add a value to this TInput's TaskInputHistory
            self.histories[tinput_id].add_value(val)

