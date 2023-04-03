

from typing import Any
from dataclasses import dataclass, field

from .. import task_inputs

from janis_core import Workflow, Tool, CommandTool, PythonTool
from janis_core import translation_utils as utils
from janis_core.workflow.workflow import InputNode
from janis_core.types import File, Filename, Directory, DataType


def register_minimal_task_inputs(tool: Workflow | CommandTool | PythonTool) -> None:
    if isinstance(tool, Workflow):
        do_register_minimal_task_inputs_wf(tool)
    else:
        do_register_minimal_task_inputs_tool(tool)

def do_register_minimal_task_inputs_wf(wf: Workflow) -> None:
    collector = TaskInputCollector()
    collector.collect(wf)
    minifier = TaskInputMinifier(collector.catalogue)
    minifier.minify()
    for tool_id, tinput_ids in minifier.task_inputs_dict.items():
        task_inputs.add(tool_id, tinput_ids)

def do_register_minimal_task_inputs_tool(tool: CommandTool | PythonTool) -> None:
    tinput_ids: set[str] = set()
    for tinput in tool.tool_inputs():
        if not tinput.intype.optional:
            tinput_ids.add(tinput.id())
    task_inputs.add(tool.id(), tinput_ids)


@dataclass
class TaskInputHistory:
    tinput_id: str
    dtype: DataType
    values: set[Any] = field(default_factory=set)

    @property
    def is_file(self) -> bool:
        basetype = utils.get_base_type(self.dtype)
        basetype = utils.ensure_single_type(basetype)
        if isinstance(basetype, (File, Filename, Directory)):
            return True
        return False
    
    @property
    def is_optional(self) -> bool:
        if self.dtype.optional == True:
            return True
        return False
    
    @property
    def non_null_values(self) -> set[Any]:
        return set([x for x in self.values if x is not None])
    
    def add_value(self, val: Any) -> None:
        if isinstance(val, list):
            val = tuple(val)
        self.values.add(val)


class TaskInputMinifier:
    def __init__(self, catalogue: dict[str, dict[str, TaskInputHistory]]) -> None:
        self.catalogue = catalogue
        self.task_inputs_dict: dict[str, set[str]] = {}

    def minify(self) -> None:
        for uuid, input_dict in self.catalogue.items():
            ids_to_remove = self.filter(input_dict)
            tinput_ids = self.collapse_task_inputs(ids_to_remove, input_dict)
            self.task_inputs_dict[uuid] = tinput_ids

    def filter(self, input_dict: dict[str, TaskInputHistory]) -> set[str]:
        ids_to_remove: set[str] = set()

        for tinput_id, history in input_dict.items():
            # non-optional files must stay (mandatory task inputs)
            if history.is_file and not history.is_optional:  
                continue

            elif history.is_file:
                if len(history.non_null_values) == 0:  
                    ids_to_remove.add(tinput_id)

            # optional
            elif len(history.non_null_values) <= 1:
                ids_to_remove.add(tinput_id)
        
        return ids_to_remove
            
    def collapse_task_inputs(self, ids_to_remove: set[str], input_dict: dict[str, TaskInputHistory]) -> set[str]:
        all_input_ids = set(input_dict.keys())
        collapsed_input_ids = all_input_ids - ids_to_remove
        return collapsed_input_ids



class TaskInputCollector:
    def __init__(self) -> None:
        self.catalogue: dict[str, dict[str, TaskInputHistory]] = {}

    def collect(self, wf: Workflow) -> None:
        for step in wf.step_nodes.values():
            # the uuid is a unique identifier for this Tool. use as key. 
            identifier = step.tool.id()

            # set all the tool inputs to None
            inputs_dict = {tinput_id: None for tinput_id in step.tool.inputs_map().keys()}
            
            # update task inputs for tinputs with item in sources
            inputs_dict = self.replace_sources(step.sources, inputs_dict)

            # trace each value in task_inputs, if is InputNode with default and is not file type, 
            # update the step inputs with the default.
            # this is an example of a static value (eg inStr='hello')
            inputs_dict = self.replace_static_values(inputs_dict)
            
            # update the dict
            self.update_catalogue(identifier, inputs_dict, step.tool)

            # recursive for subworkflows
            if isinstance(step.tool, Workflow):
                self.collect(step.tool)

    def replace_sources(self, sources: dict[str, Any], inputs_dict: dict[str, Any]) -> dict[str, Any]:
        for tinput_id in inputs_dict.keys():
            if tinput_id in sources:
                src = sources[tinput_id]
                node = utils.resolve_node(src)
                inputs_dict[tinput_id] = node
        return inputs_dict
    
    def replace_static_values(self, inputs_dict: dict[str, Any]) -> dict[str, Any]:
        for tid, src in inputs_dict.items():
            node = utils.resolve_node(src)
            if isinstance(node, InputNode):
                if node.default is not None:  # type: ignore
                    inputs_dict[tid] = node.default  # type: ignore
        return inputs_dict
    
    def update_catalogue(self, identifier: str, inputs_dict: dict[str, Any], tool: Tool) -> None:
        if identifier not in self.catalogue:
            self.catalogue[identifier] = {}

        for tinput_id, val in inputs_dict.items():
            if tinput_id not in self.catalogue[identifier]:
                tinput = [x for x in tool.tool_inputs() if x.id() == tinput_id][0]
                history = TaskInputHistory(tinput_id, tinput.intype)
                self.catalogue[identifier][tinput_id] = history
            
            self.catalogue[identifier][tinput_id].add_value(val)





