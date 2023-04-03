

from typing import Any
from dataclasses import dataclass, field

from .. import task_inputs

from janis_core import Workflow, Tool
from janis_core import translation_utils as utils
from janis_core.workflow.workflow import InputNode
from janis_core.types import File, Filename, Directory, DataType


def register_minimal_task_inputs(wf: Workflow) -> None:
    collector = TaskInputCollector()
    collector.collect(wf)
    minifier = TaskInputMinifier(collector.catalogue)
    minifier.minify()
    for tool_id, task_input_ids in minifier.task_inputs.items():
        task_inputs.add(tool_id, task_input_ids)


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
        self.values.add(val)
    
    


class TaskInputMinifier:
    def __init__(self, catalogue: dict[str, dict[str, TaskInputHistory]]) -> None:
        self.catalogue = catalogue
        self.task_inputs: dict[str, set[str]] = {}

    def minify(self) -> None:
        for uuid, task_inputs in self.catalogue.items():
            ids_to_remove = self.filter(task_inputs)
            tinput_ids = self.collapse_task_inputs(ids_to_remove, task_inputs)
            self.task_inputs[uuid] = tinput_ids

    def filter(self, task_inputs: dict[str, TaskInputHistory]) -> set[str]:
        ids_to_remove: set[str] = set()

        for tinput_id, history in task_inputs.items():
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
            
    def collapse_task_inputs(self, ids_to_remove: set[str], task_inputs: dict[str, TaskInputHistory]) -> set[str]:
        all_input_ids = set(task_inputs.keys())
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
            task_inputs = {tinput_id: None for tinput_id in step.tool.inputs_map().keys()}
            
            # update task inputs for tinputs with item in sources
            task_inputs = self.replace_sources(step.sources, task_inputs)

            # trace each value in task_inputs, if is InputNode with default and is not file type, 
            # update the step inputs with the default.
            # this is an example of a static value (eg inStr='hello')
            task_inputs = self.replace_static_values(task_inputs)
            
            # update the dict
            self.update_catalogue(identifier, task_inputs, step.tool)

            # recursive for subworkflows
            if isinstance(step.tool, Workflow):
                self.collect(step.tool)

    def replace_sources(self, sources: dict[str, Any], task_inputs: dict[str, Any]) -> dict[str, Any]:
        for tinput_id in task_inputs.keys():
            if tinput_id in sources:
                src = sources[tinput_id]
                node = utils.resolve_node(src)
                task_inputs[tinput_id] = node
        return task_inputs
    
    def replace_static_values(self, task_inputs: dict[str, Any]) -> dict[str, Any]:
        for tid, src in task_inputs.items():
            node = utils.resolve_node(src)
            if isinstance(node, InputNode):
                if node.default is not None:  # type: ignore
                    task_inputs[tid] = node.default  # type: ignore
        return task_inputs
    
    def update_catalogue(self, identifier: str, task_inputs: dict[str, Any], tool: Tool) -> None:
        if identifier not in self.catalogue:
            self.catalogue[identifier] = {}

        for tinput_id, val in task_inputs.items():
            if tinput_id not in self.catalogue[identifier]:
                tinput = [x for x in tool.tool_inputs() if x.id() == tinput_id][0]
                history = TaskInputHistory(tinput_id, tinput.intype)
                self.catalogue[identifier][tinput_id] = history
            
            self.catalogue[identifier][tinput_id].add_value(val)





