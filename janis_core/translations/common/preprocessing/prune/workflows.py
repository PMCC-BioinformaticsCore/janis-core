

from dataclasses import dataclass, field

from janis_core import WorkflowBuilder, CommandToolBuilder, PythonTool
from janis_core.types import DataType, File
from janis_core.workflow.workflow import InputNode
from janis_core.operators import InputNodeSelector
from janis_core.translations.common import trace
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType


def prune_main_workflow_inputs(wf: WorkflowBuilder) -> WorkflowBuilder:
    valid_workflow_inputs: set[str] = set()
    valid_workflow_inputs = valid_workflow_inputs | get_mandatory_input_ids(wf)
    valid_workflow_inputs = valid_workflow_inputs | get_referenced_input_ids(wf)
    wf.input_nodes = {k: v for k, v in wf.input_nodes.items() if k in valid_workflow_inputs}
    return wf

@dataclass 
class InputNodeReference:
    step_id: str 
    tinput_id: str
    tinput_dtype: DataType

@dataclass 
class InputNodeReferenceCatalogue:
    node: InputNode
    references: list[InputNodeReference] = field(default_factory=list)
    
    def add(self, step_id: str, tinput_id: str, tinput_dtype: DataType):
        self.references.append(InputNodeReference(step_id, tinput_id, tinput_dtype))

def get_mandatory_input_ids(wf: WorkflowBuilder) -> set[str]:
    collected_ids: set[str] = set()
    
    for tinput_id, tinput in wf.input_nodes.items():
        if tinput.datatype.optional == False:
            collected_ids.add(tinput_id)
    
    return collected_ids

def get_referenced_input_ids(wf: WorkflowBuilder) -> set[str]:
    collected_ids: set[str] = set()
    
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            selector = src.source_map[0].source
            if isinstance(selector, InputNodeSelector):
                collected_ids.add(selector.input_node.id())
    
    return collected_ids

# def get_placeholder_inputs(wf: WorkflowBuilder) -> set[str]:
#     collected_ids: set[str] = set()
    
#     for tinput_id, input_node in wf.input_nodes.items():
#         catalogue = get_reference_catalogue(wf, input_node)
#         if seems_unused(catalogue):
#             collected_ids.add(tinput_id)
#         elif seems_like_placeholder(catalogue):
#             collected_ids.add(tinput_id)
    
#     return collected_ids

# def get_reference_catalogue(wf: WorkflowBuilder, query_node: InputNode) -> InputNodeReferenceCatalogue:
#     reference_store = InputNodeReferenceCatalogue(query_node)

#     # iterate through each step_input for each step in workflow
#     for step in wf.step_nodes.values():
#         for tinput_id, src in step.sources.items():
#             referenced_vars = trace.trace_referenced_variables(src, wf)
#             for var in referenced_vars:
#                 if var == query_node.id():
#                     if isinstance(step.tool, CommandToolBuilder):
#                         tinput_dtype = [x for x in step.tool.inputs() if x.id() == tinput_id][0].input_type
#                     elif isinstance(step.tool, PythonTool):
#                         tinput_dtype = [x for x in step.tool.inputs() if x.id() == tinput_id][0].intype
#                     elif isinstance(step.tool, WorkflowBuilder):
#                         tinput_dtype = [x for x in step.tool.input_nodes.values() if x.id() == tinput_id][0].datatype
#                     else:
#                         raise RuntimeError
                    
#                     reference_store.add(step.id(), tinput_id, tinput_dtype)
    
#     return reference_store

# def seems_unused(catalogue: InputNodeReferenceCatalogue) -> bool:
#     if len(catalogue.references) == 0:
#         if catalogue.node.datatype.optional == True:
#             return True
#     return False

# def seems_like_placeholder(catalogue: InputNodeReferenceCatalogue) -> bool:
#     # single reference with dummy format name implies placeholder
#     if len(catalogue.references) == 1:
#         node = catalogue.node
#         ref = catalogue.references[0]
#         step_id = ref.step_id
#         tinput_id = ref.tinput_id
#         tinput_dtype = ref.tinput_dtype

#         if utils.looks_like_placeholder_node(node, step_id, tinput_id, tinput_dtype):
#             return True
    
#     return False