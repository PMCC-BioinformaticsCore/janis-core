


# from typing import Optional, Any
# from copy import deepcopy

# from janis_core import (
#     DataType,
#     TInput,
#     Workflow,
#     CommandTool, 
#     PythonTool, 
# )

# from janis_core.workflow.workflow import Workflow, InputNode
# from janis_core import translation_utils as utils
# from janis_core import settings

# from .. import task_inputs
# from .. import params
# from .. import naming 

# from ..scope import Scope
# from ..casefmt import to_case

# from .common import get_all_workflow_inputs


# def register_data_sources(entity: Workflow | CommandTool | PythonTool) -> None:
#     """
#     MAIN ENTRY POINT for this preprocessing task.

#     for each CommandTool / PythonTool, decides the variable_name (a task_input, a param_input, or None ) 
#     which feeds data for the tinput.
#     the variable_name is how the particular tinput_id will be referenced inside its process. 
#     """
    
#     scope = Scope()

#     if isinstance(entity, Workflow):
#         register_ds_variables_workflow(scope, entity)
#     elif isinstance(entity, (CommandTool | PythonTool)):  # type: ignore
#         sources: dict[str, Any] = {}
#         register_ds_variables_tool(scope, entity, sources)
#     else:
#         raise RuntimeError

# def register_ds_variables_workflow(scope: Scope, wf: Workflow) -> None:
#     for step in wf.step_nodes.values():
#         current_scope = deepcopy(scope)
#         current_scope.update(step)

#         register_ds_variables_tool(current_scope, step.tool, step.sources)
#         if isinstance(step.tool, Workflow):
#             register_ds_variables_workflow(current_scope, step.tool)

# def register_ds_variables_tool(scope: Scope, tool: CommandTool | PythonTool, sources: dict[str, Any]) -> None:
#     registerer = DataSourceRegistrationManager(scope, tool, sources)
#     registerer.register()


# class DataSourceRegistrationManager:
#     """
#     for each CommandTool / PythonTool, determines the internal variable_name for each tinput.id()
#     the variable_name is how the particular tinput_id will be referenced inside the process. 
#     """
#     def __init__(self, scope: Scope, tool: Workflow | CommandTool | PythonTool, sources: dict[str, Any]) -> None:
#         self.scope = scope
#         self.tool = tool
#         self.sources = sources

#     ### main method
#     def register(self) -> None:
#         for tinput in self.tinputs:
#             self.update_data_sources(tinput)
    
#     def update_data_sources(self, inp: TInput) -> None:
#         if inp.id() in self.task_inputs:
#             return self.update_dss_task_input(inp)
#         elif inp.id() in self.param_inputs:
#             return self.update_dss_param_input(inp)
#         elif inp.id() in self.internal_inputs:
#             return self.update_dss_internal_input(inp)
#         else:
#             pass
        
#     def update_dss_task_input(self, inp: TInput) -> None:
#         """TASK_INPUT"""
#         dtype: DataType = inp.intype  # type: ignore
#         is_duplicate = self.duplicate_datatype_exists(inp)
        
#         if utils.is_array_secondary_type(dtype):
#             name = naming.process.secondaries_array(inp, duplicate_datatype_exists=is_duplicate)
#         elif utils.is_secondary_type(dtype):
#             name = naming.process.secondaries(inp, duplicate_datatype_exists=is_duplicate)
#         else:
#             name = naming.process.generic(inp)  
#         data_sources.update(self.scope, dstype_str='task_input', tinput_id=inp.id(), value=name)
        
#     def update_dss_param_input(self, inp: TInput) -> None: 
#         """PARAM"""
#         src = self.sources[inp.id()]
#         sel = src.source_map[0].source
#         param = params.get(sel.input_node.uuid)
#         pname = to_case(param.name, case=settings.translate.nextflow.NF_PROCESS_INPUT_CASE)
#         pname = f'params.{pname}'
#         data_sources.update(self.scope, dstype_str='param', tinput_id=inp.id(), value=pname)
    
#     def update_dss_internal_input(self, inp: TInput) -> None:
#         """STATIC, IGNORED"""
#         value = None
        
#         # update value if InputNode with default 
#         if inp.id() in self.sources:
#             src = self.sources[inp.id()]
#             node = utils.resolve_node(src)
#             if isinstance(node, InputNode):
#                 value = node.default
        
#         # is this something with static value (inStr='hello'), 
#         # or ignored (tinput doesnt appear in janis step call, or is None in step call)?
#         if value is not None:
#             dstype_str = 'static'
#         else:
#             dstype_str = 'ignored'

#         data_sources.update(self.scope, dstype_str=dstype_str, tinput_id=inp.id(), value=value)

#     ### helper properties
#     @property
#     def tinputs(self) -> list[TInput]:
#         return self.tool.tool_inputs()
        
#     @property
#     def task_inputs(self) -> set[str]:
#         """the tool inputs which will become nextflow process inputs"""
#         return task_inputs.get(self.tool.id())
    
#     @property
#     def param_inputs(self) -> set[str]:
#         """gets the tool inputs which will be fed values via params"""
#         if settings.translate.nextflow.MODE == 'workflow':
#             return self.param_inputs_workflowmode
#         elif settings.translate.nextflow.MODE == 'tool':  # type: ignore
#             return self.param_inputs_toolmode
#         raise RuntimeError('DEV: settings.translate.nextflow.MODE must be either "workflow" or "tool"')

#     @property
#     def param_inputs_workflowmode(self) -> set[str]:
#         """get the inputs which have a param, but are not process inputs"""
#         if not isinstance(self.tool, Workflow):
#             return set()
        
#         # get param inputs
#         out: set[str] = set()
#         for tag, src in self.sources.items():
#             node = utils.resolve_node(src)
#             if isinstance(node, InputNode):
#                 if params.exists(node.uuid):
#                     out.add(tag)
        
#         # remove process inputs
#         out = out - self.task_inputs
#         return out

#     @property
#     def param_inputs_toolmode(self) -> set[str]:
#         """no param inputs for toolmode. """
#         return set()

#     @property
#     def internal_inputs(self) -> set[str]:
#         """
#         get the tool inputs which will not be exposed to the outside world in the nextflow process
#         internal_inputs = all_inputs - task_inputs - param_inputs
#         """
#         all_ids = set([x.id() for x in self.tinputs])
#         surviving_ids = all_ids
#         surviving_ids = surviving_ids - self.task_inputs
#         surviving_ids = surviving_ids - self.param_inputs
#         return surviving_ids
    
#     # @property
#     # def all_input_ids(self) -> set[str]:
#     #     if isinstance(self.tool, (CommandTool, PythonTool)):
#     #         assert(self._all_input_ids_tool is not None)
#     #         return self._all_input_ids_tool
#     #     else:
#     #         assert(self._all_input_ids_workflow is not None)
#     #         return self._all_input_ids_workflow
        
#     # @property
#     # def _all_input_ids_tool(self) -> Optional[set[str]]:
#     #     return {x.id() for x in self.tinputs} 
    
#     # @property
#     # def _all_input_ids_workflow(self) -> Optional[set[str]]:
#     #     assert(isinstance(self.tool, Workflow))
#     #     return get_all_workflow_inputs(self.tool)


#     ### helper methods
#     def duplicate_datatype_exists(self, inp: TInput) -> bool:
#         """
#         check if another TInput has the same dtype as this TInput.
#         only checking for secondary and array secondary datatype duplicates.
#         """
#         rtype = inp.intype  # type: ignore
#         rbasetype = utils.get_base_type(rtype)  # type: ignore
        
#         for tinput in self.tinputs:
#             # dont check tinput against itself
#             if tinput.id() == inp.id():
#                 continue
            
#             # only check against other task inputs
#             if tinput.id() in self.task_inputs:
#                 # check if types match
#                 qtype = tinput.intype  # type: ignore
#                 qbasetype = utils.get_base_type(qtype)  # type: ignore

#                 if utils.is_array_secondary_type(rtype) and utils.is_array_secondary_type(qtype):  # type: ignore
#                     if type(rbasetype) == type(qbasetype):
#                         return True
                
#                 elif utils.is_secondary_type(rtype) and utils.is_secondary_type(qtype):  # type: ignore
#                     if type(rtype) == type(qtype):  # type: ignore
#                         return True
#         return False
               