

# from copy import deepcopy
# from typing import Any

# from janis_core.workflow.workflow import Workflow, InputNode
# from janis_core import translation_utils as utils

# # from ... import params
# # from ... import task_inputs

# # from ..common import get_file_wf_inputs
# # from ..common import get_scatter_wf_inputs
# # from ..common import get_all_workflow_inputs

# from .ChannelRegistrationManager import ChannelRegistrationManager
# from .ParamRegistrationManager import ParamRegistrationManager


# def register_params(wf: Workflow) -> None:
#     """
#     register params & channels for workflow inputs. 
#     """
#     scope = Scope()
#     sources = {}
#     do_register_params_channels(wf, sources, scope)

# def do_register_params_channels(wf: Workflow, sources: dict[str, Any], scope: Scope):
    
#     # link wf InputNodes to params where driven by param in enclosing scope
#     # makes params global
#     if len(scope.labels) > 1:
#         linkable_params = get_linkable_params(wf, sources)
#         for uuid, param in linkable_params.items():
#             params.create_link(uuid, param)

#     param_inputs = get_param_inputs_to_register(wf, scope)
#     channel_inputs = get_channel_inputs_to_register(wf, scope)

#     for inp in wf.input_nodes.values():
#         # params go first
#         if inp.id() in param_inputs:
#             ParamRegistrationManager(inp, scope).register()
#         if inp.id() in channel_inputs:
#             ChannelRegistrationManager(inp, scope).register()
    
#     # repeat for nested workflows (subworkflows)
#     for step in wf.step_nodes.values():
#         current_scope = deepcopy(scope)
#         current_scope.update(step)
#         if isinstance(step.tool, Workflow):
#             do_register_params_channels(step.tool, step.sources, current_scope)


# ### params
# def get_linkable_params(wf: Workflow, sources: dict[str, Any]) -> dict[str, params.Param]:
#     all_inputs = get_all_workflow_inputs(wf)
#     input_nodes = [wf.input_nodes[x] for x in all_inputs]

#     out: dict[str, params.Param] = {}
#     for inp in input_nodes:
#         if inp.id() not in sources:
#             continue
#         else:
#             src = sources[inp.id()]
#             node = utils.resolve_node(src)
#             if isinstance(node, InputNode):
#                 if params.exists(node.uuid):
#                     param = params.get(node.uuid)
#                     out[inp.uuid] = param
#     return out

# def get_param_inputs_to_register(wf: Workflow, scope: Scope) -> set[str]:
#     if len(scope.labels) > 1:
#         return get_param_inputs_to_register_sub(wf)
#     else:
#         return get_param_inputs_to_register_main(wf)

# def get_param_inputs_to_register_main(wf: Workflow) -> set[str]:
#     return get_all_workflow_inputs(wf)

# def get_param_inputs_to_register_sub(wf: Workflow) -> set[str]:
#     """
#     which sub workflow InputNodes qualify as params we should register?
#     in sub workflow, every data source gets a param EXCEPT:
#     - if the data source will be a nextflow channel
#     - if the data source can be backtraced to a param which already exists. 
#     - if the data source is static (ie 1 or 'conservative' etc)
    
#     eg we have a global param, which is accessed within a subworkflow.
#     """
#     all_inputs = get_all_workflow_inputs(wf)
#     channel_inputs = get_channel_inputs_to_register_sub(wf)
#     global_param_inputs = get_global_param_inputs(wf)
    
#     surviving_inputs = all_inputs
#     surviving_inputs = surviving_inputs - channel_inputs
#     surviving_inputs = surviving_inputs - global_param_inputs
#     return surviving_inputs

# def get_global_param_inputs(wf: Workflow) -> set[str]:
#     all_inputs = get_all_workflow_inputs(wf)
#     input_nodes = [wf.input_nodes[x] for x in all_inputs]

#     out: set[str] = set()
#     for inp in input_nodes:
#         if params.exists(inp.uuid):
#             out.add(inp.id())
#     return out



# ### channels
# def get_channel_inputs_to_register(wf: Workflow, scope: Scope) -> set[str]:
#     if len(scope.labels) > 1:
#         items: set[str] = get_channel_inputs_to_register_sub(wf)
#     else:
#         items: set[str] = get_channel_inputs_to_register_main(wf)
#     return items

# def get_channel_inputs_to_register_main(wf: Workflow) -> set[str]:
#     """
#     Get the wf inputs for which we will create a nf channel in the main wf.
#     inputs which are:
#     - files
#     - filenames (in some cases)
#     - scattered on 
#     will become nextflow channels. all other inputs will become global scope params. 
#     """
#     all_inputs = get_all_workflow_inputs(wf)
#     file_inputs = all_inputs & get_file_wf_inputs(wf)
#     scatter_inputs = all_inputs & get_scatter_wf_inputs(wf)

#     final_inputs = file_inputs | scatter_inputs
#     return final_inputs

# def get_channel_inputs_to_register_sub(wf: Workflow) -> set[str]:
#     """
#     Get the wf inputs for which we will create a nf channel in a sub wf.
#     inputs which:
#     - are fed from a channel
#     - are fed from a step output
#     will be channels. others will be params, or have a static value which will be used in the subworkflow. 
#     """
#     return task_inputs.get(wf.id())




# #     channel_inputs = get_channel_inputs_sub(wf, sources, scope)
# #     connection_inputs = get_connection_inputs_sub(wf, sources)
# #     final_inputs = channel_inputs | connection_inputs
# #     return final_inputs

# # def get_channel_inputs_sub(wf: Workflow, sources: dict[str, Any], scope: Scope) -> set[str]:
# #     out: set[str] = set()
# #     for name, inp in sources.items():
# #         node = utils.resolve_node(inp)
# #         if isinstance(node, InputNode):
# #             if channels.exists(scope, node.uuid, for_parent=True): 
# #                 out.add(name)
# #     return out

# # def get_connection_inputs_sub(wf: Workflow, sources: dict[str, Any]) -> set[str]:
# #     out: set[str] = set()
# #     for name, inp in sources.items():
# #         node = utils.resolve_node(inp)
# #         if isinstance(node, StepNode):
# #             out.add(name)
# #     return out