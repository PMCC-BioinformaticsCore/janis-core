
from copy import deepcopy
from typing import Any, Optional

import janis_core.translation_utils as utils
from janis_core.workflow.workflow import Workflow, InputNode, StepNode
from janis_core import CommandTool, PythonTool, TInput
from janis_core import settings

from ..params_channels.helpers_common import get_all_workflow_inputs

from ...scope import Scope
from ... import data_sources
from ... import params
from ... import channels
from ... import task_inputs


def register_ds_categories(entity: Workflow | CommandTool | PythonTool) -> None:
    """
    MAIN ENTRY POINT for this preprocessing task.

    for each nextflow task to create, decides which ToolInputs are:
        - process inputs
        - accessed internally via param
        - accessed internally via static value
    """
    
    scope = Scope()

    if isinstance(entity, Workflow):
        register_ds_categories_workflow(entity, scope)
    elif isinstance(entity, (CommandTool | PythonTool)):  # type: ignore
        sources: dict[str, Any] = {}
        do_register_ds_categories(entity, sources, scope)
    else:
        raise RuntimeError

def register_ds_categories_workflow(wf: Workflow, scope: Scope) -> None:
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        do_register_ds_categories(step.tool, step.sources, current_scope)
        if isinstance(step.tool, Workflow):
            register_ds_categories_workflow(step.tool, current_scope)

def do_register_ds_categories(tool: Workflow | CommandTool | PythonTool, sources: dict[str, Any], scope: Scope) -> None:
    generator = ProcessDSCategoryGenerator(scope, tool, sources)
    generator.generate()
    data_sources.update_categories(
        scope, 
        generator.process_inputs, 
        generator.param_inputs, 
        generator.internal_inputs
    )


class ProcessDSCategoryGenerator:
    """
    for a nextflow task to create, decides which of the ToolInputs / TInputs which will be:
        - process inputs
        - accessed internally via param
        - accessed internally via static value
    """
    def __init__(self, scope: Scope, tool: Workflow | CommandTool | PythonTool, sources: dict[str, Any]) -> None:
        self.scope = scope
        self.tool = tool
        self.sources = sources
        self.process_inputs: set[str] = set()
        self.param_inputs: set[str] = set()
        self.internal_inputs: set[str] = set()

    # @property
    # def tinputs(self) -> list[TInput]:
    #     if isinstance(self.tool, (CommandTool, Workflow)):
    #         return self.tool.tool_inputs()
    #     elif isinstance(self.tool, PythonTool):  # type: ignore
    #         return self.tool.inputs()
    #     else:  # workflow
    #         raise NotImplementedError
    
    @property
    def all_input_ids(self) -> set[str]:
        if isinstance(self.tool, (CommandTool, PythonTool)):
            assert(self.all_input_ids_tool is not None)
            return self.all_input_ids_tool
        else:
            assert(self.all_input_ids_workflow is not None)
            return self.all_input_ids_workflow
        
    @property
    def all_input_ids_tool(self) -> Optional[set[str]]:
        return {x.id() for x in self.tool.tool_inputs()} 
    
    @property
    def all_input_ids_workflow(self) -> Optional[set[str]]:
        return get_all_workflow_inputs(self.tool)
        
    @property
    def param_input_ids(self) -> set[str]:
        """get tool inputs (ids) which are linked to a param"""
        out: set[str] = set()
        for tag, src in self.sources.items():
            node = utils.resolve_node(src)
            if isinstance(node, InputNode):
                if params.exists(node.uuid):
                    out.add(tag)
        return out
    

    
    # @property
    # def channel_input_ids(self) -> set[str]:
    #     """get tool inputs (ids) which are linked to a channel"""
    #     out: set[str] = set()
    #     for tag, src in self.sources.items():
    #         node = utils.resolve_node(src)
    #         if isinstance(node, InputNode):
    #             if channels.exists(self.scope, node.uuid, for_parent=True):
    #                 out.add(tag)
    #     return out

    # @property
    # def connection_input_ids(self) -> set[str]:
    #     """get tool inputs (ids) which are being fed a value from a step connection"""
    #     out: set[str] = set()
    #     for tag, src in self.sources.items():
    #         node = utils.resolve_node(src)
    #         if isinstance(node, StepNode):
    #             out.add(tag)
    #     return out

    # @property
    # def scatter_input_ids(self) -> set[str]:
    #     """get tool inputs (ids) which are being scattered on"""
    #     out: set[str] = set()
    #     for inname, src in self.sources.items():
    #         should_scatter = src.source_map[0].should_scatter
    #         node = utils.resolve_node(src)
    #         if should_scatter and isinstance(node, InputNode):
    #             out.add(inname)
    #     return out

    # @property
    # def complex_expression_input_ids(self) -> set[str]:
    #     """get tool inputs (ids) which are fed data using a complex janis expression"""
    #     out: set[str] = set()
    #     for inname, src in self.sources.items():
    #         node = utils.resolve_node(src)
    #         if not isinstance(node, InputNode) and not isinstance(node, StepNode):
    #             out.add(inname)
    #     return out
    
    ### main method
    def generate(self) -> None:
        self.process_inputs = self.get_process_inputs()
        self.param_inputs = self.get_param_inputs()
        self.internal_inputs = self.get_internal_inputs()

    ### process inputs
    def get_process_inputs(self) -> set[str]:
        """gets the tool inputs which will become nextflow process inputs"""
        return task_inputs.get(self.tool.id())

    #     if settings.translate.nextflow.MODE == 'workflow':
    #         return self.get_process_inputs_workflowmode()
    #     elif settings.translate.nextflow.MODE == 'tool':  # type: ignore
    #         return self.get_process_inputs_toolmode()
    #     else:
    #         raise RuntimeError

    # def get_process_inputs_workflowmode(self) -> set[str]:
    #     """
    #     inputs which are fed (via step inputs) using a file type workflow input
    #     inputs which are fed (via step inputs) using a connection
    #     inputs which are involved in scatter
    #     """
    #     return self.channel_input_ids | self.connection_input_ids | self.scatter_input_ids | self.complex_expression_input_ids

    # def get_process_inputs_toolmode(self) -> set[str]:
    #     """all CommandTool inputs. in toolmode, we have no greater scope than the tool itself."""
    #     return self.all_input_ids
    
    ### param inputs
    def get_param_inputs(self) -> set[str]:
        """gets the tool inputs which will be fed values via params"""
        if settings.translate.nextflow.MODE == 'workflow':
            return self.get_param_inputs_workflowmode()
        elif settings.translate.nextflow.MODE == 'tool':  # type: ignore
            return self.get_param_inputs_toolmode()
        raise RuntimeError('DEV: settings.translate.nextflow.MODE must be either "workflow" or "tool"')

    def get_param_inputs_workflowmode(self) -> set[str]:
        """get the inputs which have a param, but are not process inputs"""
        return self.param_input_ids - self.process_inputs

    def get_param_inputs_toolmode(self) -> set[str]:
        """no param inputs for toolmode. """
        return set()


    ### internal inputs
    def get_internal_inputs(self) -> set[str]:
        """
        get the tool inputs which will not be exposed to the outside world in the nextflow process
        internal_inputs = all_inputs - process_inputs - param_inputs
        """
        surviving_ids = self.all_input_ids
        surviving_ids = surviving_ids - self.process_inputs
        surviving_ids = surviving_ids - self.param_inputs
        return surviving_ids

