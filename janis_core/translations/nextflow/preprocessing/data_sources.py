
from copy import deepcopy
from typing import Any

import janis_core.translation_utils as utils
from janis_core.workflow.workflow import Workflow, InputNode, StepNode
from janis_core import CommandTool, PythonTool, TInput
from janis_core import settings


from ..scope import Scope
from ..process.inputs import data_sources
from .. import params
from .. import channels


def register_data_sources(wf: Workflow) -> None:
    """
    register param(s) for each workflow input. 
    channel(s) may also be registered if necessary.
    """
    scope = Scope()
    do_register_data_sources(wf, scope)

def do_register_data_sources(wf: Workflow, scope: Scope):
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        
        if isinstance(step.tool, CommandTool) or isinstance(step.tool, PythonTool):
            generator = ProcessDataSourceGenerator(step.tool, step.sources)
            generator.generate()
            data_sources.update(
                current_scope, 
                generator.process_inputs, 
                generator.param_inputs, 
                generator.internal_inputs
            )
        elif isinstance(step.tool, Workflow):
            do_register_data_sources(step.tool, current_scope)
        else:
            raise NotImplementedError

        


class ProcessDataSourceGenerator:
    """
    for a CommandTool / PythonTool, decides which of the ToolInputs / TInputs which will be:
        - process inputs
        - param inputs
        - internal inputs
    """
    def __init__(self, tool: CommandTool | PythonTool, sources: dict[str, Any]) -> None:
        self.tool = tool
        self.sources = sources
        self.process_inputs: set[str] = set()
        self.param_inputs: set[str] = set()
        self.internal_inputs: set[str] = set()

    @property
    def tinputs(self) -> list[TInput]:
        if isinstance(self.tool, CommandTool):
            return self.tool.tool_inputs()
        elif isinstance(self.tool, PythonTool):  # type: ignore
            return self.tool.inputs()
        else:
            raise NotImplementedError
    
    @property
    def all_input_ids(self) -> set[str]:
        return {x.id() for x in self.tinputs} 
    
    @property
    def channel_input_ids(self) -> set[str]:
        """get tool inputs (ids) which are linked to a channel"""
        out: set[str] = set()
        for tag, src in self.sources.items():
            node = utils.resolve_node(src)
            if isinstance(node, InputNode):
                if channels.exists(node.uuid):
                    out.add(tag)
        return out

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

    @property
    def connection_input_ids(self) -> set[str]:
        """get tool inputs (ids) which are being fed a value from a step connection"""
        out: set[str] = set()
        for tag, src in self.sources.items():
            node = utils.resolve_node(src)
            if isinstance(node, StepNode):
                out.add(tag)
        return out

    @property
    def scatter_input_ids(self) -> set[str]:
        """get tool inputs (ids) which are being scattered on"""
        out: set[str] = set()
        for inname, src in self.sources.items():
            should_scatter = src.source_map[0].should_scatter
            node = utils.resolve_node(src)
            if should_scatter and isinstance(node, InputNode):
                out.add(inname)
        return out

    @property
    def complex_expression_input_ids(self) -> set[str]:
        """get tool inputs (ids) which are fed data using a complex janis expression"""
        out: set[str] = set()
        for inname, src in self.sources.items():
            node = utils.resolve_node(src)
            if not isinstance(node, InputNode) and not isinstance(node, StepNode):
                out.add(inname)
        return out
    
    ### main method
    def generate(self) -> None:
        self.process_inputs = self.get_process_inputs()
        self.param_inputs = self.get_param_inputs()
        self.internal_inputs = self.get_internal_inputs()

    ### process inputs
    def get_process_inputs(self) -> set[str]:
        """gets the tool inputs which will become nextflow process inputs"""
        if settings.translate.nextflow.MODE == 'workflow':
            return self.get_process_inputs_workflowmode()
        elif settings.translate.nextflow.MODE == 'tool':  # type: ignore
            return self.get_process_inputs_toolmode()
        else:
            raise RuntimeError

    def get_process_inputs_workflowmode(self) -> set[str]:
        """
        inputs which are fed (via step inputs) using a file type workflow input
        inputs which are fed (via step inputs) using a connection
        inputs which are involved in scatter
        """
        return self.channel_input_ids | self.connection_input_ids | self.scatter_input_ids | self.complex_expression_input_ids

    def get_process_inputs_toolmode(self) -> set[str]:
        """all CommandTool inputs. in toolmode, we have no greater scope than the tool itself."""
        return self.all_input_ids
    

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

