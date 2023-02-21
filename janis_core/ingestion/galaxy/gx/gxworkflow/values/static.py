
from __future__ import annotations
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    from janis_core.ingestion.galaxy.model.workflow import WorkflowStep
    from janis_core.ingestion.galaxy.model.workflow import Workflow
    from janis_core.ingestion.galaxy.model.workflow import WorkflowInput

from janis_core.ingestion.galaxy.logs import logging

from janis_core.ingestion.galaxy.gx.command.cmdstr import gen_command_string
from janis_core.ingestion.galaxy.gx.command.components import Flag
from janis_core.ingestion.galaxy.gx.command.components import InputComponent
from janis_core.ingestion.galaxy.gx.command.components import Option

from janis_core.ingestion.galaxy.gx.gxtool import load_xmltool
from janis_core.ingestion.galaxy.gx.gxtool.text import load_partial_cheetah_command

from . import factory
from janis_core.ingestion.galaxy.model.workflow import InputValue

from janis_core.ingestion.galaxy import mapping
from janis_core.ingestion.galaxy import datatypes
from janis_core.ingestion.galaxy import settings
from janis_core.ingestion.galaxy import expressions


def handle_step_static_inputs(janis: Workflow, galaxy: dict[str, Any]) -> None:
    """supplied values in step 'tool_state'"""
    for g_step in galaxy['steps'].values():
        if g_step['type'] == 'tool':
            j_step = mapping.step(g_step['id'], janis, galaxy)
            settings.tool.set(from_wrapper=j_step.metadata.wrapper)
            ingest_values_cheetah(g_step, j_step, janis)
            ingest_values_inputs(g_step, j_step, janis)

def ingest_values_cheetah(g_step: dict[str, Any], j_step: WorkflowStep, janis: Workflow) -> None:
    ingestor = CheetahInputIngestor(g_step, j_step, janis)
    ingestor.ingest()

def ingest_values_inputs(g_step: dict[str, Any], j_step: WorkflowStep, janis: Workflow) -> None:
    ingestor = StaticInputIngestor(g_step, j_step, janis)
    ingestor.ingest()


class CheetahInputIngestor:
    def __init__(self, g_step: dict[str, Any], j_step: WorkflowStep, janis: Workflow):
        self.g_step = g_step
        self.j_step = j_step
        self.janis = janis

    def ingest(self) -> None:
        cmdstr = self.prepare_command()
        for component in self.get_linkable_components():
            match component:
                case Flag():
                    self.link_flag(component, cmdstr)
                case Option():
                    self.link_option(component, cmdstr)
                case _:
                    pass

    def prepare_command(self) -> str:
        xmltool = load_xmltool(settings.tool.tool_path)
        command = load_partial_cheetah_command(inputs_dict=self.g_step['tool_state'])
        cmdstr = gen_command_string(source='xml', the_string=command, xmltool=xmltool)
        stmtstr = cmdstr.main.cmdline
        # logging.runtime_data(command)
        # logging.runtime_data(stmtstr)
        return stmtstr
    
    def link_flag(self, flag: Flag, cmdstr: str) -> None:
        """
        links a flag component value as None if not in cmdstr
        should only detect the flag's absense, nothing else
        """
        if not expressions.is_present(flag.prefix, cmdstr):
            self.handle_not_present_flag(flag)

    def link_option(self, option: Option, cmdstr: str) -> None:
        """gets the value for a specific tool argument"""
        value = expressions.get_next_word(option.prefix, option.delim, cmdstr)
        value = None if value == '' else value
        if value is None:
            self.handle_not_present_opt(option)
        elif self.is_param(value):
            self.handle_gxvar_opt(option, value)
        elif expressions.is_var(value) or expressions.has_var(value):
            self.handle_envvar_opt(option, value)
        else:
            self.handle_value_opt(option, value)

    # TODO upgrade for pre/post task section
    def is_param(self, text: Optional[str]) -> bool:
        # how does this even work? 
        # $ can also mean env var? 
        if text:
            if text[0] == '$':
                return True
            elif len(text) > 1 and text[1] == '$':  # WTF 
                return True
        return False

    def handle_not_present_flag(self, flag: Flag) -> None:
        self.update_tool_values_static(component=flag, value=False)

    def handle_not_present_opt(self, option: Option) -> None:
        self.update_tool_values_static(component=option, value=None)
    
    def handle_gxvar_opt(self, option: Option, value: str) -> None:
        # TODO future: attach the identified param if not attached? 
        # should always be attached tho? 
        pass
    
    def handle_envvar_opt(self, option: Option, value: Any) -> None:
        self.update_tool_values_runtime(component=option)
    
    def handle_value_opt(self, option: Option, value: Any) -> None:
        self.update_tool_values_static(component=option, value=value)

    def update_tool_values_static(self, component: Flag | Option, value: Any) -> None:
        # create & add value 
        is_default = True if component.default_value == value else False
        inputval = factory.static(component, value, default=is_default)
        self.j_step.inputs.add(inputval)
    
    def update_tool_values_runtime(self, component: Flag | Option) -> None:
        # create & add new workflow input
        winp = self.create_workflow_input(component)
        self.janis.add_input(winp)
        # create & add value 
        inputval = factory.workflow_input(component, winp.uuid, is_runtime=True)
        self.j_step.inputs.add(inputval)

    def create_workflow_input(self, component: Flag | Option) -> WorkflowInput:
        """creates a workflow input for the tool input component"""
        return WorkflowInput(
            _name=component.tag,
            array=component.array,
            is_runtime=True,
            datatype=datatypes.get(component),
            optional=component.optional
        )

    def get_linkable_components(self) -> list[InputComponent]:
        out: list[InputComponent] = []
        for component in self.j_step.tool.inputs:
            if not self.j_step.inputs.get(component.uuid):
                out.append(component)
        return out



class StaticInputIngestor:
    def __init__(self, g_step: dict[str, Any], j_step: WorkflowStep, janis: Workflow):
        self.g_step = g_step
        self.j_step = j_step
        self.janis = janis

    def ingest(self) -> None:
        for component in self.get_linkable_components():
            if self.is_directly_linkable(component):
                value = self.create_value(component)
                self.j_step.inputs.add(value)
    
    def get_linkable_components(self) -> list[InputComponent]:
        out: list[InputComponent] = []
        # tool components which don't yet appear in register
        tool_inputs = self.j_step.tool.inputs
        tool_values = self.j_step.inputs
        for component in tool_inputs:
            if component.gxparam and not tool_values.get(component.uuid):
                out.append(component)
        return out

    def is_directly_linkable(self, component: InputComponent) -> bool:
        """
        checks whether a janis tool input can actually be linked to a value in the 
        galaxy workflow step.
        only possible if the component has a gxparam, and that gxparam is referenced as a
        ConnectionStepInput, RuntimeStepInput or StaticStepInput
        """
        if component.gxparam:
            query = component.gxparam.name 
            if query in self.g_step['tool_state']:
                return True
        return False

    def create_value(self, component: InputComponent) -> InputValue:
        """
        create an InputValue for this tool component using supplied step 'tool_state' input values.
        we know that the component has a galaxy param, and that the same galaxy param
        has a supplied value (or connection or runtime value specified) in the step input values. 
        this function grabs that supplied value, then creates a formalised InputValue. 
        it also does other necessary functions - in the case of a galaxy 'runtime value', 
        this needs to become a WorkflowInput in janis world. a WorkflowInput would be created, 
        added to the Workflow, and also added to the InputValue (WorkflowInputInputValue). 
        this says 'for step x using tool y, the tool input component z gets its value from 
        our new WorkflowInput'
        """
        # pull value from 'tool_state'
        # should only be static values left
        # this should be really ez?
        
        g_value = self.g_step['tool_state'][component.gxparam.name] # type: ignore
        is_default = True if component.default_value == g_value else False
        return factory.static(component, value=g_value, default=is_default)
