
from __future__ import annotations
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowStep
    from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowInput


from janis_core.ingestion.galaxy.gxtool.parsing import load_xmltool
from janis_core.ingestion.galaxy.gxworkflow import load_tool_state

from janis_core.ingestion.galaxy.gxtool.command.cmdstr import gen_command_string
from janis_core.ingestion.galaxy.gxtool.command.cmdstr.CommandString import CommandStringSource
from janis_core.ingestion.galaxy.gxtool.command.components import Flag
from janis_core.ingestion.galaxy.gxtool.command.components import InputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import Option
from janis_core.ingestion.galaxy.gxtool.command import load_templated_command_str

from . import factory
from janis_core.ingestion.galaxy.internal_model.workflow import InputValue

from janis_core.ingestion.galaxy import internal_mapping
from janis_core.ingestion.galaxy import datatypes
from janis_core.ingestion.galaxy import runtime
from janis_core.ingestion.galaxy import expressions


def handle_step_static_inputs(i_workflow: Workflow, galaxy: dict[str, Any]) -> None:
    """supplied values in step 'tool_state'"""
    for g_step in galaxy['steps'].values():
        if g_step['type'] == 'tool':
            i_step = internal_mapping.step(g_step['id'], i_workflow, galaxy)
            runtime.tool.set(from_wrapper=i_step.metadata.wrapper)
            ingest_values_cheetah(g_step, i_step, i_workflow)
            ingest_values_inputs(g_step, i_step, i_workflow)

def ingest_values_cheetah(g_step: dict[str, Any], i_step: WorkflowStep, i_workflow: Workflow) -> None:
    ingestor = CheetahInputIngestor(g_step, i_step, i_workflow)
    ingestor.ingest()

def ingest_values_inputs(g_step: dict[str, Any], i_step: WorkflowStep, i_workflow: Workflow) -> None:
    ingestor = StaticInputIngestor(g_step, i_step, i_workflow)
    ingestor.ingest()


class CheetahInputIngestor:
    """
    This is the first method of ingesting workflow step input values. 
    CheetahInputIngestor exists because we want to identify which tool flags & options don't 
    appear in the templated <command> section. 
    For those which don't appear in the <command> section, we can set their value to None.
    For those which do appear, we check if they have a galaxy param attached.
    If they have a galaxy param attached, they're handled in StaticInputIngestor.
    If they dont have a galaxy param attached, we set their value to the following token. 
    """

    def __init__(self, g_step: dict[str, Any], i_step: WorkflowStep, i_workflow: Workflow):
        self.g_step = g_step
        self.i_step = i_step
        self.i_workflow = i_workflow

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
        xmltool = load_xmltool(runtime.tool.tool_path)
        tool_state = load_tool_state(
            self.g_step, 
            additional_filters=[
                'ReplaceNullWithVarname'
                'ReplaceConnectedWithVarname',
                'ReplaceRuntimeWithVarname',
            ]
        )
        command = load_templated_command_str(inputs_dict=tool_state)
        cmdstr = gen_command_string(source=CommandStringSource.TOOL_STATE, text=command, xmltool=xmltool)
        stmtstr = cmdstr.main.cmdline
        # logging.runtime_data(command)
        # logging.runtime_data(stmtstr)
        return stmtstr
    
    def get_linkable_components(self) -> list[InputComponent]:
        out: list[InputComponent] = []
        for component in self.i_step.tool.inputs:
            if not self.i_step.inputs.get(component.uuid):
                out.append(component)
        return out
    
    def link_flag(self, flag: Flag, cmdstr: str) -> None:
        """
        links a flag component value as None if not in cmdstr
        should only detect the flag's absense, nothing else
        """
        if not expressions.is_present(flag.prefix, cmdstr):
            self.update_tool_values_static(component=flag, value=False)

    def link_option(self, option: Option, cmdstr: str) -> None:
        """gets the value for a specific tool argument"""
        if option.gxparam is not None:
            return
        value = expressions.get_next_word(option.prefix, option.separator, cmdstr)
        value = None if value == '' else value
        if value is None:
            self.update_tool_values_static(component=option, value=None)
        # elif self.is_param(value):
        #     pass
        elif expressions.is_var(value) or expressions.has_var(value):
            self.update_tool_values_runtime(component=option)
        else:
            self.update_tool_values_static(component=option, value=value)

    # # TODO upgrade for pre/post task section
    # def is_param(self, text: Optional[str]) -> bool:
    #     # how does this even work? 
    #     # $ can also mean env var? 
    #     if text:
    #         if text.strip('"\'')[0] == '$':
    #             return True
    #         elif len(text) > 1 and text[1] == '$':  # WTF 
    #             return True
    #     return False

    # def handle_gxvar_opt(self, option: Option, value: str) -> None:
    #     # TODO future: attach the identified param if not attached? 
    #     # should always be attached tho? 
    #     pass
    
    def update_tool_values_static(self, component: Flag | Option, value: Any) -> None:
        # create & add value 
        is_default = True if component.default_value == value else False
        inputval = factory.static(component, value, default=is_default)
        self.i_step.inputs.add(inputval)
    
    def update_tool_values_runtime(self, component: Flag | Option) -> None:
        # create & add new workflow input
        winp = self.create_workflow_input(component)
        self.i_workflow.add_input(winp)
        # create & add value 
        inputval = factory.workflow_input(component, winp.uuid, is_runtime=True)
        self.i_step.inputs.add(inputval)

    def create_workflow_input(self, component: Flag | Option) -> WorkflowInput:
        """creates a workflow input for the tool input component"""
        return WorkflowInput(
            _name=f'{self.i_step.tag}_{component.tag}',
            array=component.array,
            is_runtime=True,
            datatype=datatypes.get(component),
            optional=component.optional
        )



class StaticInputIngestor:
    """
    This is the fallback method of ingesting workflow step input values. 
    Sometimes the CheetahInputIngestor method results in the <command> section missing 
    some tool arguments, so we cant link them. 
    In this case, we use the full <command> section to link values. 
    get_linkable_components() checks to see which tool input components are not yet linked,
    then for each of those tries to link the tool_state value to a tool input component. 
    If the component already has a value (ie its a runtime value, connected value, or was linked using 
    CheetahInputIngestor), then we skip it.
    """

    def __init__(self, g_step: dict[str, Any], i_step: WorkflowStep, i_workflow: Workflow):
        self.g_step = g_step
        self.i_step = i_step
        self.i_workflow = i_workflow

    def ingest(self) -> None:
        for component in self.get_linkable_components():
            if self.is_directly_linkable(component):
                self.update_tool_values_static(component)
    
    def get_linkable_components(self) -> list[InputComponent]:
        out: list[InputComponent] = []
        # tool components which don't yet appear in register
        tool_inputs = self.i_step.tool.inputs
        step_inputs = self.i_step.inputs
        for component in tool_inputs:
            if component.gxparam and not step_inputs.get(component.uuid):
                out.append(component)
        return out

    def is_directly_linkable(self, component: InputComponent) -> bool:
        """
        checks whether a i_workflow tool input can actually be linked to a value in the 
        galaxy workflow step.
        only possible if the component has a gxparam, and that gxparam is referenced as a
        ConnectionStepInput, RuntimeStepInput or StaticStepInput
        """
        if component.gxparam:
            query = component.gxparam.name 
            tool_state = load_tool_state(self.g_step, additional_filters=['Flatten', 'DeNestClass'])
            if query in tool_state:
                return True
        return False

    def update_tool_values_static(self, component: InputComponent) -> None:
        """
        create an InputValue for this tool component using supplied step 'tool_state' input values.
        we know that the component has a galaxy param, and that the same galaxy param
        has a supplied value (or connection or runtime value specified) in the step input values. 
        this function grabs that supplied value, then creates a formalised InputValue. 
        it also does other necessary functions - in the case of a galaxy 'runtime value', 
        this needs to become a WorkflowInput in i_workflow world. a WorkflowInput would be created, 
        added to the Workflow, and also added to the InputValue (WorkflowInputInputValue). 
        this says 'for step x using tool y, the tool input component z gets its value from 
        our new WorkflowInput'
        """
        # pull value from 'tool_state'
        # should only be static values left
        tool_state = load_tool_state(self.g_step, additional_filters=['Flatten', 'DeNestClass'])
        g_value = tool_state[component.gxparam.name] # type: ignore
        is_default = True if component.default_value == g_value else False
        value = factory.static(component, value=g_value, default=is_default)
        self.i_step.inputs.add(value)
