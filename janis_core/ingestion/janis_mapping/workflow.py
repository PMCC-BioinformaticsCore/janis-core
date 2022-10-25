

from typing import Any, Optional
from datetime import datetime

from galaxy2janis.entities.workflow import Workflow
from galaxy2janis.entities.workflow import WorkflowMetadata
from galaxy2janis.entities.workflow import InputValue
from galaxy2janis.entities.workflow import WorkflowStep
from galaxy2janis.entities.workflow import StaticInputValue
from galaxy2janis.entities.workflow import ConnectionInputValue
from galaxy2janis.entities.workflow import WorkflowInputInputValue
from galaxy2janis.runtime.dates import JANIS_DATE_FMT
from galaxy2janis import tags

from janis_core import WorkflowBuilder
from janis_core import CommandToolBuilder
from janis_core import ScatterDescription
from janis_core import ScatterMethods
from janis_core import WorkflowMetadata as JanisWorkflowMetadata
from janis_core import InputNodeSelector
from janis_core import StepOutputSelector
from janis_core import ToolInput
from janis_core.types import File
from janis_core.types import Int
from janis_core.types import Float
from janis_core.types import String
from janis_core.types import Boolean
from janis_core.types import ParseableType

from .general import to_janis_datatype
from .tool import to_janis_tool


### MODULE EXPORTS

def to_janis_workflow(internal: Workflow) -> WorkflowBuilder:
    """
    maps internal model workflow to janis model workflow
    missing the following (unnecessary?):
        friendly_name
        tool_provider
        tool_module
    """
    return Mapper(internal).map()

def to_janis_inputs_dict(internal: Workflow) -> dict[str, Any]:
    out_dict: dict[str, Any] = {}
    for internal_inp in internal.inputs:
        key = internal_inp.tag
        val = internal_inp.value if internal_inp.value else None
        out_dict[key] = val
    return out_dict


class Mapper:
    def __init__(self, gx_model: Workflow):
        self.gx_model = gx_model
        self.janis_model = WorkflowBuilder(
            identifier=gx_model.tag,
            version=gx_model.metadata.version,
            metadata=self.map_metadata(gx_model.metadata),
            doc=gx_model.metadata.annotation
        )
    
    def map(self) -> WorkflowBuilder:
        self.add_inputs()
        self.add_steps()
        self.add_outputs()
        return self.janis_model

    def map_metadata(self, meta: WorkflowMetadata) -> JanisWorkflowMetadata:
        contributors = ['gxtool2janis']
        return JanisWorkflowMetadata(
            short_documentation=meta.annotation,
            contributors=contributors,
            keywords=meta.tags,
            dateCreated=datetime.today().strftime(JANIS_DATE_FMT),
            dateUpdated=datetime.today().strftime(JANIS_DATE_FMT),
            version=meta.version
        )
    
    def add_inputs(self) -> None:
        """missing wf.input.default (not supported in galaxy)"""
        for inp in self.gx_model.inputs:
            self.janis_model.input(
                identifier=inp.tag,
                datatype=to_janis_datatype(inp),
                value=inp.value,
                doc=inp.docstring, # type: ignore
            )
    
    def add_steps(self) -> None:
        """
        missing the following (unnecessary / feature not in galaxy):
            _foreach
            when
            ignore_missing
        """
        for step in self.gx_model.steps:
            step_tag = step.tag
            scatter_names = _get_scatter_targets(step)
            scatter_obj = _get_scatter_object(scatter_names)
            tool = to_janis_tool(step.tool)
            self.set_tool_input_values(step, tool)
            self.janis_model.step(
                identifier=step_tag,
                tool=tool,
                scatter=scatter_obj,        # type: ignore
                doc=step.metadata.label
            )
    
    def add_outputs(self) -> None:
        """
        missing the following (unnecessary):
            output_folder
            output_name
            extension
        """
        for outp in self.gx_model.outputs:
            step_tag = tags.get(outp.step_uuid)
            step_node = self.janis_model.nodes[step_tag]
            stepout_tag = outp.tool_output.tag
            self.janis_model.output(
                identifier=f'{step_tag}_{stepout_tag}',
                datatype=to_janis_datatype(outp),
                source=(step_node, stepout_tag),
                doc=outp.tool_output.docstring,  # type: ignore
            )

    def set_tool_input_values(self, gstep: WorkflowStep, jtool: CommandToolBuilder) -> None:
        """
        cast galaxy2janis InputValues to the janis model.
        janis stores these step input values on the tool itself (tool.connections).
        """
        unknown_count: int = 0
        invalues = gstep.inputs.all
        invalues = [x for x in invalues if not isinstance(x, StaticInputValue)]
        for gvalue in invalues:
            # tool input is known (normal case)
            if gvalue.component:       
                tag = gvalue.input_tag

            # tool input is unknown
            # create a dummy 'unknown' tool input and add to tool. 
            else:   
                unknown_count += 1
                tag = f'unknown{unknown_count}'
                inp = self.create_unknown_tool_input(tag, gvalue)
                jtool._inputs.append(inp)  # yes, I know this is protected. I need to add an input - sue me.
            
            value = self.get_input_value(gvalue)
            jtool.connections[tag] = value
    
    def create_unknown_tool_input(self, tag: str, gvalue: InputValue) -> ToolInput:
        """
        create a dummy 'unknown' tool input. 
        Used in situations where we know a step is supplied this InputValue, 
        but don't know which underlying tool input it maps to. 
        """
        return ToolInput(
            tag=tag,
            input_type=self.infer_type_from_value(gvalue),
            position=None,
            prefix=None,
            default=gvalue.raw_value if isinstance(gvalue, StaticInputValue) else None,
            doc='ATTENTION: unknown input. please address.'
        )

    def infer_type_from_value(self, gvalue: InputValue) -> ParseableType:
        """guess a type using an InputValue (used when creating 'unknown' tool input)"""
        jtype = File  # fallback 
        if isinstance(gvalue, StaticInputValue):
            jtype = _infer_type_from_static(gvalue)
        elif isinstance(gvalue, WorkflowInputInputValue):
            jtype = self.infer_type_from_input(gvalue)
        elif isinstance(gvalue, ConnectionInputValue):
            jtype = self.infer_type_from_connection(gvalue)
        return jtype
    
    def infer_type_from_input(self, gvalue: WorkflowInputInputValue) -> ParseableType:
        """get the relevant workflow input which supplies the value, then return its type"""
        ginput = self.gx_model.get_input(gvalue.input_uuid)
        return to_janis_datatype(ginput)

    def infer_type_from_connection(self, gvalue: ConnectionInputValue) -> ParseableType:
        """get the relevant step output which supplies the value, then return its type"""
        step = self.gx_model.get_step(gvalue.step_uuid)
        for out in step.outputs.list():
            if out.uuid == gvalue.step_uuid:
                return to_janis_datatype(out)
        return File

    def get_input_value(self, gvalue: InputValue) -> Any | InputNodeSelector | StepOutputSelector:
        """casts a galaxy2janis InputValue to janis model"""
        if isinstance(gvalue, StaticInputValue):
            return gvalue.raw_value
        elif isinstance(gvalue, WorkflowInputInputValue):
            str_value = gvalue.wrapped_value
            node_tag = str_value.split('.')[-1]
            return InputNodeSelector(self.janis_model.input_nodes[node_tag])
        elif isinstance(gvalue, ConnectionInputValue):
            str_value = gvalue.wrapped_value
            node_tag = str_value.split('.')[-2]
            out_tag = str_value.split('.')[-1]
            return StepOutputSelector(self.janis_model.step_nodes[node_tag], out_tag)
        else:
            raise RuntimeError("gvalue must be a galaxy2janis 'InputValue'")



### helper functions

def _infer_type_from_static(gvalue: StaticInputValue) -> ParseableType:
    """
    use the static value to guess a type. 
    If '1', will return Int. if 'hello', will return String etc. 
    """
    type_map = {
        'boolean': Boolean,
        'string': String,
        'float': Float,
        'int': Int,
        'none': File
    }
    return type_map[gvalue.value_type]

def _get_scatter_targets(gstep: WorkflowStep) -> list[str]:
    unknown_count: int = 0
    targets: list[str] = []
    
    invalues = gstep.inputs.all
    invalues = [x for x in invalues if not isinstance(x, StaticInputValue)]
    for v in invalues:
        if v.scatter:
            if v.component:
                name = v.input_tag
            else:
                unknown_count += 1
                name = f'unknown{unknown_count}'
            targets.append(name)
    return targets

def _get_scatter_object(names: list[str]) -> Optional[str | ScatterDescription]:
    if len(names) == 0:
        return None
    elif len(names) == 1:
        return names[0]
    else:
        return ScatterDescription(fields=names, method=ScatterMethods.dot)

