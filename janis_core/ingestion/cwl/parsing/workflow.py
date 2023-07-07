


from typing import Any, Optional
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

from janis_core.workflow.workflow import Workflow, OutputNode, InputNodeSelector, StepNode, InputNode, StepOutputSelector
from janis_core import ScatterDescription, ScatterMethod

from ..types import ingest_cwl_type
from ..identifiers import get_id_entity
from ..identifiers import remove_output_name_from_output_source
from ..graph import get_janis_wf_sources

from janis_core import settings
from janis_core.types import File
from janis_core.messages import log_warning

parsing_error_count = 0



@dataclass
class WorkflowEntityParser(ABC):
    cwl_utils: Any
    entity: Any
    wf: Workflow
    success: bool = False
    uuid: Optional[str] = None
    error_msgs: list[str] = field(default_factory=list)
    failure_message: str = ''
    should_log_messages: bool = True

    def parse(self) -> Any:
        # normal mode
        if settings.ingest.SAFE_MODE:
            try:
                j_entity = self.do_parse()
                self.success = True
            except Exception:
                j_entity = self.fallback()
        
        # dev mode
        else:
            j_entity = self.do_parse()
            self.success = True
        
        self.log_messages()
        return j_entity

    def log_messages(self) -> None:
        if self.success == False:
            self.error_msgs.append(self.failure_message)
        
        if self.should_log_messages:
            for msg in self.error_msgs:
                log_warning(self.uuid, msg)

    @abstractmethod
    def do_parse(self) -> Any:
        ...
    
    @abstractmethod
    def fallback(self) -> Any:
        ...



@dataclass
class WorkflowInputParser(WorkflowEntityParser):
    """parses a cwl WorkflowInputParameter to add an input to the janis Workflow"""

    failure_message: str = 'error parsing workflow input. returned generic optional File input as fallback'

    def fallback(self) -> InputNodeSelector:
        global parsing_error_count
        parsing_error_count += 1 
        
        identifier = f'unknown_{parsing_error_count}'
        inp = self.wf.input(
            identifier=identifier,
            datatype=File,
        )
        # this must be set for error messages to be associated with this entity
        self.uuid = inp.input_node.uuid
        return inp

    def do_parse(self) -> InputNodeSelector:

        identifier = get_id_entity(self.entity.id)
        dtype, error_msgs = ingest_cwl_type(self.entity.type, self.cwl_utils, secondaries=self.entity.secondaryFiles)
        self.error_msgs += error_msgs
        
        inp = self.wf.input(
            identifier=identifier,
            datatype=dtype,
            default=self.entity.default,
            doc=self.entity.doc,
        )
        # this must be set for error messages to be associated with this workflow input
        self.uuid = inp.input_node.uuid
        return inp



@dataclass
class WorkflowOutputParser(WorkflowEntityParser):
    """parses a cwl WorkflowInputParameter to add an input to the janis Workflow"""

    failure_message: str = "error parsing workflow output. \
                            created a generic File workflow output as fallback. \
                            created a generic File workflow input as the fallback's data source."

    def fallback(self) -> OutputNode:
        global parsing_error_count
        parsing_error_count += 1
        
        identifier = f'unknown_{parsing_error_count}'

        inp = self.wf.input(
            identifier,
            datatype=File,
        )

        out = self.wf.output(
            identifier=identifier,
            datatype=File,
            source=inp
        )
        # this must be set for error messages to be associated with this workflow output
        self.uuid = out.uuid
        return out

    def do_parse(self) -> OutputNode: 
        out_identifier = get_id_entity(self.entity.id)
        
        if isinstance(self.entity.outputSource, list):
            cwl_source = [remove_output_name_from_output_source(x) for x in self.entity.outputSource]
        else:
            cwl_source = remove_output_name_from_output_source(self.entity.outputSource)
        
        janis_sources = get_janis_wf_sources(self.wf, cwl_source)
        if len(janis_sources) == 1:
            source = janis_sources[0]
        else:
            source = janis_sources

        dtype, error_msgs = ingest_cwl_type(self.entity.type, self.cwl_utils, secondaries=self.entity.secondaryFiles)
        self.error_msgs += error_msgs

        out = self.wf.output(
            identifier=out_identifier,
            datatype=dtype,
            source=source,
        )
        # this must be set for error messages to be associated with this workflow output
        self.uuid = out.uuid
        return out
    


@dataclass
class WorkflowStepInputsParser(WorkflowEntityParser):

    failure_message: str = "error parsing task invocation inputs. \
                            task invocation will have no inputs as fallback. "

    def fallback(self) -> dict[str, Any]:
        step_identifier = get_id_entity(self.entity.id)
        jstep = self.wf.step_nodes[step_identifier]
        # this must be set for error messages to be associated with this workflow output
        self.uuid = jstep.uuid
        return {}

    def do_parse(self) -> dict[str, Any]: 
        step_identifier = get_id_entity(self.entity.id)
        jstep = self.wf.step_nodes[step_identifier]

        # this must be set for error messages to be associated with this workflow output
        self.uuid = jstep.uuid
        
        valid_step_inputs = self.get_valid_step_inputs(self.entity, jstep)

        inputs_dict = {}
        for inp in valid_step_inputs:
            inp_identifier = get_id_entity(inp.id)
            parser = WorkflowStepInputParser(cwl_utils=self.cwl_utils, entity=inp, wf=self.wf)
            source = parser.parse()
            if source is not None:
                inputs_dict[inp_identifier] = source
            else:
                raise Exception("should this be checked?")
                print()
            # collect error messages from step inputs, log them as errors parsing this step instead
            self.error_msgs += parser.error_msgs  

        return inputs_dict

    def get_valid_step_inputs(self, cwlstp: Any, jstep: StepNode) -> list[Any]:
        return [x for x in cwlstp.in_ if self.is_valid_step_input(x, jstep)]

    def is_valid_step_input(self, inp: Any, jstep: StepNode) -> bool:
        inp_identifier = get_id_entity(inp.id)
        if inp_identifier in jstep.tool.inputs_map():
            return True
        else:
            msg = f'{jstep.tool.id()} task has no input named "{inp_identifier}". Ignored as fallback.'
            self.error_msgs.append(msg)
            return False
        
    

@dataclass
class WorkflowStepInputParser(WorkflowEntityParser):

    failure_message: str = ''
    should_log_messages: bool = False
    
    def fallback(self) -> Optional[InputNode | StepOutputSelector]:
        identifier = get_id_entity(self.entity.id)
        self.failure_message = f"error parsing data source for '{identifier}'. Ignored as fallback. "
        return None

    def do_parse(self) -> Optional[InputNode | StepOutputSelector]: 
        identifier = get_id_entity(self.entity.id)
        self.failure_message = f"error parsing data source for '{identifier}'. Ignored as fallback. "
        
        source = self.resolve_source()
        value_from = self.resolve_value_from(source)
        
        if source is not None and value_from is not None:
            value = source
        
        elif source is not None:
            value = source

        elif value_from is not None:
            value = value_from

        elif self.entity.default is not None:
            value = self.entity.default

        else:
            value = None
            # print(f"Source is None from object: {entity.save()}")
        
        return value

    def resolve_source(self) -> Any:
        inp = self.entity
        source = None
        if inp.source is not None:
            sources = get_janis_wf_sources(self.wf, inp.source)
            if len(sources) == 1:
                source = sources[0]
            else:
                source = sources
        return source

    def resolve_value_from(self, source: Any) -> Any:
        inp = self.entity
        value = None

        if inp.valueFrom is not None:
            if 'self.' in inp.valueFrom:
                if isinstance(source, InputNodeSelector):
                    replacement = f'inputs.{source.input_node.id()}.'
                    inp.valueFrom = inp.valueFrom.replace('self.', replacement)
                elif isinstance(source, StepOutputSelector):
                    replacement = f'steps.{source.node.id()}.{source.tag}.'
                    inp.valueFrom = inp.valueFrom.replace('self.', replacement)
                else:
                    raise NotImplementedError

            # we do not want to unwrap inp.valueFrom. doesnt translate to nextflow correctly.
            # just mark it as untranslated javascript after the 'self.' has been changed to something
            # more meaningful
            value = f'<js>{inp.valueFrom}</js>'
            
            inp_identifier = get_id_entity(inp.id)
            msg = f"'{inp_identifier}' input value contains untranslated javascript expression: {value}"
            self.error_msgs.append(msg)

        return value
    
    

@dataclass
class WorkflowStepScatterParser(WorkflowEntityParser):
    """parses janis step to return a ScatterDescription."""

    failure_message: str = 'error parsing task parallelisation. ignored parallelisation as fallback'
    
    def fallback(self) -> Optional[ScatterDescription]:
        return None

    def do_parse(self) -> Optional[ScatterDescription]:
        scatter = None

        if self.entity.scatter:
            if isinstance(self.entity.scatter, list):
                scatter_fields_raw: list[str] = self.entity.scatter
            else:
                scatter_fields_raw: list[str] = [self.entity.scatter]

            scatter_fields: list[str] = []
            for field in scatter_fields_raw:
                input_to_scatter = field.rsplit("/", 1)[1]
                scatter_fields.append(input_to_scatter)

            scatter_method = self.ingest_scatter_method(self.entity.scatterMethod)
            scatter = ScatterDescription(fields=scatter_fields, method=scatter_method)
        
        return scatter
        
    def ingest_scatter_method(self, scatter_method: Optional[str]) -> Optional[ScatterMethod]:
        if scatter_method is None or scatter_method == "":
            return None
        
        elif scatter_method == "dotproduct":
            return ScatterMethod.dot
        
        elif scatter_method == "nested_crossproduct":
            msg = 'task parallelisation method was nested crossproduct, but Janis only supports flat crossproduct. Used flat crossproduct as fallback.'
            self.error_msgs.append(msg)
            # j.Logger.warn(
            #     "Requesting nested_crossproduct, but Janis only supports flat_crossproduct. Will fallback to flat_crossproduct"
            # )
            return ScatterMethod.cross
        
        elif scatter_method == "flat_crossproduct":
            return ScatterMethod.cross
        
        else:
            msg = f"Unsupported scatter method '{scatter_method}'. Used flat crossproduct as fallback."
            self.error_msgs.append(msg)