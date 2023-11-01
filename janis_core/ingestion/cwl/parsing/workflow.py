


from typing import Any, Optional
from dataclasses import dataclass
from abc import ABC, abstractmethod

from janis_core.workflow.workflow import Workflow, OutputNode, InputNodeSelector, StepNode, InputNode, StepOutputSelector
from janis_core.operators import FirstOperator, FilterNullOperator
from janis_core import ScatterDescription, ScatterMethod

from ..types import ingest_cwl_type
from ..identifiers import get_id_entity
from ..identifiers import remove_output_name_from_output_source
from ..graph import get_janis_wf_sources
from ..expressions import parse_expression

from janis_core import settings
from janis_core.types import File
from janis_core.messages import log_error
from janis_core.messages import ErrorCategory

parsing_error_count = 0


@dataclass
class WorkflowEntityParser(ABC):
    cwl_utils: Any
    entity: Any
    wf: Workflow
    tool_uuid: str
    success: bool = False

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
        
        return j_entity

    @abstractmethod
    def do_parse(self) -> Any:
        ...
    
    @abstractmethod
    def fallback(self) -> Any:
        ...



@dataclass
class WorkflowInputParser(WorkflowEntityParser):
    """parses a cwl WorkflowInputParameter to add an input to the janis Workflow"""

    def fallback(self) -> InputNodeSelector:
        # log message
        msg = f'error parsing {get_id_entity(self.entity.id)}. returned generic optional File input as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK, subsection='inputs')

        # fallback
        global parsing_error_count
        parsing_error_count += 1 
        
        identifier = f'unknown_{parsing_error_count}'
        inp = self.wf.input(
            identifier=identifier,
            datatype=File(optional=True),
        )
        return inp

    def do_parse(self) -> InputNodeSelector:
        identifier = get_id_entity(self.entity.id)
        dtype = ingest_cwl_type(self.entity.type, self.cwl_utils, self.entity, self.tool_uuid, secondaries=self.entity.secondaryFiles)
        
        inp = self.wf.input(
            identifier=identifier,
            datatype=dtype,
            default=self.entity.default,
            doc=self.entity.doc,
        )
        return inp



@dataclass
class WorkflowOutputParser(WorkflowEntityParser):
    """parses a cwl WorkflowInputParameter to add an input to the janis Workflow"""

    def fallback(self) -> OutputNode:
        # log message
        msg = f'error parsing {get_id_entity(self.entity.id)}. returned generic File output (and corresponding input) as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK, subsection='outputs')
        
        # fallback
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

        dtype = ingest_cwl_type(self.entity.type, self.cwl_utils, self.entity, self.tool_uuid, secondaries=self.entity.secondaryFiles)

        out = self.wf.output(
            identifier=out_identifier,
            datatype=dtype,
            source=source,
        )
        self.ingest_workflow_output_pickvalue(out)
        return out

    def ingest_workflow_output_pickvalue(self, out: OutputNode) -> None:
        # new
        operator_map = {
            'first_non_null': FirstOperator,
            'the_only_non_null': FirstOperator,
            'all_non_null': FilterNullOperator,
        }
            
        # must have more than one source
        if not isinstance(out.source, list) or len(out.source) <= 1:
            return 

        # must have pickValue set on cwl step input
        if not hasattr(self.entity, 'pickValue') or self.entity.pickValue is None:
            return 

        operator = operator_map[self.entity.pickValue]
        out.operator = operator(out.source)
    


@dataclass
class WorkflowStepInputsParser(WorkflowEntityParser):

    def fallback(self) -> dict[str, Any]:
        # log message
        step_identifier = get_id_entity(self.entity.id)
        jstep = self.wf.step_nodes[step_identifier]
        msg = f"error parsing step inputs"
        log_error(self.tool_uuid, msg, ErrorCategory.PLUMBING, subsection=f'step:{jstep.id()}')
        
        # fallback
        return {}

    def do_parse(self) -> dict[str, Any]: 
        step_identifier = get_id_entity(self.entity.id)
        jstep = self.wf.step_nodes[step_identifier]

        # valid_step_inputs = self.get_valid_step_inputs(self.entity, jstep)

        inputs_dict = {}
        for inp in self.entity.in_:
            inp_identifier = get_id_entity(inp.id)
            parser = WorkflowStepInputParser(cwl_utils=self.cwl_utils, entity=inp, wf=self.wf, tool_uuid=self.wf.uuid)
            parser.step_name = jstep.id()
            source = parser.parse()
            inputs_dict[inp_identifier] = source

        return inputs_dict

    # this may have been a mistake to comment out
    # def get_valid_step_inputs(self, cwlstp: Any, jstep: StepNode) -> list[Any]:
    #     return [x for x in cwlstp.in_ if self.is_valid_step_input(x, jstep)]

    # def is_valid_step_input(self, inp: Any, jstep: StepNode) -> bool:
    #     inp_identifier = get_id_entity(inp.id)
    #     if inp_identifier in jstep.tool.inputs_map():
    #         return True
    #     else:
    #         raise RuntimeError
    #         # msg = f'{jstep.tool.id()} task has no input named "{inp_identifier}". Ignored as fallback.'
    #         # self.error_msgs.append(msg)
    #         # return False
        
    

@dataclass
class WorkflowStepInputParser(WorkflowEntityParser):
    step_name: str = ''

    def fallback(self) -> Optional[InputNode | StepOutputSelector]:
        # log message
        identifier = get_id_entity(self.entity.id)
        msg = f"error parsing data source for '{identifier}'. Returned None as fallback."
        log_error(self.tool_uuid, msg, ErrorCategory.PLUMBING, subsection=f'step:{self.step_name}')
        
        # fallback
        return None

    def do_parse(self) -> Optional[InputNode | StepOutputSelector]: 
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
        # TODO HERE 
        # raise NotImplementedError
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
class WorkflowStepAttributesParser(WorkflowEntityParser):
    """parses janis step to return a ScatterDescription."""

    scatter: Optional[ScatterDescription] = None
    when: Optional[Any] = None

    def fallback(self) -> None:
        # log message
        step_identifier = get_id_entity(self.entity.id)
        jstep = self.wf.step_nodes[step_identifier]
        msg = 'error parsing step modifiers (scatter | conditional execution). ignored as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.PLUMBING, subsection=f'step:{jstep.id()}')
        
        # fallback
        return None

    def do_parse(self) -> None:
        self.scatter = self.parse_scatter()
        self.when = self.parse_when()

    def parse_when(self) -> Any:
        if hasattr(self.entity, 'when') and self.entity.when is not None:
            res, success = parse_expression(self.entity.when, self.tool_uuid)
            return res
        return None

    def parse_scatter(self) -> Optional[ScatterDescription]:
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
        step_identifier = get_id_entity(self.entity.id)
        jstep = self.wf.step_nodes[step_identifier]

        if scatter_method is None or scatter_method == "":
            return None
        
        elif scatter_method == "dotproduct":
            return ScatterMethod.dot
        
        elif scatter_method == "nested_crossproduct":
            msg = 'task parallelisation method was nested crossproduct, but Janis only supports flat crossproduct. Used flat crossproduct as fallback.'
            log_error(self.tool_uuid, msg, ErrorCategory.PLUMBING, subsection=f'step:{jstep.id()}')
            return ScatterMethod.cross
        
        elif scatter_method == "flat_crossproduct":
            return ScatterMethod.cross
        
        else:
            msg = f"Unsupported scatter method '{scatter_method}'. Used flat crossproduct as fallback."
            log_error(self.tool_uuid, msg, ErrorCategory.PLUMBING, subsection=f'step:{jstep.id()}')