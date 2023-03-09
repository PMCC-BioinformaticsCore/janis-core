


from typing import Any, Optional

from janis_core.workflow.workflow import Workflow, OutputNode, InputNodeSelector, StepNode
from janis_core import ScatterDescription, ScatterMethod, StringFormatter

from .common import EntityParser
from ..types import ingest_cwl_type
from ..types import cast_cwl_type_to_python
from ..identifiers import get_id_entity
from ..identifiers import remove_output_name_from_output_source
from ..graph import get_janis_wf_sources
from ..expressions import parse_basic_expression






class WorkflowInputParser(EntityParser):

    def parse(self, entity: Any, wf: Workflow) -> InputNodeSelector:  # type: ignore
        identifier = get_id_entity(entity.id)
        # datatype
        dtype, error_msgs = ingest_cwl_type(entity.type, self.cwl_utils, secondary_files=entity.secondaryFiles)
        self.error_msgs += error_msgs
        
        inp = wf.input(
            identifier=identifier,
            datatype=dtype,
            default=entity.default,
            doc=entity.doc,
        )

        self.log_messages(inp.input_node.uuid)
        return inp


class WorkflowOutputParser(EntityParser):

    def parse(self, entity: Any, wf: Workflow) -> OutputNode:  # type: ignore
        out_identifier = get_id_entity(entity.id)
        if isinstance(entity.outputSource, list):
            cwl_source = [remove_output_name_from_output_source(x) for x in entity.outputSource]
        else:
            cwl_source = remove_output_name_from_output_source(entity.outputSource)
        
        sources = get_janis_wf_sources(wf, cwl_source)
        if len(sources) == 1:
            source = sources[0]
        else:
            source = sources

        dtype, error_msgs = ingest_cwl_type(entity.type, self.cwl_utils, secondary_files=entity.secondaryFiles)
        self.error_msgs += error_msgs

        out = wf.output(
            identifier=out_identifier,
            datatype=dtype,
            source=source,
        )
        
        self.log_messages(out.uuid)
        return out
    

class WorkflowStepInputsParser(EntityParser):

    def parse(self, entity: Any, wf: Workflow) -> dict[str, Any]:  # type: ignore
        step_identifier = get_id_entity(entity.id)
        jstep = wf.step_nodes[step_identifier]

        valid_step_inputs = self.get_valid_step_inputs(entity, jstep)

        inputs_dict = {}
        for inp in valid_step_inputs:
            inp_identifier = get_id_entity(inp.id)
            source = self.get_input_source(wf, inp)
            inputs_dict[inp_identifier] = source

        self.log_messages(jstep.uuid)
        return inputs_dict

    def get_valid_step_inputs(self, cwlstp: Any, jstep: StepNode) -> list[Any]:
        return [x for x in cwlstp.in_ if self.is_valid_step_input(x, jstep)]

    def is_valid_step_input(self, inp: Any, jstep: StepNode) -> bool:
        inp_identifier = get_id_entity(inp.id)
        if inp_identifier in jstep.tool.inputs_map():
            return True
        else:
            print(f'WARNING: [STEP: {jstep.id()}] [TOOL: {jstep.tool.id()}] - Tool has no input named "{inp_identifier}". Ignoring.')
            return False
    
    def get_input_source(self, wf: Workflow, inp: Any) -> Any:
        source = self.resolve_source(wf, inp)
        value_from = self.resolve_value_from(inp, source)
        
        if source is not None and value_from is not None:
            value = value_from
            # src_name = str(source)
            
            # if isinstance(value_from, str):
            #     if 'self' not in value_from:
            #         raise NotImplementedError
            #     value = value_from.replace('self', src_name)
            
            # elif isinstance(value_from, StringFormatter):
            #     for key, val in value_from.kwargs.items():
            #         if 'self' not in val:
            #             raise NotImplementedError
            #         value_from.kwargs[key] = val.replace('self', src_name)
            #     value = value_from

            # else:
            #     raise NotImplementedError
        
        elif source is not None:
            value = source

        elif value_from is not None:
            value = value_from

        elif inp.default is not None:
            value = cast_cwl_type_to_python(inp.default)

        else:
            value = None
            print(f"Source is None from object: {inp.save()}")
        
        return value

    def resolve_source(self, wf: Workflow, inp: Any) -> Any:
        source = None
        if inp.source is not None:
            sources = get_janis_wf_sources(wf, inp.source)
            if len(sources) == 1:
                source = sources[0]
            else:
                source = sources
        return source

    def resolve_value_from(self, inp: Any, source: Any) -> Any:
        value = None

        if inp.valueFrom is not None:
            if 'self.' in inp.valueFrom:
                inp.valueFrom = inp.valueFrom.replace('self.', f'{source}.')

            value, success = parse_basic_expression(inp.valueFrom)
            if not success:
                inp_identifier = get_id_entity(inp.id)
                msg = f"'{inp_identifier}' input value contains untranslated javascript expression"
                self.error_msgs.append(msg)

        return value
        


class WorkflowStepScatterParser(EntityParser):

    def parse(self, entity: Any, jstep: StepNode) -> Optional[ScatterDescription]:  # type: ignore
        
        scatter = None
        if entity.scatter:
            if isinstance(entity.scatter, list):
                scatter_fields_raw: list[str] = entity.scatter
            else:
                scatter_fields_raw: list[str] = [entity.scatter]

            scatter_fields: list[str] = []
            for field in scatter_fields_raw:
                input_to_scatter = field.rsplit("/", 1)[1]
                scatter_fields.append(input_to_scatter)

            scatter_method = self.ingest_scatter_method(entity.scatterMethod)
            scatter = ScatterDescription(fields=scatter_fields, method=scatter_method)
        
        self.log_messages(jstep.uuid)
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