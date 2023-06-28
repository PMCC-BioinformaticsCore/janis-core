

import copy

from janis_core.utils.errors import UnsupportedError

from janis_core.utils import first_value
from janis_core.workflow.workflow import Workflow, StepNode, InputNode, StepOutputSelector
from janis_core.workflow.workflow import verify_or_try_get_source
from janis_core.types.data_types import is_python_primitive
from janis_core.types import Filename
from janis_core.tool.documentation import (
    InputDocumentation,
    InputQualityType,
)

from .identifiers import get_id_entity
from .identifiers import get_id_path


def get_janis_wf_sources(wf: Workflow, sources: str | list[str]) -> list[InputNode | StepOutputSelector]:
    """
    each source is a workflow input, step output, or complex expression.
    input nodes will all be parsed into janis wf at this stage.
    we can check if the source is an input on the janis wf, then if not, must be a step output.
    """
    out: list[InputNode | StepOutputSelector] = []
    
    if isinstance(sources, str):
        sources = [sources]

    for src in sources:
        # get the wfinp / step output identifier
        identifier = get_id_entity(src)

        # is complex expression?
        if identifier.startswith("$("):
            raise UnsupportedError(
                f"This script can't parse expressions in the step input {step_input}"
            )
        
        # is step output?
        # if referencing step output, that step will have already been parsed into the janis wf
        if get_id_path(src) and get_id_path(src) in wf.step_nodes:
            stp_id = get_id_path(src)
            out_id = get_id_entity(src)
            stp = wf.step_nodes[stp_id]
            selector = stp.get_item(out_id)
            out.append(selector)
        
        # is wf input?
        elif identifier in wf.input_nodes:
            resolved_src = wf[identifier]
            out.append(resolved_src) 
        
        else:
            raise NotImplementedError
        
    return out


def add_step_edges_to_graph(jstep: StepNode, wf: Workflow) -> None:
    connections = jstep.tool.connections
    tinputs = jstep.tool.inputs_map()
    jstep.sources = {}
    
    added_edges = []
    for (k, v) in connections.items():
        
        # static values provided when creating step.
        # janis wants to create a workflow input for these.
        isfilename = isinstance(v, Filename)
        if is_python_primitive(v) or isfilename:
            inp_identifier = f"{jstep.id()}_{k}"
            referencedtype = copy.copy(tinputs[k].intype) if not isfilename else v

            referencedtype.optional = True

            indoc = tinputs[k].doc
            indoc.quality = InputQualityType.configuration

            v = wf.input(
                inp_identifier,
                referencedtype,
                default=v.generated_filename() if isfilename else v,
                doc=indoc,
            )

        if v is None:
            inp_identifier = f"{jstep.id()}_{k}"
            doc = copy.copy(InputDocumentation.try_parse_from(tinputs[k].doc))
            doc.quality = InputQualityType.configuration
            v = wf.input(inp_identifier, tinputs[k].intype, default=v, doc=doc)

        verifiedsource = verify_or_try_get_source(v)
        if isinstance(verifiedsource, list):
            for vv in verifiedsource:
                added_edges.append(jstep._add_edge(k, vv))
        else:
            added_edges.append(jstep._add_edge(k, verifiedsource))

    for e in added_edges:
        si = e.finish.sources[e.ftag] if e.ftag else first_value(e.finish.sources)
        wf.has_multiple_inputs = wf.has_multiple_inputs or si.multiple_inputs
    