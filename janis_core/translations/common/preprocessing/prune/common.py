
from collections import Counter

from janis_core.workflow.workflow import Workflow, InputNode
from janis_core.types import File, Filename
from janis_core import translation_utils as utils
from janis_core.translations.common import trace



def get_true_workflow_inputs(wf: Workflow) -> set[str]:
    true_inputs = get_referenced_workflow_inputs(wf)
    file_inputs = get_file_wf_inputs(wf)
    filename_inputs = get_filename_wf_inputs(wf)
    scatter_inputs = get_scatter_wf_inputs(wf)
    final_inputs = true_inputs | file_inputs | filename_inputs | scatter_inputs
    return final_inputs

def get_referenced_workflow_inputs(wf: Workflow) -> set[str]:
    """
    identifies which workflow InputNodes are 'true' workflow inputs, and which 
    are simply static values provided in the step call. 

    'true' workflow inputs:
    1. if an InputNode is referenced in > 1 step inputs among all steps
    2. if an InputNode is referenced in 1 step input among all steps, but has no default value

    other inputs are either not referenced in workflow (will be ignored), or 
    are referenced in 1 step input with a default (static)
    """
    out: set[str] = set()
    counts = _get_input_node_reference_counts(wf)
    for input_node_id, count in counts.items():
        if count > 1:
            out.add(input_node_id)
        elif count == 1:
            inp = wf.input_nodes[input_node_id]
            if inp.default is None:
                out.add(input_node_id)

    return out

def _get_input_node_reference_counts(wf: Workflow) -> dict[str, int]:
    """
    identifies the number of times each workflow input is referenced in step inputs. 
    to facilitate: 
        - for each step, each step input expression is explored
        - any references to workflow inputs in the expression are noted. 
    """
    counts = Counter()

    for step in wf.step_nodes.values():
        for tinput_id, src in step.sources.items():
            refs = trace.trace_referenced_variables(src, wf)
            refs = _filter_to_input_node_refs(refs, wf)
            for ref in refs:
                counts[ref] += 1
    return counts

def _filter_to_input_node_refs(references: set[str], wf: Workflow) -> set[str]:
    out: set[str] = set()
    for ref in references:
        if ref in wf.input_nodes:
            out.add(ref)
    return out

def get_file_wf_inputs(wf: Workflow) -> set[str]:
    # wf inputs with file type
    out: set[str] = set()
    for name, inp in wf.input_nodes.items():
        basetype = utils.get_base_type(inp.datatype)
        basetype = utils.ensure_single_type(basetype)
        # main file types
        if isinstance(basetype, File):
            out.add(name)
        # file pairs
        elif basetype.name() in ['FastqPair', 'FastqGzPair']:
            out.add(name)
    return out

def get_filename_wf_inputs(wf: Workflow) -> set[str]:
    # wf inputs with filename type
    out: set[str] = set()
    for name, inp in wf.input_nodes.items():
        basetype = utils.get_base_type(inp.datatype)
        basetype = utils.ensure_single_type(basetype)
        if isinstance(basetype, Filename):
            out.add(name)
    return out

def get_scatter_wf_inputs(wf: Workflow) -> set[str]:
    # wf inputs which are scattered on in 1+ steps
    out: set[str] = set()
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            should_scatter = src.source_map[0].should_scatter
            node = utils.resolve_node(src)
            if should_scatter and isinstance(node, InputNode):
                out.add(node.id())
    return out

