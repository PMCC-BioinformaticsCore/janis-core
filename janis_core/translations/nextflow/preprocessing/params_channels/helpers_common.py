

from typing import Any
from collections import Counter

from janis_core.workflow.workflow import Workflow, InputNode, CommandTool
from janis_core.types import File, Filename, Directory
from janis_core import translation_utils as utils

from ... import params
from ... import trace
from ... import task_inputs


def get_linkable_params(wf: Workflow, sources: dict[str, Any]) -> dict[str, params.Param]:
    all_inputs = get_all_workflow_inputs(wf)
    input_nodes = [wf.input_nodes[x] for x in all_inputs]

    out: dict[str, params.Param] = {}
    for inp in input_nodes:
        if inp.id() not in sources:
            continue
        else:
            src = sources[inp.id()]
            node = utils.resolve_node(src)
            if isinstance(node, InputNode):
                if params.exists(node.uuid):
                    param = params.get(node.uuid)
                    out[inp.uuid] = param
    return out

def get_all_workflow_inputs(wf: Workflow) -> set[str]:
    true_inputs = get_task_referenced_wf_inputs(wf)
    file_inputs = get_file_wf_inputs(wf)
    non_optional_inputs = get_non_optional_wf_inputs(wf)
    scatter_inputs = get_scatter_wf_inputs(wf)
    final_inputs = true_inputs | file_inputs | non_optional_inputs | scatter_inputs
    return final_inputs

def get_task_referenced_wf_inputs(wf: Workflow) -> set[str]:
    """
    identifies which workflow InputNodes are referenced in minimal task inputs. 
    """
    out: set[str] = set()

    for step in wf.step_nodes.values():
        identifier = step.tool.id()
        minimal_task_inputs = task_inputs.get(identifier)

        for tinput_id in minimal_task_inputs:
            if tinput_id in step.sources:
                src = step.sources[tinput_id]
                refs = _get_input_node_refs(src, wf)
                out = out | refs

    return out

def _get_input_node_refs(src: Any, wf: Workflow) -> set[str]:
    """
    identifies the number of times each workflow input is referenced in step inputs. 
    to facilitate: 
        - for each step, each step input expression is explored
        - any references to workflow inputs in the expression are noted. 
    """
    refs = trace.trace_referenced_variables(src, wf)
    refs = _filter_to_input_node_refs(refs, wf)
    return refs

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
        if isinstance(basetype, (File, Filename, Directory)):
            out.add(name)
        # file pairs
        elif basetype.name() in ['FastqPair', 'FastqGzPair']:
            out.add(name)
    return out

def get_non_optional_wf_inputs(wf: Workflow) -> set[str]:
    # wf inputs with file type
    out: set[str] = set()
    for name, inp in wf.input_nodes.items():
        if not inp.datatype.optional:
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





def get_filename_wf_inputs_dep(wf: Workflow) -> set[str]:
    """
    Edge case!
    ToolInputs which have Filename DataType may require channel.
    
    For a ToolInput which uses InputSelector:
        - Assume it derives name using the InputSelector (another ToolInput)
        - Therefore ToolInput is internal to the (future) process
        - Don't create channel
    
    For a ToolInput which does not use InputSelector:
        - String value does actually need to be supplied to ToolInput
        - Therefore param or channel needed to feed value to (future) process input
        - Create channel  (because Filenames move similarly to Files in nextflow workflow)

    [eg InputSelector]

    inside BwaMem_SamToolsView:
        ToolInput(
            "outputFilename",
            Filename(prefix=InputSelector("sampleName"), extension=".bam"),
            position=8,
            shell_quote=False,
            prefix="-o",
            doc="output file name [stdout]",
        ),

    [eg no InputSelector]

    step call:
        BcfToolsNorm(
            vcf=self.sortSomatic1.out,
            reference=self.reference,
            outputType="v",
            outputFilename="normalised.vcf",
        ),

    inside BcfToolsNorm:
        ToolInput(
            "outputFilename",
            Filename(extension=".vcf.gz"),
            prefix="-o",
            doc="--output: When output consists of a single stream, "
            "write it to FILE rather than to standard output, where it is written by default.",
        ),
    """
    out: set[str] = set()
    for step in wf.step_nodes.values():
        
        # CommandTools
        if isinstance(step.tool, CommandTool):
            # get all tool inputs with Filename type
            filename_inputs = [x for x in step.tool.inputs() if isinstance(x.input_type, Filename)]
            # if fed value from wf input, mark wf input for channel creation
            for inp in filename_inputs:
                if inp.id() in step.sources:
                    src = step.sources[inp.id()]
                    node = utils.resolve_node(src)
                    if isinstance(node, InputNode):
                        out.add(node.id())
    return out


