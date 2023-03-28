

from typing import Any

from janis_core.workflow.workflow import Workflow, InputNode, CommandTool
from janis_core.types import File, Filename
from janis_core import translation_utils as utils

from ... import params


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

# def get_true_workflow_inputs(wf: Workflow) -> set[str]:
#     """
#     Get the wf inputs for which we will create a nf channel in the main wf.
#     inputs which are:
#     - files
#     - filenames (in some cases)
#     - scattered on 
#     will become nextflow channels. all other inputs will become global scope params. 
#     """
#     file_inputs = get_file_wf_inputs(wf)
#     filename_inputs = get_filename_wf_inputs(wf)
#     scatter_inputs = get_scatter_wf_inputs(wf)
#     final_inputs = file_inputs | filename_inputs | scatter_inputs
#     return final_inputs

def get_all_workflow_inputs(wf: Workflow) -> set[str]:
    true_inputs = get_true_workflow_inputs(wf)
    file_inputs = get_file_wf_inputs(wf)
    filename_inputs = get_filename_wf_inputs(wf)
    scatter_inputs = get_scatter_wf_inputs(wf)
    final_inputs = true_inputs | file_inputs | filename_inputs | scatter_inputs
    return final_inputs
    
def get_true_workflow_inputs(wf: Workflow) -> set[str]:
    """
    workflow inputs are the set of InputNodes which have are referred to in 
    step.sources among all steps in the workflow
    """
    out: set[str] = set()
    for step in wf.step_nodes.values():
        for tinput_id, src in step.sources.items():
            src_node = utils.resolve_node(src)
            if isinstance(src_node, InputNode):  # anything else is a static value, ie 1 or 'conservative' etc
                if tinput_id in step.tool.connections:
                    connection = step.tool.connections[tinput_id]
                    connection_node = utils.resolve_node(connection)
                    if isinstance(connection_node, InputNode):
                        out.add(src_node.id())
                else:
                    out.add(src_node.id())

    return out

def get_file_wf_inputs(wf: Workflow) -> set[str]:
    # wf inputs with file type are fed via channels.
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
        - Create channel  (because Filenames move similarly to Files in workflow)

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

def get_scatter_wf_inputs(wf: Workflow) -> set[str]:
    # scattered inputs of steps are fed via channels.
    out: set[str] = set()
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            should_scatter = src.source_map[0].should_scatter
            node = utils.resolve_node(src)
            if should_scatter and isinstance(node, InputNode):
                out.add(node.id())
    return out


