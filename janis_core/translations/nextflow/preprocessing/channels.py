

from copy import deepcopy

from janis_core.workflow.workflow import Workflow, InputNode, CommandTool
from janis_core.types import File, Filename, Directory, DataType
from janis_core import translation_utils as utils
from janis_core import settings

from .. import params
from .. import channels
from ..scope import Scope


def register_channels(wf: Workflow) -> None:
    """
    register param(s) for each workflow input. 
    channel(s) may also be registered if necessary.
    """
    scope = Scope()
    do_register_channels(wf, scope)

def do_register_channels(wf: Workflow, scope: Scope):
    channel_inputs = get_channel_inputs_to_register(wf, scope)

    for inp in wf.input_nodes.values():
        if inp.id() in channel_inputs:
            ChannelRegistrationHelper(inp, scope).register()
    
    # repeat for nested workflows (subworkflows)
    for step in wf.step_nodes.values():
        current_scope = deepcopy(scope)
        current_scope.update(step)
        if isinstance(step.tool, Workflow):
            do_register_channels(step.tool, current_scope)


class ChannelRegistrationHelper:
    def __init__(self, inp: InputNode, scope: Scope) -> None:
        self.inp = inp
        self.scope = scope

    @property
    def basetype(self) -> DataType:
        basetype = utils.get_base_type(self.inp.datatype)
        basetype = utils.ensure_single_type(basetype)
        return basetype

    @property
    def is_subworkflow(self) -> bool:
        if self.scope.labels != [settings.translate.nextflow.NF_MAIN_NAME]:
            return True
        return False
    
    @property
    def method(self) -> str:
        if isinstance(self.basetype, File) or isinstance(self.basetype, Filename) or isinstance(self.basetype, Directory):
            method = 'fromPath'
        else: 
            method = 'of'
        return method

    @property
    def source(self) -> str:
        if self.is_subworkflow:
            src = ''
        else:
            param_name = params.getall(self.inp.uuid)[0].name
            # @secondaryarrays
            if utils.is_array_secondary_type(self.inp.datatype):
                src = f'params.{param_name}.flatten()'
            else:
                src = f'params.{param_name}'
        return src

    def register(self) -> None:
        operations = self.get_operations()
        channels.add(
            janis_tag=self.inp.id(),
            method=self.method,
            source=self.source,
            operations=operations,
            janis_dtype=self.inp.datatype,
            janis_uuid=self.inp.uuid,
            define=False if self.is_subworkflow else True
        )

    def get_operations(self) -> str:
        if utils.is_array_secondary_type(self.inp.datatype):
            ops = self.get_operations_secondary_array()
        
        elif utils.is_secondary_type(self.inp.datatype):
            ops = self.get_operations_secondary()
        
        elif isinstance(self.basetype, File) or isinstance(self.basetype, Filename) or isinstance(self.basetype, Directory):
            if self.inp.datatype.is_array():
                ops = self.get_operations_file_array()
            else:
                ops = self.get_operations_generic()
        
        elif self.inp.datatype.is_array():
            ops = self.get_operations_nonfile_array()
        
        else:
            ops = self.get_operations_generic()
        return ops

    def get_operations_secondary_array(self) -> str:
        exts = utils.get_extensions(self.basetype)
        size = len(exts)
        
        ops: str = ''
        ops += f'.collate( {size} )'
        if self.inp.datatype.optional:
            ops += '.ifEmpty( null )'
        return ops

    def get_operations_secondary(self) -> str:
        ops: str = ''
        ops += '.toList()'
        if self.inp.datatype.optional:
            ops += '.ifEmpty( null )'
        return ops

    def get_operations_file_array(self) -> str:
        ops: str = ''
        ops += '.toList()'
        if self.inp.datatype.optional:
            ops += '.ifEmpty( null )'
        return ops

    def get_operations_nonfile_array(self) -> str:
        ops: str = ''
        ops += '.toList()'
        if self.inp.datatype.optional:
            ops += '.ifEmpty( null )'
        return ops

    def get_operations_generic(self) -> str:
        ops: str = ''
        if self.inp.datatype.optional:
            ops += '.ifEmpty( null )'
        return ops
    


def get_channel_inputs_to_register(wf: Workflow, scope: Scope) -> set[str]:
    if scope.labels == [settings.translate.nextflow.NF_MAIN_NAME]:
        items: set[str] = _get_channel_input_ids(wf)
    else:
        items: set[str] = set(wf.connections.keys())
    return items

def _get_channel_input_ids(wf: Workflow) -> set[str]:
    """
    Get the wf inputs for which we will create a nf channel.
    """
    subworkflow_inputs = _get_subworkflow_inputs(wf)
    file_inputs = _get_file_wf_inputs(wf)
    filename_inputs = _get_filename_wf_inputs(wf)
    null_value_to_subworkflow_inputs = _get_null_value_to_subworkflow_inputs(wf)
    scatter_inputs = _get_scatter_wf_inputs(wf)

    channel_inputs: list[InputNode] = []
    for name, inp in wf.input_nodes.items():
        if name in subworkflow_inputs:
            channel_inputs.append(inp)
        elif name in file_inputs or name in filename_inputs or name in scatter_inputs or name in null_value_to_subworkflow_inputs:
            channel_inputs.append(inp)
    
    # final ordering
    return {x.id() for x in channel_inputs}

def _get_subworkflow_inputs(wf: Workflow) -> set[str]:
    # for subworkflows. 
    # for a given subworkflow, ensures each wf input which was specified
    # in the step call becomes a channel.
    subworkflow_inputs: set[str] = set()
    if wf.connections:
        subworkflow_inputs = subworkflow_inputs | set(wf.connections.keys())
    return subworkflow_inputs

def _get_file_wf_inputs(wf: Workflow) -> set[str]:
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

def _get_filename_wf_inputs(wf: Workflow) -> set[str]:
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

def _get_null_value_to_subworkflow_inputs(wf: Workflow) -> set[str]:
    # bit of an edge case. 
    # nextflow doesnt allow null values being passed to subworkflows. 
    # any param which is fed to a subworkflow which has a default=None value
    # must have a channel (ifEmpty(null)) created for this param.
    # the channel is fed to subworkflow, not the param. 
    out: set[str] = set()
    for step in wf.step_nodes.values():
        # subworkflows
        if isinstance(step.tool, Workflow):
            for inp in step.tool.input_nodes.values():
                if inp.datatype.optional and inp.default is None:
                    out.add(inp.id())
                    print('DEV: optional param with null default passed to subworkflow!')
    return out

def _get_scatter_wf_inputs(wf: Workflow) -> set[str]:
    # scattered inputs of steps are fed via channels.
    out: set[str] = set()
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            should_scatter = src.source_map[0].should_scatter
            node = utils.resolve_node(src)
            if should_scatter and isinstance(node, InputNode):
                out.add(node.id())
    return out