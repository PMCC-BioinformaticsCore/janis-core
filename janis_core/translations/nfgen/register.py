

import os
from typing import Any

from janis_core.workflow.workflow import Workflow, InputNode, CommandTool
from janis_core.types import File, Array, Filename
from janis_core import PythonTool

from . import params
from . import channels
from . import nfgen_utils
from . import naming
from . import settings

from copy import deepcopy
    

def register_params_channels(wf: Workflow, scope: list[str]) -> None:
    # register param(s) for each workflow input. 
    # channel(s) may also be registered if necessary.

    # handle this workflow
    handler = ParamChannelRegisterer(wf, scope)
    handler.register()
    
    # handle nested workflows (subworkflows)
    for step in wf.step_nodes.values():
        if isinstance(step.tool, Workflow):
            current_scope = deepcopy(scope)
            current_scope.append(step.id())
            register_params_channels(step.tool, scope=current_scope)


def get_code_file_path(tool: PythonTool) -> str:
    basedir = settings.BASE_OUTDIR
    subfolder = settings.CODE_FILES_OUTDIR
    filename = tool.id()
    filepath = os.path.join(basedir, subfolder, filename)
    filepath += '.py'
    return filepath


class ParamChannelRegisterer:
    # sorry about horrid name
    def __init__(self, wf: Workflow, scope: list[str]) -> None:
        self.wf = wf
        self.scope = scope

    @property
    def is_subworkflow(self) -> bool:
        if self.scope != [settings.NF_MAIN_NAME]:
            return True
        return False

    @property
    def channels_to_register_wfinps(self) -> set[str]:
        if self.scope == [settings.NF_MAIN_NAME]:
            items: set[str] = get_channel_input_ids(self.wf)
        else:
            items: set[str] = set(self.wf.connections.keys())
        return items
    
    @property
    def params_to_register_wfinps(self) -> set[str]:
        if self.scope == [settings.NF_MAIN_NAME]:
            items: set[str] = {x.id() for x in self.wf.input_nodes.values()}
        else:
            items: set[str] = {x.id() for x in self.wf.input_nodes.values()} - self.channels_to_register_wfinps
        return items
    
    def register(self) -> None:
        self.register_wf_inputs()
        self.register_python_tools()

    def register_python_tools(self) -> None:
        # A param will be registered for the code_file of each PythonTool.
        for step in self.wf.step_nodes.values():
            current_scope = deepcopy(self.scope)
            current_scope.append(step.id())
            if isinstance(step.tool, PythonTool):
                default = get_code_file_path(step.tool)
                params.add(
                    var_name='code_file',
                    var_scope=current_scope,
                    dtype=File(),
                    default=default,
                    is_channel_input=False,
                    janis_uuid=None,
                )



    def register_wf_inputs(self) -> None:
        # registers param for each wf input which requires a param.
        # registers channel for each wf input which requires a channel.
        for inp in self.wf.input_nodes.values():
            # secondaries
            if isinstance(inp.datatype, File) and inp.datatype.has_secondary_files():
                self.register_wfinp_secondaries(inp)
                continue
            
            # array secondaries
            elif isinstance(inp.datatype, Array):
                basetype = nfgen_utils.get_base_type(inp.datatype)
                if isinstance(basetype, File) and basetype.has_secondary_files():
                    self.register_wfinp_secondaries_array(inp)
                    continue

            # anything else
            self.register_wfinp(inp)
            continue
    
    def register_wfinp(self, inp: InputNode) -> None:
        is_channel_input = True if inp.id() in self.channels_to_register_wfinps else False
        is_param_input = True if inp.id() in self.params_to_register_wfinps else False
        default: Any = inp.default if inp.default is not None else None

        # register a param for the wf input
        if is_param_input:
            params.add(
                var_name=inp.id(),
                var_scope=self.scope,
                dtype=inp.datatype,
                default=default,
                is_channel_input=is_channel_input,
                janis_uuid=inp.uuid,
            )
        
        # register a channel for the wf input if required
        if is_channel_input:
            channels.add(
                janis_tag=inp.id(),
                params=params.getall(inp.uuid),
                method=channels.get_channel_method(inp),
                collect=channels.should_collect(inp),
                allow_null=channels.should_allow_null(inp),
                janis_uuid=inp.uuid,
                define=False if self.is_subworkflow else True
            )

    def register_wfinp_secondaries(self, inp: InputNode) -> None:
        # get the extensions. each extension will create individual param. 
        is_channel_input = True if inp.id() in self.channels_to_register_wfinps else False
        is_param_input = True if inp.id() in self.params_to_register_wfinps else False
        names: list[str] = []
        names = naming.get_varname_secondaries(inp.datatype)
        
        # register a param for each individual file
        if is_param_input:
            for name in names:
                params.add(
                    var_name=inp.id(),
                    var_scope=self.scope,
                    dtype=inp.datatype,
                    is_channel_input=True,
                    # name_override=name,
                    name_override=f'{inp.id()}_{name}',
                    janis_uuid=inp.uuid
                )

        # register channel for the workflow input
        if is_channel_input:
            channels.add(
                janis_tag=inp.id(),
                params=params.getall(inp.uuid),
                method='fromPath',
                collect=True,
                allow_null=channels.should_allow_null(inp),
                janis_uuid=inp.uuid,
                define=False if self.is_subworkflow else True
            )

    def register_wfinp_secondaries_array(self, inp: InputNode) -> None:    
        # get the extensions. 
        # each extension will create individual param and individual channel.
        is_param_input = True if inp.id() in self.params_to_register_wfinps else False
        is_channel_input = True if inp.id() in self.channels_to_register_wfinps else False
        basetype = nfgen_utils.get_base_type(inp.datatype)
        names = naming.get_varname_secondaries(basetype)

        for name in names:
            secondary_params: list[params.Param] = []
            # should we register a param?
            if is_param_input:
                new_param = params.add(
                    var_name=inp.id(),
                    var_scope=self.scope,
                    dtype=inp.datatype,
                    is_channel_input=True,
                    name_override=f'{inp.id()}_{name}s',
                    janis_uuid=inp.uuid
                )
                secondary_params.append(new_param)
            
            # should we register a channel?
            if is_channel_input:
                channels.add(
                    janis_tag=inp.id(),
                    params=secondary_params,
                    method='fromPath',
                    collect=True,
                    allow_null=channels.should_allow_null(inp),
                    name_override=f'{inp.id()}_{name}s',
                    janis_uuid=inp.uuid,
                    define=False if self.is_subworkflow else True
                )


# helper functions

def get_channel_input_ids(wf: Workflow) -> set[str]:
    """
    Get the wf inputs for which we will create a nf channel.
    """
    subworkflow_inputs = get_subworkflow_inputs(wf)
    file_inputs = get_file_wf_inputs(wf)
    filename_inputs = get_filename_wf_inputs(wf)
    scatter_inputs = get_scatter_wf_inputs(wf)

    channel_inputs: list[InputNode] = []
    for name, inp in wf.input_nodes.items():
        if name in subworkflow_inputs:
            channel_inputs.append(inp)
        elif name in file_inputs or name in filename_inputs or name in scatter_inputs:
            channel_inputs.append(inp)
    
    # final ordering
    return {x.id() for x in channel_inputs}

def get_subworkflow_inputs(wf: Workflow) -> set[str]:
    # for subworkflows. 
    # for a given subworkflow, ensures each wf input which was specified
    # in the step call becomes a channel.
    subworkflow_inputs: set[str] = set()
    if wf.connections:
        subworkflow_inputs = subworkflow_inputs | set(wf.connections.keys())
    return subworkflow_inputs

def get_file_wf_inputs(wf: Workflow) -> set[str]:
    # wf inputs with file type are fed via channels.
    out: set[str] = set()
    for name, inp in wf.input_nodes.items():
        dtype = nfgen_utils.get_base_type(inp.datatype)
        if isinstance(dtype, File):
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
                    node = nfgen_utils.resolve_node(src)
                    if isinstance(node, InputNode):
                        out.add(node.id())
    return out

def get_scatter_wf_inputs(wf: Workflow) -> set[str]:
    # scattered inputs of steps are fed via channels.
    out: set[str] = set()
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            scatter = src.source_map[0].scatter
            node = nfgen_utils.resolve_node(src)
            if scatter and isinstance(node, InputNode):
                out.add(node.id())
    return out



