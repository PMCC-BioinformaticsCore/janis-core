

from typing import Any

from janis_core.workflow.workflow import Workflow, InputNode
from janis_core.types import File, Array

from . import params
from . import channels
from . import nfgen_utils

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



class ParamChannelRegisterer:
    # horrid name, I know
    def __init__(self, wf: Workflow, scope: list[str]) -> None:
        self.wf = wf
        self.scope = scope

    @property
    def is_subworkflow(self) -> bool:
        if self.scope:
            return True
        return False

    @property
    def channels_to_register_wfinps(self) -> set[str]:
        if self.scope:
            items: set[str] = set(self.wf.connections.keys())
        else:
            items: set[str] = get_channel_input_ids(self.wf)
        return items
    
    @property
    def params_to_register_wfinps(self) -> set[str]:
        if self.scope:
            items: set[str] = {x.id() for x in self.wf.input_nodes.values()} - self.channels_to_register_wfinps
        else:
            items: set[str] = {x.id() for x in self.wf.input_nodes.values()}
        return items

    # @property
    # def params_to_register_toolouts(self) -> set[str]:
    #     if self.scope:
    #         items: set[str] = {x.id() for x in self.wf.input_nodes.values()} - self.channels_to_register_wfinps
    #     else:
    #         items: set[str] = {x.id() for x in self.wf.input_nodes.values()}
    #     return items
    
    def register(self) -> None:
        for inp in self.wf.input_nodes.values():
            if (not params.exists(inp.uuid)):
                
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
            
            else:
                raise NotImplementedError
    
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
                var_name=inp.id(),
                params=params.getall(inp.uuid),
                method=channels.get_channel_method(inp),
                collect=channels.should_collect(inp),
                allow_null=channels.should_allow_null(inp),
                var_scope=self.scope,
                janis_uuid=inp.uuid,
                define=False if self.is_subworkflow else True
            )

    def register_wfinp_secondaries(self, inp: InputNode) -> None:
        # get the extensions. each extension will create individual param. 
        is_channel_input = True if inp.id() in self.channels_to_register_wfinps else False
        is_param_input = True if inp.id() in self.params_to_register_wfinps else False
        exts: list[str] = []
        exts = nfgen_utils.get_extensions(inp.datatype)
        
        # register a param for each individual file
        if is_param_input:
            for ext in exts:
                params.add(
                    var_name=inp.id(),
                    var_scope=self.scope,
                    dtype=inp.datatype,
                    is_channel_input=True,
                    name_override=f'{inp.id()}_{ext}',
                    janis_uuid=inp.uuid
                )

        # register channel for the workflow input
        if is_channel_input:
            channels.add(
                var_name=inp.id(),
                params=params.getall(inp.uuid),
                method='fromPath',
                collect=True,
                allow_null=channels.should_allow_null(inp),
                var_scope=self.scope,
                janis_uuid=inp.uuid,
                define=False if self.is_subworkflow else True
            )

    def register_wfinp_secondaries_array(self, inp: InputNode) -> None:    
        # get the extensions. 
        # each extension will create individual param and individual channel.
        is_param_input = True if inp.id() in self.params_to_register_wfinps else False
        is_channel_input = True if inp.id() in self.channels_to_register_wfinps else False
        basetype = nfgen_utils.get_base_type(inp.datatype)
        exts = nfgen_utils.get_extensions(basetype)
        
        if is_param_input:
            # register a param for each .ext
            for ext in exts:
                params.add(
                    var_name=inp.id(),
                    var_scope=self.scope,
                    dtype=inp.datatype,
                    is_channel_input=True,
                    name_override=f'{inp.id()}_{ext}s',
                    janis_uuid=inp.uuid
                )
            
        if is_channel_input:
            # register a channel for each .ext
            for ext in exts:
                channels.add(
                    var_name=inp.id(),
                    params=params.getall(inp.uuid),
                    method='fromPath',
                    collect=True,
                    allow_null=channels.should_allow_null(inp),
                    var_scope=self.scope,
                    name_override=f'{inp.id()}_{ext}s',
                    janis_uuid=inp.uuid,
                    define=False if self.is_subworkflow else True
                )


# helper functions

def get_channel_input_ids(wf: Workflow) -> set[str]:
    """
    Get the (assumed) true workflow inputs. 
    Assume that a workflow input is an InputNode which:
        - has the 'File' datatype
        - is referenced in a step input
    Everything else are static step inputs, or non-exposed tool inputs. 
    """
    source_inputs = get_source_inputs(wf)
    file_inputs = get_file_wf_inputs(wf)
    scatter_inputs = get_scatter_wf_inputs(wf)

    channel_inputs: list[InputNode] = []
    for name, inp in wf.input_nodes.items():
        if name in source_inputs:
            channel_inputs.append(inp)
        elif name in file_inputs or name in scatter_inputs:
            channel_inputs.append(inp)
    
    # final ordering
    return {x.id() for x in channel_inputs}

def get_source_inputs(wf: Workflow) -> set[str]:
    source_inputs: set[str] = set()
    if wf.connections:
        source_inputs = source_inputs | set(wf.connections.keys())
    return source_inputs

def get_file_wf_inputs(wf: Workflow) -> set[str]:
    # wf inputs with file type
    out: set[str] = set()
    for name, inp in wf.input_nodes.items():
        dtype = nfgen_utils.get_base_type(inp.datatype)
        if isinstance(dtype, File):
            out.add(name)
    return out

def get_scatter_wf_inputs(wf: Workflow) -> set[str]:
    # scattered inputs are always fed via channels
    out: set[str] = set()
    for step in wf.step_nodes.values():
        for src in step.sources.values():
            scatter = src.source_map[0].scatter
            node = nfgen_utils.resolve_node(src)
            if scatter and isinstance(node, InputNode):
                out.add(node.id())
    return out



