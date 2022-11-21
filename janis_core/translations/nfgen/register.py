

from typing import Any

from janis_core.workflow.workflow import Workflow, InputNode
from janis_core.types import File, Array

from . import params
from . import channels
from . import utils

from .params import Param

        
def register_workflow_inputs(wf: Workflow, scope: list[str]) -> None:
    # get the workflow input ids for which channel(s) should be created
    if scope:
        # subworkflow
        channel_input_ids: set[str] = set()
    else:
        # main workflow
        channel_input_ids = utils.get_channel_input_ids(wf)

    # register param(s) for each workflow input. 
    # channel(s) may also be registered if necessary.
    for inp in wf.input_nodes.values():
        if (not params.exists(inp.id(), scope)):
            
            # secondaries
            if isinstance(inp.datatype, File) and inp.datatype.has_secondary_files():
                register_wfinp_secondaries(inp, scope)
                continue
            
            # array secondaries
            elif isinstance(inp.datatype, Array):
                subtype = inp.datatype.subtype()
                if isinstance(subtype, File) and subtype.has_secondary_files():
                    register_wfinp_secondaries_array(inp, scope)
                    continue

            # anything else
            register_wfinp(inp, scope, channel_input_ids)
            continue
        
        else:
            raise NotImplementedError

# workflow inputs
def register_wfinp(    
    inp: InputNode, 
    scope: list[str], 
    channel_input_ids: set[str]
    ) -> None:
    # param
    default: Any = inp.default if inp.default is not None else None
    is_channel_input = True if inp.id() in channel_input_ids else False
    p = params.add(
        ref_name=inp.id(),
        ref_scope=scope,
        dtype=inp.datatype,
        default=default,
        is_wf_input=is_channel_input
    )
    # channel
    if is_channel_input:
        channels.add(
            ref_name=inp.id(),
            params=[p],
            method=channels.get_channel_method(inp),
            collect=channels.should_collect(inp),
            allow_null=channels.should_allow_null(inp),
            ref_scope=scope
        )
    print()

def register_wfinp_secondaries(    
    inp: InputNode, 
    scope: list[str], 
    ) -> None:
    # get the extensions. each extension will create individual param. 
    exts: list[str] = []
    exts.append(inp.datatype.extension)
    exts += inp.datatype.secondary_files()
    exts = [x.split('.')[-1] for x in exts]
    new_params: list[Param] = []
    
    # register a param for each individual file
    for ext in exts:
        p = params.add(
            ref_name=inp.id(),
            ref_scope=scope,
            dtype=inp.datatype,
            is_wf_input=True,
            name_override=f'{inp.id()}_{ext}'
        )
        new_params.append(p)
    
    # register channel for the workflow input
    channels.add(
        ref_name=inp.id(),
        params=new_params,
        method='fromPath',
        collect=True,
        allow_null=channels.should_allow_null(inp),
        ref_scope=scope
    )

def register_wfinp_secondaries_array(    
    inp: InputNode, 
    scope: list[str], 
    ) -> None:
    # get the extensions. 
    # each extension will create individual param and individual channel.
    exts: list[str] = []
    exts.append(inp.datatype.subtype().extension)
    exts += inp.datatype.subtype().secondary_files()
    exts = [x.split('.')[-1] for x in exts]

    for ext in exts:
        # register a param for file array
        p = params.add(
            ref_name=inp.id(),
            ref_scope=scope,
            dtype=inp.datatype,
            is_wf_input=True,
            name_override=f'{inp.id()}_{ext}s'
        )
        
        # register a channel for file array
        channels.add(
            ref_name=inp.id(),
            params=[p],
            method='fromPath',
            collect=True,
            allow_null=channels.should_allow_null(inp),
            ref_scope=scope,
            name_override=f'{inp.id()}_{ext}s'
        )
    print()




### OLD

"""
### DELEGATING FUNCTIONS

def register(
    the_entity: Optional[Workflow | CommandTool]=None,
    the_dict: Optional[dict[str, Any]]=None, 
    sources: Optional[dict[str, Any]]=None, 
    scope: Optional[list[str]]=None, 
    override: bool=True
    ) -> None:

    if not isinstance(scope, list):
        scope = []
    if isinstance(the_entity, Workflow):
        register_params_for_wf_inputs(the_entity, scope, override)
    elif isinstance(the_entity, CommandTool):
        register_params_for_tool(the_entity, scope, override, sources)
    elif the_dict is not None:
        register_params_for_dict(the_dict, scope, override)
    else:
        raise RuntimeError("nothing to register: please supply 'the_entity' or 'the_dict'")

def register_params_for_wf_inputs(
    workflow: Workflow,
    scope: list[str], 
    override: bool
    ) -> None:
    # param_ids = utils.get_channel_input_ids(workflow)
    # param_inputs = utils.items_with_id(list(workflow.input_nodes.values()), param_ids)
    # for inp in param_inputs:
    for inp in workflow.input_nodes.values():
        default: Any = inp.default if inp.default is not None else None
        # if sources is not None and inp.id() in sources:
        #     src = sources[inp.id()]
        #     default = utils.get_source_value(src)
        register(inp, scope=scope, default=default, override=override)

def register_params_for_tool(
    tool: CommandTool, 
    scope: list[str], 
    override: bool,
    sources: Optional[dict[str, Any]]=None, 
    ) -> None:
    # Registers a param tool or workflow inputs.
    # Workflow / tool inputs which are exposed to the user must to be listed
    # as part of the global params object. 
    raise NotImplementedError
    sources = sources if sources is not None else {}
    param_ids = utils.get_param_input_ids(tool, sources=sources)
    param_inputs = utils.items_with_id(tool.inputs(), param_ids)
    for inp in param_inputs:
        default = None
        if sources is not None and inp.id() in sources:
            src = sources[inp.id()]
            default = utils.get_source_value(src)
        register_toolinp_param(inp, scope=scope, default=default, override=override)

def register_params_for_dict(
    the_dict: dict[str, Any],
    scope: list[str],
    override: bool
    ) -> None:
    for k, v in the_dict.items():
        if (not params.exists(k, scope)) or (params.exists(k, scope) and override):
            add(
                ref_name=k,
                ref_scope=scope,
                default=v,
                is_wf_input=False
            )


"""