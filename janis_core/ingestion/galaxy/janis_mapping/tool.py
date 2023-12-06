

from datetime import datetime
from typing import Optional

from janis_core import CommandToolBuilder
from janis_core import ToolInput
from janis_core import ToolOutput
from janis_core import ToolMetadata

from janis_core.ingestion.galaxy.runtime.dates import JANIS_DATE_FMT
from janis_core.ingestion.galaxy.gxtool.model import XMLMetadata
from janis_core.ingestion.galaxy.internal_model.tool import ITool
from janis_core.ingestion.galaxy.gxtool.command.components import (
    InputComponent,
    Flag,
    Option,
    OutputComponent
)

from .general import to_janis_datatype
from .general import to_janis_selector


### MODULE EXPORTS

def to_janis_tool(internal: ITool) -> CommandToolBuilder:
    """
    maps internal model tool to janis model tool
    missing the following (unnessesary):
        friendly_name
        arguments
        env_vars
        tool_module
        tool_provider
        cpus
        memory
        time
        disk
        directories_to_create
        files_to_create
    """
    outs = [to_janis_tool_output(gout) for gout in internal.outputs]
    inps = [to_janis_tool_input(ginp) for ginp in internal.inputs]
    return CommandToolBuilder(
        tool=internal.tag,
        base_command=internal.base_command,
        inputs=inps,
        outputs=outs,
        container=internal.container,                   # TODO check None type is ok
        version=internal.metadata.version,
        metadata=to_janis_metadata(internal.metadata),
        files_to_create=to_janis_files_to_create(internal),  # type: ignore
        doc=internal.metadata.help
    )

def to_janis_files_to_create(internal: ITool) -> dict[str, str]:
    files_to_create: dict[str, str] = {}
    if internal.configfiles:
        for configfile in internal.configfiles:
            files_to_create[configfile.varname] = configfile.contents
    if internal.scripts:
        for script in internal.scripts:
            files_to_create[script.filename] = script.contents
    return files_to_create

def to_janis_tool_input(internal_inp: InputComponent) -> ToolInput:
    """ 
    maps internal model tool input to janis model tool input
    missing the following (some are unnessesary, others may need to be implemented):
        prefix_applies_to_all_elements
        presents_as
        secondaries_present_as
        separator
        shell_quote
        localise_file
    """
    # these should be the janis ToolInput defaults
    prefix: Optional[str] = None
    separate: Optional[bool] = None

    # derive special attributes in case of flag tool input
    if isinstance(internal_inp, Flag):
        prefix = internal_inp.prefix
    
    # derive special attributes in case of option tool input
    elif isinstance(internal_inp, Option):
        if internal_inp.separator != ' ':
            prefix = f'{internal_inp.prefix}{internal_inp.separator}'
            separate = False
        else:
            prefix = internal_inp.prefix
        
    tinp = ToolInput(
        tag=internal_inp.tag,
        input_type=to_janis_datatype(internal_inp),
        position=internal_inp.cmd_pos,
        prefix=prefix,
        separate_value_from_prefix=separate,
        default=internal_inp.default_value,
        doc=internal_inp.docstring
    )
    return tinp

def to_janis_tool_output(internal_out: OutputComponent) -> ToolOutput:
    """ 
    maps internal model tool output to janis model tool output
    missing the following (unnessesary):
        glob
        presents_as
        secondaries_present_as
        _skip_output_quality_check
    """
    return ToolOutput(
        tag=internal_out.tag,
        output_type=to_janis_datatype(internal_out),
        selector=to_janis_selector(internal_out),      
        doc=internal_out.docstring
    )

def to_janis_metadata(internal_meta: XMLMetadata) -> ToolMetadata:
    """ maps internal model tool metadata to janis model tool metadata"""
    return ToolMetadata(
        short_documentation=internal_meta.description,
        keywords=[],
        contributors=_get_contributors(internal_meta),
        dateCreated=datetime.today().strftime(JANIS_DATE_FMT),
        dateUpdated=datetime.today().strftime(JANIS_DATE_FMT),
        version=internal_meta.version,  
        doi=internal_meta.doi_citation,
        citation=internal_meta.main_citation,
        documentationUrl=None,
        documentation=f'"""{internal_meta.help}"""'
    )

def _get_contributors(internal_meta: XMLMetadata) -> list[str]:
    contributors: list[str] = ['gxtool2janis']
    if internal_meta.owner:
        contributors += [f'Wrapper owner: galaxy toolshed user {internal_meta.owner}']
    if internal_meta.creator:
        contributors += [f'Wrapper creator: {internal_meta.creator}']
    return contributors