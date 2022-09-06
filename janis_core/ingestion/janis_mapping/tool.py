

from datetime import datetime
from typing import Optional

from janis_core import CommandToolBuilder
from janis_core import ToolInput
from janis_core import ToolOutput
from janis_core import ToolMetadata

from galaxy2janis.runtime.dates import JANIS_DATE_FMT
from galaxy2janis.gx.gxtool.metadata import ToolXMLMetadata
from galaxy2janis.entities.tool import Tool
from galaxy2janis.gx.command.components import (
    InputComponent,
    Flag,
    Option,
    OutputComponent
)

from .general import to_janis_datatype
from .general import to_janis_selector


### MODULE EXPORTS

def to_janis_tool(internal: Tool) -> CommandToolBuilder:
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
    return CommandToolBuilder(
        tool=internal.tag,
        base_command=internal.base_command,
        inputs=[to_janis_tool_input(inp) for inp in internal.inputs],
        outputs=[to_janis_tool_output(out) for out in internal.outputs],
        container=internal.container,                   # TODO check None type is ok
        version=internal.metadata.version,
        metadata=to_janis_metadata(internal.metadata),
        doc=internal.metadata.help
    )


### HELPER METHODS ###

def to_janis_tool_input(internal_inp: InputComponent) -> ToolInput:
    """ 
    maps internal model tool input to janis model tool input
    missing the following (unnessesary):
        prefix_applies_to_all_elements
        presents_as
        secondaries_present_as
        separator
        shell_quote
        localise_file
    """
    # these should be the janis ToolInput defaults
    prefix: Optional[str] = None
    separate: bool = True

    # derive special attributes in case of flag tool input
    if isinstance(internal_inp, Flag):
        prefix = internal_inp.prefix
    
    # derive special attributes in case of option tool input
    elif isinstance(internal_inp, Option):
        if internal_inp.delim != ' ':
            prefix = f'{internal_inp.prefix}{internal_inp.delim}'
            separate = False
        
    return ToolInput(
        tag=internal_inp.tag,
        input_type=to_janis_datatype(internal_inp),
        position=internal_inp.cmd_pos,
        prefix=prefix,
        separate_value_from_prefix=separate,
        default=internal_inp.default_value,
        doc=internal_inp.docstring
    )

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

def to_janis_metadata(internal_meta: ToolXMLMetadata) -> ToolMetadata:
    """ maps internal model tool metadata to janis model tool metadata"""
    return ToolMetadata(
        short_documentation=internal_meta.description,
        keywords=[],
        contributors=_get_contributors(internal_meta),
        dateCreated=datetime.today().strftime(JANIS_DATE_FMT),
        dateUpdated=datetime.today().strftime(JANIS_DATE_FMT),
        version=internal_meta.version,  
        doi=internal_meta.get_doi_citation(),
        citation=internal_meta.get_main_citation(),
        documentationUrl=None,
        documentation=f'"""{internal_meta.help}"""'
    )

def _get_contributors(internal_meta: ToolXMLMetadata) -> list[str]:
    contributors: list[str] = ['gxtool2janis']
    if internal_meta.owner:
        contributors += [f'Wrapper owner: galaxy toolshed user {internal_meta.owner}']
    if internal_meta.creator:
        contributors += [f'Wrapper creator: {internal_meta.creator}']
    return contributors