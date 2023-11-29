
from typing import Optional

from janis_core import CommandToolBuilder, CodeTool
from janis_core import translation_utils as utils
from janis_core.messages import load_loglines, ErrorCategory

from ...model.files import NFFile
from ...model.files import NFImportsBlock
from ...model.files import NFFunctionsBlock
from ...model.files import NFMessageBlock
from ...model.process import NFProcess

from janis_core.translations.common import trace


get_primary_files_code = """\
def get_primary_files(var, element_count) {
    def primary_files = []
    var.eachWithIndex {item, index -> 
        if (index % element_count == 0) {
            primary_files.add(item)
        }
    }
    return primary_files
}"""


def generate_file_process(process: NFProcess, tool: CommandToolBuilder | CodeTool) -> NFFile:
    """generates nextflow file for nextflow process derived from CommandToolBuilder"""
    nf_file = NFFile(subtype='process', name=process.name)

    # item: messages
    message_items = gen_message_items_for_process_file(tool, process)
    if message_items:
        nf_file.items += message_items

    # groovy library imports & groovy functions used in process
    # item: imports
    imports_item = gen_imports_for_process_file(tool)
    if imports_item:
        nf_file.items.append(imports_item)
    
    # item: functions
    functions_item = gen_functions_for_process_file(tool)
    if functions_item:
        nf_file.items.append(functions_item)

    # item: process
    nf_file.items.append(process)

    return nf_file

def gen_message_items_for_process_file(tool: CommandToolBuilder | CodeTool, process: NFProcess) -> list[NFMessageBlock]:
    message_items = []

    # fallbacks
    loglines = load_loglines(category=ErrorCategory.FALLBACK, tool_uuid=tool.uuid)
    if loglines:
        message_items.append(NFMessageBlock(
            heading='ERROR: PARSING FALLBACKS',
            lines=[f'// {x.message}' for x in loglines]
        ))
    
    # datatypes
    loglines = load_loglines(category=ErrorCategory.DATATYPE, tool_uuid=tool.uuid)
    if loglines:
        # TODO handle updated variable names here?
        message_items.append(NFMessageBlock(
            heading='WARNING: DATATYPES',
            lines=[f'// {x.message}' for x in loglines]
        ))
    
    # version
    loglines = load_loglines(category=ErrorCategory.VERSION, tool_uuid=tool.uuid)
    if loglines:
        message_items.append(NFMessageBlock(
            heading='WARNING: VERSION',
            lines=[f'// {x.message}' for x in loglines]
        ))
    
    # experimental
    loglines = load_loglines(category=ErrorCategory.EXPERIMENTAL, tool_uuid=tool.uuid)
    if loglines:
        message_items.append(NFMessageBlock(
            heading='WARNING: EXPERIMENTAL FEATURES',
            lines=[f'// {x.message}' for x in loglines]
        ))
    
    # untranslated expressions
    loglines = load_loglines(category=ErrorCategory.SCRIPT, tool_uuid=tool.uuid)
    if loglines:
        messages = [x.message for x in loglines]
        filtered = []
        
        # filter to messages where __TOKEN__ still is in the process text
        for msg in messages:
            token = msg.split(' = ', 1)[0]
            if token in process.get_string():
                filtered.append(msg)
        
        # create block & append
        if filtered:
            message_items.append(NFMessageBlock(
                heading='ERROR: UNTRANSLATED EXPRESSIONS',
                lines=filtered
            ))

    return message_items

def gen_imports_for_process_file(tool: CommandToolBuilder | CodeTool) -> Optional[NFImportsBlock]:
    imports: list[str] = []
    declarations: list[str] = []

    if _should_add_math(tool):
        imports.append('import java.lang.Math')
    
    if _should_add_json_slurper(tool):
        imports.append('import groovy.json.JsonSlurper')
        declarations.append('jsonSlurper = new JsonSlurper()')
    
    if imports:
        # return NFImportsBlock(methods, imports, declarations)
        return NFImportsBlock(imports, declarations)
    else:
        return None

def _should_add_math(tool: CommandToolBuilder | CodeTool) -> bool:
    if isinstance(tool, CommandToolBuilder):
        entity_counts = {}
        # runtime
        for item in [tool._cpus, tool._disk, tool._memory, tool._time]:
            if item is not None:
                entity_counts = entity_counts | trace.trace_entity_counts(item, tool)
        # inputs
        if tool._inputs is not None:
            for tinput in tool._inputs:
                entity_counts = entity_counts | trace.trace_entity_counts(tinput, tool)
        # arguments
        if tool._arguments is not None:
            for targ in tool._arguments:
                entity_counts = entity_counts | trace.trace_entity_counts(targ, tool)
        # outputs
        if tool._outputs is not None:
            for tout in tool._outputs:
                entity_counts = entity_counts | trace.trace_entity_counts(tout, tool)
        
        for op in ["CeilOperator", "FloorOperator", "RoundOperator"]:
            if op in entity_counts:
                return True
        return False
    return False


def _should_add_json_slurper(tool: CommandToolBuilder | CodeTool) -> bool:
    if isinstance(tool, CommandToolBuilder):
        for toutput in tool.outputs():
            entity_counts = trace.trace_entity_counts(toutput.selector, tool)
            if 'ReadJsonOperator' in entity_counts:
                return True
    return False

def gen_functions_for_process_file(tool: CommandToolBuilder | CodeTool) -> Optional[NFFunctionsBlock]:
    funcs: list[str] = []

    if _should_add_get_primary_files(tool):
        funcs.append(get_primary_files_code)
    
    if funcs:
        return NFFunctionsBlock(funcs)
    else:
        return None
    
def _should_add_get_primary_files(tool: CommandToolBuilder | CodeTool) -> bool:
    for tinput in tool.tool_inputs():
        if utils.is_secondary_array_type(tinput.intype):
            return True
    return False