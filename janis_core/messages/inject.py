
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from janis_core import CommandToolBuilder, CodeTool, WorkflowBuilder

from Levenshtein import distance as levenshtein_distance
import regex as re
from typing import Optional, Any, Tuple

from janis_core import settings
from .logfile import ErrorCategory
from .main import load_loglines

### COMMENTERS ###

commenter_map = {
    'nextflow': '//',
    'cwl': '#',
    'wdl': '#',
}

### PUBLIC ###

def inject_messages_tool(internal: CommandToolBuilder | CodeTool, translated: str) -> str:
    """
    Injects messages into the translated workflow file.
    Messages need to be logged during ingest to use this feature.
    
    When logging messages during ingest you need to provide the unique tool_uuid. 
    This enables the right messages to be written to the right output files. 
    Examples:
        log_message(tool.uuid, msg=="couldn't parse datatype...", category=ErrorCategory.DATATYPE)

    """
    messages = _gen_block(internal, translated)
    if messages:
        translated = messages + '\n\n' + translated
    return translated

def inject_messages_workflow(internal: WorkflowBuilder, translated: str) -> str:
    """
    Injects messages into the translated workflow file.
    Default behaviour is to dump workflow messages at the top of the file, and step specific  
    messages above the relevant step call. 

    If split_categories=True, the messages at the top of the file are split into sections based on their category.
    There are three main sections - general, inputs, outputs. 
    If you want to use this feature, make sure the ingest unit is providing the subsection=[section] argument 
    when logging messages. Eg:
        log_message(tool_uuid, msg, category, subsection='general')
        log_message(tool_uuid, msg, category, subsection='inputs')
        log_message(tool_uuid, msg, category, subsection='outputs')
    
    Relies on workflow step names being almost equivalent to the step call in the translated file.
    """
    if settings.messages.USE_SUBSECTIONS:
        translated = _inject_workflow_header_messages_subsections(internal, translated)
    else:
        translated = _inject_workflow_header_messages_default(internal, translated)
    translated = _inject_workflow_step_messages(internal, translated)
    return translated


### PRIVATE ###

def _inject_workflow_header_messages_default(internal: WorkflowBuilder, translated: str) -> str:
    messages = _gen_block(internal, translated, subsections=[None, 'general', 'inputs', 'outputs'])
    if messages:
        translated = messages + '\n\n' + translated
    return translated

def _inject_workflow_header_messages_subsections(internal: WorkflowBuilder, translated: str) -> str:
    general = _gen_block(internal, translated, heading='general', subsections=['general'])
    inputs = _gen_block(internal, translated, heading='inputs', subsections=['inputs'])
    outputs = _gen_block(internal, translated, heading='outputs', subsections=['outputs'])
    messages = '\n\n'.join([x for x in [general, inputs, outputs] if x is not None])
    if messages:
        translated = messages + '\n\n' + translated
    return translated

def _inject_workflow_step_messages(internal: WorkflowBuilder, translated: str) -> str:
    """
    This is horrendous.
    For each step: 
    - get messages
    - if messages, find step call
    - if found, get indent
    - inject message block above step call with correct indent
    """

    for sname in internal.step_nodes.keys():
        # get message block for this step if step has messages
        messages = _gen_block(internal, translated, subsections=[sname])
        if not messages:
            continue
        
        # get step location in translated text
        loc = _get_step_loc(sname, translated)
        if loc is None:
            raise NotImplementedError
        
        # get indent for the step call - awful in general
        indent = translated[:loc].split('\n')[-1]
        lines = messages.split('\n')
        first = lines[0]
        other = [indent + ln for ln in lines[1:]]
        lines = [first] + other
        messages = '\n'.join(lines)

        # inject message block above step call with correct indent
        translated = translated[:loc] + messages + '\n' + indent + translated[loc:]
    
    return translated

def _get_step_loc(step_name: str, translated: str) -> Optional[int]:

    # get loc, step call name for all step calls in translated text
    if settings.translate.DEST == 'cwl':
        matches = _get_step_call_starts_cwl(translated)
    elif settings.translate.DEST == 'nextflow':
        matches = _get_step_call_starts_nextflow(translated)
    elif settings.translate.DEST == 'wdl':
        matches = _get_step_call_starts_wdl(translated)
    else:
        raise NotImplementedError(f'No step location pattern for {settings.translate.DEST}')
    
    # select the correct step call from the list of step calls
    # need to standardise the step name and call name to properly compare. 
    step_symbol = _standardise_symbol(step_name)
    best_score, best_loc = 999, None
    for loc, call_symbol in matches:
        call_symbol = _standardise_symbol(call_symbol)
        score = levenshtein_distance(step_symbol, call_symbol)
        if score < best_score:
            best_score = score
            best_loc = loc
    return best_loc

def _standardise_symbol(symbol: str) -> str:
    return symbol.replace('_', '').replace('-', '').upper()

def _get_step_call_starts_cwl(translated: str) -> list[Tuple[int, str]]:
    # gather things which look like entity names
    # [ \t]+ included to get the location for the start of the line.
    PATTERN_FMT1 = r'(?<!(in:|out:|run:|requirements:|hints:|label:|doc:|scatter:|scatterMethod:)[\s]+)[ \t]*- ?id: ?([\w]+)[ \t]*\n'
    PATTERN_FMT2 = r'[ \t]+-? ?(?!in:|out:|run:|requirements:|hints:|label:|doc:|scatter:|scatterMethod:)([\w]+:[ \t]*\n)'
    PATTERN = f'{PATTERN_FMT1}|{PATTERN_FMT2}'
    step_call_matches = re.finditer(PATTERN, translated)
    
    # find where the steps section is
    PATTERN = r'steps:[\s\S]*?(?=outputs:|inputs:|$)'
    step_sec_match = list(re.finditer(PATTERN, translated))[0]
    steps_start, steps_end = step_sec_match.start(), step_sec_match.end()
    
    # filter to entity name matches in the steps section
    out = []
    for m in step_call_matches:
        if m.start() > steps_start and m.end() < steps_end:
            if m.group(2) is not None:
                name = m.group(2)
            elif m.group(3) is not None:
                name = m.group(3)
            else:
                raise RuntimeError('No group found') 
            # line start and correct group for each match
            out.append((m.start(), name))
    return out

def _get_step_call_starts_nextflow(translated: str) -> list[Tuple[int, str]]:
    # this is weak compared to CWL 
    # [ \t]+ included to get the location for the start of the line. 
    PATTERN = r'[ \t]+([A-Z_]+)\('
    matches = re.finditer(PATTERN, translated)
    
    # return the line start and correct group for each match
    return [(m.start(), m.group(1)) for m in matches]

def _get_step_call_starts_wdl(translated: str) -> list[Tuple[int, str]]:
    # [ \t]+ included to get the location for the start of the line.
    PATTERN = r'[ \t]+call ([\w.]+)( as ([\w.]+))? ?{'
    matches = re.finditer(PATTERN, translated)
    
    # get the line start and correct group for each match
    out = []
    for m in matches:
        if len(m.groups()) == 3:
            # if there is an 'as' clause, use what follows the 'as' clause
            name = m.group(3)
            out.append((m.start(), name))
        else:
            name = m.group(1)
            if '.' in name:
                # handling "call reports.align_and_count_summary {" case
                name = name.split('.')[-1]
            out.append((m.start(), name))
    return out

def _gen_block(
    internal: WorkflowBuilder| CommandToolBuilder | CodeTool, 
    translated: str,
    heading: Optional[str]=None,
    subsections: Optional[list[Any]]=None
    ) -> Optional[str]:

    # TODO HERE USE THESE FORMATS INSTEAD ----
    # Tool / Workflow default:
        
        ### MESSAGES ###

        # [WARNING][DATATYPES] test message1
        # [WARNING][METADATA] test message2
        # [WARNING][EXPERIMENTAL] test message4
        # [ERROR][FALLBACKS] test message3

        # UNTRANSLATED EXPRESSIONS
        # __TOKEN1__ = 'hello'
    
    # Workflow split:

        ### MESSAGES ###

        # GENERAL
        # [WARNING][DATATYPES] test message1

        # INPUTS
        # [WARNING][METADATA] test message2

        # OUTPUTS
        # [ERROR][FALLBACKS] test message3

        # UNTRANSLATED EXPRESSIONS
        # __TOKEN1__ = 'hello'
        
    # text for this block
    components = _gen_block_components(internal, translated, subsections)
    if components:
        # add extra component (the heading) for this block when using subsections
        commenter = commenter_map[settings.translate.DEST]
        if heading is not None: 
            if commenter == '#':
                heading = f'### {heading.upper()} ###'
            else:
                heading = f'{commenter} ### {heading.upper()} ###'
            components = [heading] + components
        return '\n\n'.join(components)
    return None

def _gen_block_components(    
    internal: WorkflowBuilder| CommandToolBuilder | CodeTool, 
    translated: str,
    subsections: Optional[list[Any]]=None
    ) -> list[str]:
    
    components: list[str] = []
    commenter = commenter_map[settings.translate.DEST]
    
    for cat in ErrorCategory:
        # we do this later - its weird. 
        if cat == ErrorCategory.SCRIPTING:
            continue
        loglines = load_loglines(category=cat, tool_uuid=internal.uuid, subsections=subsections)
        if loglines:
            lines = [f'{commenter} [{cat.value[0]}][{cat.value[1]}] {x.message}\n' for x in loglines]
            components.append('\n'.join(lines))

    cat = ErrorCategory.SCRIPTING
    loglines = load_loglines(category=cat, tool_uuid=internal.uuid, subsections=subsections)
    if loglines:
        messages = [x.message for x in loglines]
        
        # filter to messages where __TOKEN__ still is in the process text
        filtered = []
        for msg in messages:
            token = msg.split(' = ', 1)[0]
            if token in translated:
                filtered.append(msg) 
        if filtered:
            header = f'{commenter} UNTRANSLATED EXPRESSIONS'
            # have to add comment for cwl otherwise invalid syntax. 
            if settings.translate.DEST == 'cwl':
                filtered = [f'{commenter} {msg}' for msg in filtered]
            component = '\n'.join([header] + filtered)
            components.append(component)

    return components 


