
from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from janis_core import CommandToolBuilder, CodeTool, WorkflowBuilder

from Levenshtein import distance as levenshtein_distance
import regex as re
from typing import Optional, Tuple
from abc import ABC, abstractmethod

from janis_core import settings
from janis_core.tool.tool import ToolType
from .logfile import ErrorCategory
from .main import load_loglines
from .gather import gather_uuids
from .enums import FormatCategory

### MISC CONSTANTS ###

COMMENTER_MAP = {
    'nextflow': '//',
    'cwl': '#',
    'wdl': '#',
}

### PUBLIC ###

def inject_messages(internal: WorkflowBuilder | CommandToolBuilder | CodeTool, translated: str) -> str:
    """
    Injects messages into the translated file text.
    
    Note: 
    - Messages need to be logged during ingest to use this feature.
    - When logging messages during ingest you need to provide the unique entity_uuid. 
    - This enables the right messages to be written to the right output files. 
    - Example:
        log_message(tool.uuid, msg=="couldn't parse datatype...", category=ErrorCategory.DATATYPE)
    """
    if internal.type() in [ToolType.CommandTool, ToolType.CodeTool]:
        return ToolInjector(internal, translated).inject()          # type: ignore
    elif internal.type() == ToolType.Workflow and settings.messages.USE_SUBSECTIONS:
        return WorkflowInjectorSplit(internal, translated).inject() # type: ignore
    elif internal.type() == ToolType.Workflow:
        return WorkflowInjector(internal, translated).inject()      # type: ignore
    else:
        raise RuntimeError(f"Can't inject messages into {internal.type()}")


### PRIVATE ###

class MessageInjector(ABC):

    @abstractmethod
    def inject(self) -> str:
        ...

    def _gen_block(self, heading: Optional[str]=None, uuids: Optional[set[str]]=None, is_scripting: bool=False) -> Optional[str]:
        # text for this block
        msgs = self._gen_block_lines(uuids, is_scripting)
        if msgs:
            # add extra component (the heading) for this block when using subsections
            if heading is not None: 
                msgs = [heading.upper()] + msgs
            
            # comment each line & return as string
            commenter = COMMENTER_MAP[settings.translate.DEST]
            msgs = [f'{commenter} {ln}' for ln in msgs]
            return '\n'.join(msgs)
        
        return None

    def _gen_block_lines(self, uuids: Optional[set[str]]=None, is_scripting: bool=False) -> list[str]:
        msgs: list[str] = []

        # scripting block
        if is_scripting:
            cat = ErrorCategory.SCRIPTING
            loglines = load_loglines(category=cat, entity_uuids=uuids)
            if loglines:
                messages = [x.message for x in loglines]
                
                # filter to messages where __TOKEN__ still is in the process text
                filtered = []
                for msg in messages:
                    token = msg.split(' = ', 1)[0]
                    if token in self.translated: # type: ignore
                        filtered.append(msg) 
                if filtered:
                    # have to add comment for cwl otherwise invalid syntax. 
                    if settings.translate.DEST == 'cwl':
                        filtered = [f'{msg}' for msg in filtered]
                    msgs += filtered
            return msgs

        # normal block
        for cat in ErrorCategory:
            # we do this later - its weird. 
            if cat == ErrorCategory.SCRIPTING:
                continue
            loglines = load_loglines(category=cat, entity_uuids=uuids)
            if loglines:
                msgs += [f'[{cat.value[0]}][{cat.value[1]}] {x.message}' for x in loglines]
        return msgs 



class ToolInjector(MessageInjector):
    
    def __init__(self, internal: CommandToolBuilder | CodeTool, translated: str) -> None:
        self.internal = internal
        self.translated = translated
        self.entity_uuids_map = gather_uuids(internal)

    def inject(self) -> str:
        uuids = set(self.entity_uuids_map.keys())
        
        # scripting messages
        messages = self._gen_block(uuids=uuids, is_scripting=True)
        if messages:
            commenter = COMMENTER_MAP[settings.translate.DEST]
            self.translated = f'{commenter} {settings.messages.SCRIPTING_BANNER}\n{messages}\n\n{self.translated}'

        # normal messages
        messages = self._gen_block(uuids=uuids)
        if messages: 
            commenter = COMMENTER_MAP[settings.translate.DEST]
            self.translated = f'{commenter} {settings.messages.MESSAGES_BANNER}\n{messages}\n\n{self.translated}'
        return self.translated


class WorkflowInjector(MessageInjector):
    """
    Dumps messages to top of file, except for step specific messages which are dumped 
    above the relevant step call.
    """
    
    def __init__(self, internal: WorkflowBuilder, translated: str) -> None:
        self.internal = internal
        self.translated = translated
        self.entity_uuids_map = gather_uuids(internal)
    
    def inject(self) -> str:
        self.inject_header_messages()
        self.inject_step_messages()
        return self.translated

    def inject_header_messages(self) -> None:
        uuids = set([uuid for uuid, fcat in self.entity_uuids_map.items() if fcat != FormatCategory.STEP])
        
        # scripting messages
        messages = self._gen_block(uuids=uuids, is_scripting=True)
        if messages:
            commenter = COMMENTER_MAP[settings.translate.DEST]
            self.translated = f'{commenter} {settings.messages.SCRIPTING_BANNER}\n{messages}\n\n{self.translated}'
        
        # normal messages
        messages = self._gen_block(uuids=uuids)
        if messages:
            commenter = COMMENTER_MAP[settings.translate.DEST]
            self.translated = f'{commenter} {settings.messages.MESSAGES_BANNER}\n{messages}\n\n{self.translated}'
        
    def inject_step_messages(self) -> None:
        """
        Note: No headings or subheadings for step-related messages, except scripting.
        For each step: 
        - get messages
        - if messages, find step call
        - if found, get indent
        - inject message block above step call with correct indent
        """
        for sname, step in self.internal.step_nodes.items():
            step_uuids_map = gather_uuids(step)
            
            # sanity check
            assert all([fcat == FormatCategory.STEP for fcat in step_uuids_map.values()])

            # get message block for this step if step has messages
            uuids = set(step_uuids_map.keys())
            normal_messages = self._gen_block(uuids=uuids)
            scripting_messages = self._gen_block(uuids=uuids, is_scripting=True)
            if not normal_messages and not scripting_messages:
                continue

            # get step location in translated text
            loc = _get_step_loc(sname, self.translated)
            if loc is None:
                print(self.translated)
                raise NotImplementedError
            
            # get indent for the step call - awful in general but works
            # split file into top and bottom pivoting on the match location
            top = self.translated[:loc]
            bottom = self.translated[loc:]

            # get the indent level of the step call
            if not bottom.startswith(' ') and not bottom.startswith('\t'):
                indent = ''
            else:
                indent = re.findall(r'[ \t]+', bottom)[0]

            # apply indent to each line of the message block
            if normal_messages:
                msg_lines = normal_messages.split('\n')
                msg_lines = [indent + ln for ln in msg_lines]
                normal_messages = '\n'.join(msg_lines)
            
            if scripting_messages:
                # add special header
                msg_lines = ['UNTRANSLATED EXPRESSIONS']
                msg_lines += scripting_messages.split('\n')
                msg_lines = [indent + ln for ln in msg_lines]
                scripting_messages = '\n'.join(msg_lines)

            messages = ''
            if normal_messages:
                messages += normal_messages
            if scripting_messages:
                messages += f'\n\n{scripting_messages}'
            
            # inject indented message block above step call
            self.translated = top + messages + '\n' + bottom
        


class WorkflowInjectorSplit(WorkflowInjector):
    """
    Instead of dumping all header messages together at top of file, organises header into 
    GENERAL, INPUTS, OUTPUTS sections, followed by UNTRANSLATED EXPRESSIONS
    """

    # single override method for this class
    def inject_header_messages(self) -> None:
        main_uuids = set([uuid for uuid, fcat in self.entity_uuids_map.items() if fcat == FormatCategory.MAIN])
        input_uuids = set([uuid for uuid, fcat in self.entity_uuids_map.items() if fcat == FormatCategory.INPUT])
        output_uuids = set([uuid for uuid, fcat in self.entity_uuids_map.items() if fcat == FormatCategory.OUTPUT])
        
        # get the messages for header subheadings as blocks
        main = self._gen_block(heading='general', uuids=main_uuids)
        inputs = self._gen_block(heading='inputs', uuids=input_uuids) 
        outputs = self._gen_block(heading='outputs', uuids=output_uuids)
        
        # get the messages for header scripting as block
        all_uuids = main_uuids | input_uuids | output_uuids
        scripting = self._gen_block(heading='outputs', uuids=all_uuids, is_scripting=True)
        
        # write scripting
        if scripting:
            commenter = COMMENTER_MAP[settings.translate.DEST]
            self.translated = f'{commenter} {settings.messages.SCRIPTING_BANNER}\n{scripting}\n\n{self.translated}'
        
        # write normal messages
        messages = '\n\n'.join([x for x in [main, inputs, outputs] if x is not None])
        if messages:
            commenter = COMMENTER_MAP[settings.translate.DEST]
            self.translated = f'{commenter} {settings.messages.MESSAGES_BANNER}\n{messages}\n\n{self.translated}'


def _get_step_loc(step_name: str, translated: str) -> Optional[int]:
    """
    returned loc should be the first character of the target line, not the start of the step symbol
    eg nextflow process call:
        "               MINIMAP2("
         ^ loc is here, ^ not here 
    """ 

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
    # AN: this is ew
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
                # TODO downgrade this to a warning? or just ignore?
                raise RuntimeError('No group found') 
            # line start and correct group for each match
            out.append((m.start(), name))
    return out

def _get_step_call_starts_nextflow(translated: str) -> list[Tuple[int, str]]:
    # this is weak compared to CWL 
    # [ \t]+ included to get the location for the start of the line. 
    PATTERN = r'[ \t]+([A-Z0-9_]+)\('
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




