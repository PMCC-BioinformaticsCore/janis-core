

import regex as re
from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional

from ...text.simplification.simplify import simplify_cmd
from ...model import XMLParam

MULTILINE_BOOL_MATCHER = r'#if.*?__SEARCH_TERM__[\s\S]*?#end if'
LINUX_CMD_MATCHER = r'\s*?(set |ln |cp |mv |export |mkdir |tar |ls |head |wget |grep |awk |cut |sed |gzip |gunzip |cd |echo |trap |touch ).*?$'
CHEETAH_MACRO_MATCHER = r'\s*?(#set |#import |#from |#silent |#echo ).*?$'
CHEETAH_CONDITIONAL_MATCHER = r'\s*?(#if |#unless |#else if |#elif ).*?$'
CHEETAH_LOOP_MATCHER = r'\s*?(#for |#while ).*?$'


class CmdstrReferenceType(Enum):
    MULTILINE_BOOL              = auto()
    INLINE_LINUX_CMD            = auto()
    INLINE_CHEETAH_MACRO        = auto()
    INLINE_CHEETAH_CONDITIONAL  = auto()
    INLINE_CHEETAH_LOOP         = auto()
    INLINE_PLAIN_TEXT           = auto()

@dataclass
class CmdstrReference:
    text: str
    rtype: CmdstrReferenceType

        
def get_cmdstr_appearences(
    cmdstr: str, 
    entity: str | XMLParam, 
    filter_to: Optional[CmdstrReferenceType | list[CmdstrReferenceType]]=None
    ) -> list[CmdstrReference]:
    """find all appearences of a param or string in the cmdstr, according to the CmdstrReferenceTypes we can look for. """
    
    # set up base search term
    if isinstance(entity, XMLParam):
        pname = entity.name.replace(r'.', r'\.')
        search_term = fr'\${pname}'
    else:
        search_term = entity
    
    # find appearences in cmdstr
    appearences = _get_all_appearences(cmdstr, search_term)
    
    # filter if required
    if isinstance(filter_to, CmdstrReferenceType):
        filter_to = [filter_to]
    if filter_to:
        appearences = [x for x in appearences if x.rtype in filter_to]
    return appearences
    
def _get_all_appearences(cmdstr: str, search_term: str) -> list[CmdstrReference]:
    cmdstr = simplify_cmd(cmdstr, purpose='parsing')
    refs: list[CmdstrReference] = []
    refs += _get_multiline_appearences(cmdstr, search_term)
    refs += _get_inline_appearences(cmdstr, search_term)
    return refs
    
def _get_multiline_appearences(cmdstr: str, search_term: str) -> list[CmdstrReference]:
    """only checking for multiline bool pattern at this stage. may add others.""" 
    return _get_multiline_bool_pattern_refs(cmdstr, search_term)

def _get_multiline_bool_pattern_refs(cmdstr: str, search_term: str) -> list[CmdstrReference]:
    """ 
    checking for this situation:
    #if $out.filtCounts:
        -F
    #end if
    """
    pattern = MULTILINE_BOOL_MATCHER.replace('__SEARCH_TERM__', search_term)
    iterator = re.finditer(pattern, cmdstr, re.MULTILINE)
    matches = [m for m in iterator]
    
    reflist: list[CmdstrReference] = []
    # checking each match has 3 lines to avoid nested control structure matches
    for match in matches:
        text = match.group(0)
        textlines = text.split('\n')
        if len(textlines) == 3:
            ref = CmdstrReference(text, CmdstrReferenceType.MULTILINE_BOOL)
            reflist.append(ref)
    return reflist

def _get_inline_appearences(cmdstr: str, search_term: str) -> list[CmdstrReference]:
    # all follow same pattern
    pattern = rf'.*?{search_term}.*?(?=\n)'
    iterator = re.finditer(pattern, cmdstr)
    matches = [m for m in iterator]
    
    reflist: list[CmdstrReference] = []
    for match in matches:
        line = match.group(0)
        if re.search(LINUX_CMD_MATCHER, line):
            ref = CmdstrReference(line, CmdstrReferenceType.INLINE_LINUX_CMD)
        elif re.search(CHEETAH_MACRO_MATCHER, line):
            ref = CmdstrReference(line, CmdstrReferenceType.INLINE_CHEETAH_MACRO)
        elif re.search(CHEETAH_CONDITIONAL_MATCHER, line):
            ref = CmdstrReference(line, CmdstrReferenceType.INLINE_CHEETAH_CONDITIONAL)
        elif re.search(CHEETAH_LOOP_MATCHER, line):
            ref = CmdstrReference(line, CmdstrReferenceType.INLINE_CHEETAH_LOOP)
        else:
            ref = CmdstrReference(line, CmdstrReferenceType.INLINE_PLAIN_TEXT)
        reflist.append(ref)
    
    return reflist

