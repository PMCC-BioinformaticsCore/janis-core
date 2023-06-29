

import regex as re
from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional, Tuple
from janis_core.ingestion.galaxy.gxtool.model import XMLTool

from ...text.simplification.simplify import simplify_cmd
from ...model import XMLParam, XMLBoolParam, XMLSelectParam
from ..tokens import Token, TokenType


BLANK = r'^\s*$'
SIMPLE_VARIABLE = r'^([\'"\t ]?\$\w[\w\d._]+[\'"\t ]?)$'
SIMPLE_PREFIX = r'^-[\w\d_-]+$'
SIMPLE_OPTION = r'^(-[\w\d_-]+)( +|=|:)([\'"]?\$\w[\w\d._]+[\'"]?)$'

INLINE_PARAM_MATCHER = r'(^|.*?[\'"(\t ]+)\$(__PARAM_NAME__)(([\'"):\t ]+.*?(?=\n))|([\'"):\t ]*?(?=\n))|$)'
MULTILINE_BOOL_MATCHER = r'#if.*?\$(__PARAM_NAME__)[\s\S]*?#end if'
LINUX_CMD_MATCHER = r'(^|.*? )(set |ln |cp |mv |export |mkdir |tar |ls |head |wget |grep |awk |cut |sed |gzip |gunzip |cd |echo |trap |touch ).*?$'
CHEETAH_MACRO_MATCHER = r'(^|.*? )(#set |#import |#from |#silent |#echo ).*?$'
CHEETAH_CONDITIONAL_MATCHER = r'(^|.*? )(#if |#unless |#else if |#elif ).*?$'
CHEETAH_LOOP_MATCHER = r'(^|.*? )(#for |#while ).*?$'

class CmdstrReferenceType(Enum):
    MULTILINE_BOOL              = auto()
    INLINE_LINUX_CMD            = auto()
    INLINE_CHEETAH_MACRO        = auto()
    INLINE_CHEETAH_CONDITIONAL  = auto()
    INLINE_CHEETAH_LOOP         = auto()
    INLINE_PLAIN_TEXT           = auto()

@dataclass
class CmdstrReference:
    cmdstr: str
    start: int
    text: str
    rtype: CmdstrReferenceType



### PUBLIC FUNCTIONS ###

def is_bool_select(token: Token) -> bool:
    if token.ttype == TokenType.GX_INPUT:
        if isinstance(token.gxparam, XMLBoolParam | XMLSelectParam):
            return True
    return False

def is_blank(phrase: str) -> bool:
    if re.search(BLANK, phrase):
        return True
    return False

def is_prefix(phrase: str) -> bool:
    phrase = phrase.strip()
    if re.search(SIMPLE_PREFIX, phrase):
        return True
    return False

def is_simple_flags(phrase: str) -> bool:
    words = phrase.strip().split()
    for word in words:
        if not re.search(SIMPLE_PREFIX, word):
            return False
    return True

def is_simple_variable(phrase: str) -> bool:
    phrase = phrase.strip()
    if re.search(SIMPLE_VARIABLE, phrase):
        return True
    return False

def is_simple_phrases(phrase: str) -> bool:
    words = phrase.strip().split()
    for word in words:
        raise NotImplementedError
    return True

def is_simple_option(phrase: str) -> bool:
    phrase = phrase.strip()
    if re.search(SIMPLE_OPTION, phrase):
        return True
    return False

def extract_simple_option(phrase: str) -> Tuple[str, str, str]:
    phrase = phrase.strip()
    match = re.search(SIMPLE_OPTION, phrase)
    if not match:
        raise RuntimeError
    prefix = match.group(1)
    separator = match.group(2)
    value = match.group(3)
    return prefix, separator, value

def argument_resembles_prefix(argument: str, prefix: str) -> bool:
    if prefix in argument:
        return True
    # for when argument has incorrect number of dashes
    elif prefix.lstrip('-') == argument.lstrip('-'):
        return True
    return False

def single_inline_plaintext_appearence(xmltool: XMLTool, param: XMLParam) -> bool:
    appearences = get_cmdstr_appearences(xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)
    return len(appearences) == 1

def single_multiline_bool_appearence(xmltool: XMLTool, param: XMLParam) -> bool:
    appearences = get_cmdstr_appearences(xmltool.raw_command, param, filter_to=CmdstrReferenceType.MULTILINE_BOOL)
    return len(appearences) == 1

def get_cmdstr_appearences(
    cmdstr: str, 
    param: XMLParam, 
    filter_to: Optional[CmdstrReferenceType | list[CmdstrReferenceType]]=None
    ) -> list[CmdstrReference]:
    """find all appearences of a param or string in the cmdstr, according to the CmdstrReferenceTypes we can look for. """

    # find appearences in cmdstr
    appearences = _get_all_appearences(cmdstr, param)
    
    # filter if required
    if isinstance(filter_to, CmdstrReferenceType):
        filter_to = [filter_to]
    if filter_to:
        appearences = [x for x in appearences if x.rtype in filter_to]
    return appearences
    

### PRIVATE HELPERS ###

def _get_all_appearences(cmdstr: str, param: XMLParam) -> list[CmdstrReference]:
    cmdstr = simplify_cmd(cmdstr, purpose='parsing')
    refs: list[CmdstrReference] = []
    refs += _get_multiline_appearences(cmdstr, param)
    refs += _get_inline_appearences(cmdstr, param)
    return refs
    
def _get_multiline_appearences(cmdstr: str, param: XMLParam) -> list[CmdstrReference]:
    """only checking for multiline bool pattern at this stage. may add others.""" 
    return _get_multiline_bool_pattern_refs(cmdstr, param)

def _get_multiline_bool_pattern_refs(cmdstr: str, param: XMLParam) -> list[CmdstrReference]:
    """ 
    checking for this situation:
    #if $out.filtCounts:
        -F
    #end if
    """
    # set up base search term
    pname = param.name.replace(r'.', r'\.')
    pattern = MULTILINE_BOOL_MATCHER.replace('__PARAM_NAME__', pname)
    
    # find all matches
    iterator = re.finditer(pattern, cmdstr, re.MULTILINE)
    matches = [m for m in iterator]
    
    reflist: list[CmdstrReference] = []
    # checking each match has 3 lines to avoid nested control structure matches
    for match in matches:
        text = match.group(0)
        textlines = text.split('\n')
        if len(textlines) == 3:
            ref = CmdstrReference(cmdstr, match.start(), text, CmdstrReferenceType.MULTILINE_BOOL)
            reflist.append(ref)
    return reflist

def _get_inline_appearences(cmdstr: str, param: XMLParam) -> list[CmdstrReference]:
    # set up base search term
    pname = param.name.replace(r'.', r'\.')
    pattern = INLINE_PARAM_MATCHER.replace('__PARAM_NAME__', pname)

    # all follow same pattern
    iterator = re.finditer(pattern, cmdstr, re.MULTILINE)
    matches = [m for m in iterator]
    
    reflist: list[CmdstrReference] = []
    for match in matches:
        text = match.group(0)
        if re.search(LINUX_CMD_MATCHER, text):
            ref = CmdstrReference(cmdstr, match.start(), text, CmdstrReferenceType.INLINE_LINUX_CMD)
        elif re.search(CHEETAH_MACRO_MATCHER, text):
            ref = CmdstrReference(cmdstr, match.start(), text, CmdstrReferenceType.INLINE_CHEETAH_MACRO)
        elif re.search(CHEETAH_CONDITIONAL_MATCHER, text):
            ref = CmdstrReference(cmdstr, match.start(), text, CmdstrReferenceType.INLINE_CHEETAH_CONDITIONAL)
        elif re.search(CHEETAH_LOOP_MATCHER, text):
            ref = CmdstrReference(cmdstr, match.start(), text, CmdstrReferenceType.INLINE_CHEETAH_LOOP)
        else:
            ref = CmdstrReference(cmdstr, match.start(), text, CmdstrReferenceType.INLINE_PLAIN_TEXT)
        reflist.append(ref)
    
    return reflist

