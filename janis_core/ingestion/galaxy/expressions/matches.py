

from typing import Optional, Tuple
import regex as re
import numpy as np

from .patterns import QUOTED_SECTION 


def get_matches(the_string: str, expression: str) -> list[re.Match[str]]:
    matches = re.finditer(expression, the_string)
    return [m for m in matches]

def get_next_word(word: str, delim: str, text: str) -> Optional[str]:
    NEXT_WORD = r'(?<=(?:\s|^))' + f'{word}{delim}' + r'+?([\w\d\'"${}\\_.\-\:/]+)(?=\s|$)'
    matches = re.finditer(NEXT_WORD, text)
    matches = [m for m in matches]
    if matches:
        value = matches[0].group(1)
        value = value.strip('"\'')
        return value
    return None

def get_preceeding_dashes(search_term: str, text: str) -> list[str]:
    PRECEEDING_DASHES = r'(?<![$.{])(-+?)' + fr'({search_term})' + r'(?=[\s=:]|$|[\'"])'
    matches = re.finditer(PRECEEDING_DASHES, text)
    return [m.group(1) for m in matches]

def get_quoted_sections(the_string: str):
    # find the areas of the string which are quoted
    matches = re.finditer(QUOTED_SECTION, the_string)
    quoted_sections = [(m.start(), m.end()) for m in matches]

    # transform to mask
    quotes_mask = np.zeros(len(the_string)) # type: ignore
    for start, end in quoted_sections:
        quotes_mask[start: end] = 1
    
    return quotes_mask

def find_unquoted(the_string: str, pattern: str) -> Tuple[int, int]:
    """
    finds the pattern in string. ensures section is not quoted. 
    """
    # find quoted sections of input string
    quotes_mask = get_quoted_sections(the_string)
    
    # get pattern match locations
    matches = re.finditer(pattern, the_string)
    match_spans = [(m.start(), m.end()) for m in matches]

    # check each match to see if its in a quoted section
    if len(match_spans) > 0:
        for start, end in match_spans:
            if sum(quotes_mask[start: end]) == 0:
                # return position of first unquoted match
                return start, end
    return -1, -1

# def get_unpaired_quotes_start(the_string: str) -> int:
#     # find all quotes
#     matches = re.finditer(QUOTES, the_string)
#     quote_locations = [m.start() for m in matches]
#     # find paired quotes sections
#     quotes_mask = get_quoted_sections(the_string)
#     # check each quote is part of quoted section
#     for loc in quote_locations:
#         if quotes_mask[loc] != 1:
#             return loc
#     return -1