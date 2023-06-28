



import regex as re
from .patterns import *

def convert(text: str) -> str:
    text = _do_wildcard_replacements(text)
    text = _do_logical_or_replacements(text)
    text = _do_char_set_replacements(text)
    text = _do_special_char_replacements(text)
    text = _do_escape_replacements(text)
    return text

def _do_wildcard_replacements(text: str) -> str:
    patterns = [
        ZERO_OR_MORE,
        ONE_OR_MORE
    ]
    for pattern in patterns:
        while re.findall(pattern, text):
            match = list(re.finditer(pattern, text))[0]
            text = text[:match.start()] + "*" + text[match.end():]
    return text

def _do_logical_or_replacements(text: str) -> str:
    while re.findall(LOGICAL_OR, text):
        match = list(re.finditer(LOGICAL_OR, text))[0]
        inner = match.group(1).replace('|', ',')
        text = text[:match.start()] + f"{{{inner}}}" + text[match.end():]
    return text

def _do_char_set_replacements(text: str) -> str:
    if not re.findall(CHAR_SET, text):
        return text
    
    new_text = ''
    prev_match = None
    for idx, match in enumerate(re.finditer(CHAR_SET, text)):
        # for first match, add everything before the match
        if idx == 0:
            new_text += text[:match.start()]
        else:
            new_text += text[prev_match.end():match.start()]

        old_inner = match.group(1)
        new_inner_list: list[str] = []
        
        # go through each char, group by ranges or single chars & 
        i = 0
        while i < len(old_inner):
            
            # negative set
            if i == 0 and old_inner[0] == '^':
                new_inner_list.append('!')

            # looking for ranges
            elif i + 2 < len(old_inner):
                if old_inner[i].isalnum() and old_inner[i+1] == '-' and old_inner[i+2].isalnum():
                    new_inner_list.append(old_inner[i:i+3])
                    i += 2
                else:
                    new_inner_list.append(old_inner[i])
            
            else:
                new_inner_list.append(old_inner[i])

            i += 1
        
        # join the negative set marker with the first actual item
        if new_inner_list[0] == '!':
            new_inner_list = [new_inner_list[0] + new_inner_list[1]] + new_inner_list[2:] 
        new_inner = ','.join(new_inner_list)
        new_text += f"[{new_inner}]"

        prev_match = match

    # final match need to add end also 
    new_text += text[match.end():]

    return new_text

def _do_special_char_replacements(text: str) -> str:
    if not re.findall(SPECIAL_CHARS, text):
        return text

    new_text = ''
    prev_match = None
    for idx, match in enumerate(re.finditer(SPECIAL_CHARS, text)):
        # for first match, add everything before the match
        if idx == 0:
            new_text += text[:match.start()]
        else:
            new_text += text[prev_match.end():match.start()]
        
        if not inside_set(match, text):
            new_text += '?'
        else:
            new_text += match.group(0)
        
        prev_match = match
    
    # final match need to add end also 
    new_text += text[prev_match.end():]
    return new_text

def inside_set(match: re.Match, text: str) -> bool:
    bracket_locations = [(m.start(), m.end()) for m in re.finditer(BRACKETS, text)]
    for start, end in bracket_locations:
        if match.start() or match.end() in range(start, end):
            return True
    return False

def _do_escape_replacements(text: str) -> str:
    return text.replace('\\', '')