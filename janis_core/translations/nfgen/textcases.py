
from janis_core.translations.nfgen import patterns

from enum import Enum, auto
import regex as re 

# module entry
def to_case(text: str, format: str='snake_caps') -> str:
    """
    casts text to specific case format.
    format is either:\n
    - camel (camelCase)\n 
    - snake (snake_case)\n
    - snake_caps (SNAKE_CAPS)
    """
    src_case = identify_case(text)
    text_split = split_words(text, src_case)
    dest_case = case_map[format]
    return as_case(text_split, dest_case)




class Case(Enum):
    camel       = auto()
    snake       = auto()
    snake_caps  = auto()


def as_case(words: list[str], case: Case) -> str:
    func = as_map[case]
    return func(words)

def as_camel(words: list[str]) -> str:
    if len(words) == 0:
        return ''
    elif len(words) == 1:
        return words[0]
    else:
        return words[0] + ''.join([x.capitalize() for x in words[1:]])

def as_snake(words: list[str]) -> str:
    return '_'.join(words)

def as_snake_caps(words: list[str]) -> str:
    words = [x.upper() for x in words]
    return '_'.join(words)


case_map = {
    'camel': Case.camel,
    'snake': Case.snake,
    'snake_caps': Case.snake_caps
}

regex_map = {
    Case.camel: patterns.CAMEL,
    Case.snake: patterns.SNAKE,
    Case.snake_caps: patterns.SNAKE_CAPS
}

as_map = {
    Case.camel: as_camel,
    Case.snake: as_snake,
    Case.snake_caps: as_snake_caps
}

def identify_case(text: str) -> Case:
    for case, pattern in regex_map.items():
        matches = re.finditer(pattern, text)
        for m in matches:
            if m.start() == 0 and m.end() == len(text):
                return case
    
    # haven't found good evidence of case yet
    if '_' in text:
        return Case.snake 
    else:
        return Case.camel 

def split_words(text: str, case: Case) -> list[str]:
    out: list[str] = []
    matches = re.finditer(patterns.CAMEL, text)
    if case == Case.camel:
        for m in matches:
            for i, group in enumerate(m.groups()):
                out += m.captures(i + 1)
    elif case in [Case.snake, Case.snake_caps]:
        out += text.split('_')
    else:
        raise NotImplementedError
    return [x.lower() for x in out]





# tests 

# my_text1 = 'helloThere'
# my_text2 = 'helloThereFriend2'
# my_text3 = 'hello_there_friend2'
# my_text4 = 'hello_there_Friend2'
# my_text5 = 'Hello_There_Friend2'

# print()
# print(my_text1)
# print(to_case(my_text1, format='snake_caps'))

# print()
# print(my_text2)
# print(to_case(my_text2, format='snake_caps'))

# print()
# print(my_text3)
# print(to_case(my_text3, format='snake_caps'))

# print()
# print(my_text4)
# print(to_case(my_text4, format='snake_caps'))

# print()
# print(my_text5)
# print(to_case(my_text5, format='snake_caps'))

