
from . import patterns
from . import mapping
import regex as re 


# module entry
def to_case(text: str, case: str='snake') -> str:
    """
    casts text to specific case format.
    format is either:\n
    - camel (camelCase)\n 
    - pascal (PascalCase)\n 
    - snake (snake_case)\n
    - snake_caps (SNAKE_CAPS)
    """
    if text == 'NamingTestWF':
        print()
    text_split = split_words(text)
    mapping_func = mapping_func_map[case]
    final_text = mapping_func(text_split)
    return final_text

def split_words(text: str) -> list[str]:
    out: list[str] = [text]
    out = split_underscores(out)
    out = split_camel(out)
    out = split_pascal(out)
    out = [x.lower() for x in out]
    return out

def split_underscores(words: list[str]) -> list[str]:
    out: list[str] = []
    for word in words:
        split_word = word.split('_')
        out += split_word
    return out

def split_camel(words: list[str]) -> list[str]:
    out: list[str] = []
    for word in words:
        split = False
        matches = re.finditer(patterns.CAMEL, word)
        for m in matches:
            if m.start() == 0 and m.end() == len(word):
                split = True
                for i, group in enumerate(m.groups()):
                    out += m.captures(i + 1)
        if not split:
            out.append(word)
    out = [x for x in out if not x == '']
    return out

def split_pascal(words: list[str]) -> list[str]:
    out: list[str] = []
    for word in words:
        split = False
        matches = re.finditer(patterns.PASCAL, word)
        for m in matches:
            if m.start() == 0 and m.end() == len(word):
                split = True
                for i, group in enumerate(m.groups()):
                    out += m.captures(i + 1)
        if not split:
            out.append(word)
    out = [x for x in out if not x == '']
    return out


mapping_func_map = {
    'camel': mapping.as_camel,
    'pascal': mapping.as_pascal,
    'snake': mapping.as_snake,
    'snake_caps': mapping.as_snake_caps,
}

# tests 

# my_text1 = 'helloThere'
# my_text2 = 'helloThereB'
# my_text3 = 'helloThereFriend2'
# my_text4 = 'hello_there_friend2'
# my_text5 = 'hello_there_Friend2'
# my_text6 = 'Hello_There_Friend2'
# my_text7 = 'HelloThere_Friend2'
# my_text8 = 'HelloThere_friend2'
# my_text9 = 'HelloThere_friendB'

# print()
# print(my_text1)
# print(to_case(my_text1, case='snake_caps'))

# print()
# print(my_text2)
# print(to_case(my_text2, case='snake_caps'))

# print()
# print(my_text3)
# print(to_case(my_text3, case='snake_caps'))

# print()
# print(my_text4)
# print(to_case(my_text4, case='snake_caps'))

# print()
# print(my_text5)
# print(to_case(my_text5, case='snake_caps'))

# print()
# print(my_text6)
# print(to_case(my_text6, case='snake_caps'))

# print()
# print(my_text7)
# print(to_case(my_text7, case='snake_caps'))

# print()
# print(my_text8)
# print(to_case(my_text8, case='snake_caps'))

# print()
# print(my_text9)
# print(to_case(my_text9, case='snake_caps'))

