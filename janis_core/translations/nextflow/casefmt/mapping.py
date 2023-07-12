

def as_camel(words: list[str]) -> str:
    if len(words) == 0:
        return ''
    elif len(words) == 1:
        return words[0]
    else:
        return words[0] + ''.join([x.capitalize() for x in words[1:]])

def as_pascal(words: list[str]) -> str:
    if len(words) == 0:
        return ''
    else:
        return ''.join([x.capitalize() for x in words])

def as_snake(words: list[str]) -> str:
    return '_'.join(words)

def as_snake_caps(words: list[str]) -> str:
    words = [x.upper() for x in words]
    return '_'.join(words)