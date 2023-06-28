


def split_lines(cmdstr: str) -> list[str]:
    lines = cmdstr.split('\n')
    lines = [ln.strip() for ln in lines]
    lines = [ln for ln in lines if ln != '']
    return lines 

def split_lines_blanklines(cmdstr: str) -> list[str]:
    lines = cmdstr.split('\n')
    lines = [ln.strip() for ln in lines]
    return lines 

def split_lines_whitespace(cmdstr: str) -> list[str]:
    return cmdstr.split('\n')

def join_lines(cmdlines: list[str]) -> str:
    return '\n'.join(cmdlines)