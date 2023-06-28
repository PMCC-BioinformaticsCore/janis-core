

import shlex


def split_lines(cmdstr: str) -> list[str]:
    sh = shlex.shlex(cmdstr)
    sh.commenters = ''
    sh.whitespace = '\n'
    sh.whitespace_split = True
    lines = list(sh)
    lines = [ln.strip() for ln in lines]
    lines = [ln for ln in lines if ln != '']
    return lines 

def join_lines(lines: list[str]) -> str:
    return '\n'.join(lines)

def split_to_words(line: str) -> list[str]:
    return shlex.split(line, posix=False)

