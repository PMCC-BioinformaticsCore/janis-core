

def split_lines_blanklines(text: str) -> list[str]:
    lines = text.split('\n')
    lines = [ln.strip() for ln in lines]
    return lines 

def split_lines(text: str) -> list[str]:
    return text.split('\n')

def join_lines(lines: list[str]) -> str:
    return '\n'.join(lines)