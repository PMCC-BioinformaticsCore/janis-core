
import WDL
from janis_core import CommandToolBuilder

from .parsers import NativeCommandParser, ShellCommandParser, ShellScriptParser
from .cmdline import CmdLine
import regex as re


def parse_command(internal: CommandToolBuilder, task: WDL.Tree.Task) -> CommandToolBuilder:
    lines = split_newlines(task)
    cmds = split_cmdlines(lines)
    cmds = remove_comments(cmds)
    
    for parser_c in [NativeCommandParser, ShellCommandParser, ShellScriptParser]:
        try:
            parser = parser_c(internal, task, cmds)
            parser.parse()
            return parser.internal
        except Exception as e:
            pass

    return internal

def split_newlines(task: WDL.Tree.Task) -> list[str | WDL.Expr.Placeholder]:
    """does what"""
    text = []
    for part in task.command.parts:
        if isinstance(part, str):
            part = re.sub(r'\\[ \t]*?\n', '\\\n', part)
            lines = re.split(r'(?<!\\)\n', part)                    # split newlines
            lines = [ln.strip(' \t') for ln in lines]               # strip whitespace
            lines = [ln for ln in lines if not ln.startswith('#')]  # remove comment lines
            lines = [ln.split('#')[0] for ln in lines]              # strip comments
            text += lines
        else:
            text.append(part)
    return text

def split_cmdlines(text: list[str | WDL.Expr.Placeholder]) -> list[CmdLine]:
    cmds = []
    curr_cmd = CmdLine()
    
    for line in text:
        if isinstance(line, str):
            # newline: new command
            if line == '':
                if len(curr_cmd.elems) > 0:
                    cmds.append(curr_cmd)
                    curr_cmd = CmdLine()
            # escaped newline: do nothing
            elif line == '\\\n':
                continue
            # two string lines in a row means linebreak (new command)
            elif len(curr_cmd.elems) > 0 and isinstance(curr_cmd.elems[-1], str):
                cmds.append(curr_cmd)
                curr_cmd = CmdLine()
                line = line.replace('\\\n', '')
                line = line.strip(' \t')
                curr_cmd.elems.append(line)
                pass
            # add to current command
            else:
                line = line.replace('\\\n', '')
                line = line.strip(' \t')
                curr_cmd.elems.append(line)
        else:
            curr_cmd.elems.append(line)
    if len(curr_cmd.elems) > 0:
        cmds.append(curr_cmd)
    return cmds

def remove_comments(cmds: list[CmdLine]) -> list[CmdLine]:
    return [cmd for cmd in cmds if not cmd.is_comment]
    
    

