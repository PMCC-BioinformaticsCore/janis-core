
from typing import Optional 
from copy import deepcopy
import WDL
from janis_core import ToolInput
from janis_core import ToolOutput
from janis_core import ToolArgument
from janis_core import CommandToolBuilder
from dataclasses import dataclass, field

from .basic_parser import BasicCommandParser
from .shell_parser import ShellCommandParser
import regex as re


def parse_command(internal: CommandToolBuilder, task: WDL.Tree.Task) -> CommandToolBuilder:
    commands = get_text_commands(task)
    parser_c = ShellCommandParser if requires_shell_script(commands) else BasicCommandParser
    parser = parser_c(internal, task, commands)
    parser.parse()
    return parser.internal

@dataclass
class Command:
    elems: list[str | WDL.Expr.Placeholder] = field(default_factory=list)

    @property
    def is_shell_var_cmd(self) -> bool:
        if not isinstance(self.elems[0], str):
            return False
        if not re.match(r'^[ \t]*?(set|declare|env|export)', self.elems[0]):
            return False
        return True
    
    @property
    def is_comment(self) -> bool:
        if not isinstance(self.elems[0], str):
            return False
        if not re.match(r'^[ \t]*?#', self.elems[0]):
            return False
        return True

def get_text_commands(task: WDL.Tree.Task) -> list[Command]:
    cmds: list[Command] = []

    if len(task.command.parts) == 0:
        return cmds
    
    curr_cmd = Command()
    for i, part in enumerate(task.command.parts):
        
        if isinstance(part, str):
            lines = re.split(r'(?<!\\)\n', part)
            lines = [ln.strip(' \t') for ln in lines]
            if lines[0] != '' and isinstance(curr_cmd.elems[-1], WDL.Expr.Placeholder):
                



            if lines[0] != '':
                curr_cmd.elems.append(lines[0])


            for j, line in enumerate(lines):
                # first line: append to previous command
                if j == 0 and line != '' and isinstance(curr_cmd.elems[-1], WDL.Expr.Placeholder):
                    curr_cmd.elems.append(line)
                    cmds.append(curr_cmd)
                    curr_cmd = Command()
                        

                # middle lines

                # last line





                # ignore empty lines
                elif line == '':
                    continue 
                # ignore comments 
                elif re.match(r'^[ \t]*?#', line):
                    continue
                elif i != len(lines) - 1:
                    curr_cmd.elems.append(line)
                    cmds.append(curr_cmd)
                    curr_cmd = Command()
                elif 


                
                    

            lines = [ln.strip() for ln in lines]
            should_append = False if lines[-1] == '' else True
            lines = [ln for ln in lines if ln != '']

            print()
    print()


def requires_shell_script(commands: list[str]) -> bool:
    """Returns True if the task requires a shell script"""
    return True




class WDLCommandParser:
    def __init__(self, internal: CommandToolBuilder, task: WDL.Tree.Task):
        self.internal = internal
        self.task = task
        # self.base_command: list[str] = []
        # self.inputs: list[ToolInput] = []
        # self.arguments: list[ToolArgument] = []
        # self.outputs: list[ToolOutput] = []

    @property
    def requires_shell_script(self) -> bool:
        
        return True
        return False
    
    def parse(self) -> None:
        self.update_base_command()
        self.update_inputs()
        self.update_arguments()
        self.update_outputs()

    def update_base_command(self) -> None:
        text = self.extract_base_command_text()
        text = self.strip_newline_ignores(text)
        base_command = self.split_base_command(text)
        base_command = self.strip_set_e(base_command)
        self.internal.base_command = base_command # type: ignore

    def extract_base_command_text(self) -> str:
        text = ''
        for part in self.task.command.parts:
            if not isinstance(part, str):
                break 
            text += part
        return text
    
    def strip_newline_ignores(self, text: str) -> str:
        text = text.replace('\\\n', '')
        return text
    
    def split_base_command(self, text: str) -> list[str]:
        base_command = text.split()
        base_command = [x.strip() for x in base_command]
        base_command = [x for x in base_command if x != ""]
        return base_command
    
    def strip_set_e(self, base_command: list[str]) -> list[str]:
        for ptr, token in enumerate(base_command):
            if token.lower() == 'set':
                start = deepcopy(ptr)
                end = eat_set_command(ptr, base_command)
                break
        return base_command[:start] + base_command[end:]
    
    def update_inputs(self) -> None:
        # position
        # prefix
        # delim
        # item separator
        # unwrapping expression to janis
        # - wrapped in function(s)?
        # - should these become janis arguments? 
        raise NotImplementedError
    
    def update_arguments(self) -> None:
        raise NotImplementedError
    
    def update_outputs(self) -> None:
        raise NotImplementedError


def eat_set_command(ptr: int, tokens: list[str]) -> int:
    """returns the index of the next token after the set command"""
    ptr += 1
    while ptr < len(tokens) and tokens[ptr].startswith('-'):
        ptr += 1
    if ptr < len(tokens) -1 and tokens[ptr].lower() == 'pipefail': 
        ptr += 1
    return ptr