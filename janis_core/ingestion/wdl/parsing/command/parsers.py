

from abc import ABC, abstractmethod, abstractproperty
from typing import Any
import regex as re
import WDL

from janis_core import CommandToolBuilder, Selector, StringFormatter
from janis_core import ToolArgument
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory
from ..expressions import parse_expr
from ..patterns import FLAG, OPTION
from .cmdline import CmdLine
from ..EntityParser import TaskParser


class CommandParser(TaskParser):

    def __init__(self, task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> None:
        super().__init__(task, cmdtool)
        self.base_command = None
        self.env_vars: dict[str, Any] = {}
        self.files_to_create: dict[str, Any] = {}
        self.directories_to_create: list[str | Selector] = []


class NativeSimpleParser(CommandParser):

    def do_parse(self) -> None:
        raise NotImplementedError
    
    def fallback(self) -> None:
        # msg = f'Error parsing tool command natively. Fellback to using shell command.'
        # log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return None


class NativeArgumentParser(CommandParser):

    def do_parse(self) -> None:
        raise NotImplementedError
    
    def fallback(self) -> None:
        # msg = f'Error parsing tool command natively. Fellback to using shell command.'
        # log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return None
 

class ShellCommandParser(CommandParser):

    def do_parse(self) -> None:
        self.base_command = self.parse_base_command()
        self.env_vars = self.parse_env_vars()
        self.files_to_create = self.parse_files_to_create()
        self.directories_to_create = self.parse_dirs_to_create()
        # mark as shell script cmdtool if successful parse
        self.cmdtool.is_shell_script = True
    
    def parse_base_command(self) -> list[str]:
        return ["sh", "script.sh"]
    
    def parse_env_vars(self) -> dict[str, Any]:
        return {}
    
    def parse_files_to_create(self) -> dict[str, Any]:
        res, success = parse_expr(self.task.command, self.task, self.cmdtool)
        if not success:
            print(str(res))
            raise RuntimeError
        return {"script.sh": res}
    
    def parse_dirs_to_create(self) -> list[str | Selector]:
        return []
    
    def fallback(self) -> None:
        msg = f'Error parsing tool command. Output is not correct.'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return None 








input_types = (
    WDL.Type.Boolean,
    WDL.Type.Int,
    WDL.Type.Float,
    WDL.Type.String,
    WDL.Type.File,
    WDL.Type.Directory,
    WDL.Type.Array,
    WDL.Type.Map,
    WDL.Type.Pair,
    WDL.Type.StructInstance
)

class ComponentParser(ABC):
    def __init__(self, task: WDL.Tree.Task, ctoken: str | WDL.Expr.Placeholder, ntoken: str | WDL.Expr.Placeholder) -> None:
        self.task = task
        self.ctoken = ctoken
        self.ntoken = ntoken
    
    @abstractmethod
    def passes_check(self) -> bool: ...
    
    @abstractmethod
    def update_command_tool(self, internal: CommandToolBuilder, ptr: int) -> int: ...

    def is_input(self, token: Any) -> bool:
        if self.task.inputs is None:
            return False
        if not isinstance(token, WDL.Expr.Ident):
            return False
        for inp in self.task.inputs:
            if inp.name == token.name:
                return True
        return False
        

class PrefixExprOptionParser(ComponentParser):

    def passes_check(self) -> bool:
        if isinstance(self.ctoken, str) and isinstance(self.ntoken, WDL.Expr.Placeholder):
            if re.match(OPTION, self.ctoken):
                if isinstance(self.ntoken, WDL.Expr.Placeholder):
                    if isinstance(self.ntoken.expr, WDL.Expr.Get):
                        if self.is_input(self.ntoken.expr.expr):
                            return True
        return False
    
    def update_command_tool(self, internal: CommandToolBuilder, ptr: int) -> int:
        raise NotImplementedError


class ExprOptionParser(ComponentParser):

    def passes_check(self) -> bool:
        return False
    
    def update_command_tool(self, internal: CommandToolBuilder, ptr: int) -> int:
        raise NotImplementedError


class FlagParser(ComponentParser):

    def passes_check(self) -> bool:
        if isinstance(self.ctoken, WDL.Expr.Placeholder):
            if isinstance(self.ctoken.expr, WDL.Expr.Get):
                if self.is_input(self.ctoken.expr.expr):
                    if isinstance(self.ctoken.expr.type, WDL.Type.Boolean):
                        return True
        return False
    
    def update_command_tool(self, internal: CommandToolBuilder, ptr: int) -> int:
        input_name = str(self.ctoken.expr.expr.name)
        tinput = [inp for inp in internal._inputs if inp.id() == input_name][0]
        truevalue = self.ctoken.options['true'] if self.ctoken.options['true'] != '' else None
        falsevalue = self.ctoken.options['false'] if self.ctoken.options['false'] != '' else None
        if truevalue is not None and falsevalue is not None:
            raise RuntimeError
        elif truevalue is not None:
            tinput.prefix = truevalue
        else:
            tinput.prefix = falsevalue
        tinput.position = ptr + 1
        return 1


class PositionalParser(ComponentParser):

    def passes_check(self) -> bool:
        if isinstance(self.ctoken, WDL.Expr.Placeholder):
            if isinstance(self.ctoken.expr, WDL.Expr.Get):
                if self.is_input(self.ctoken.expr.expr):
                    return True
        return False
    
    def update_command_tool(self, internal: CommandToolBuilder, ptr: int) -> int:
        input_name = str(self.ctoken.expr.expr.name)
        tinput = [inp for inp in internal._inputs if inp.id() == input_name][0]
        tinput.position = ptr + 1
        return 1  # move to next token


class ArgumentParser(ComponentParser):
    # not linked to input

    def passes_check(self) -> bool:
        return False
    
    def update_command_tool(self, internal: CommandToolBuilder, ptr: int) -> int:
        raise NotImplementedError



class NativeCommandParser:
    """
    parses WDL command into a CommandToolBuilder object.
    does not use ShellCommandRequirement. 
    does not create shell script. 

    Want to identify:
    - base command
    - details of inputs (position, prefix, delim, item separator)
    - arguments (from strings not related to inputs)
    """
    def __init__(self, internal: CommandToolBuilder, task: WDL.Tree.Task, cmds: list[CmdLine]):
        self.internal = internal
        self.task = task
        self.cmds = cmds
        self.cmdline = cmds[0]

    def parse(self) -> None:
        self.validate_parser()
        self.update_base_command()
        self.greedy_parse()
        # self.detect_options()
        # self.detect_flags()
        # self.detect_positionals()
        # self.detect_arguments()

    def validate_parser(self) -> None:
        local_cmds = [c for c in self.cmds if not c.is_set_command]
        assert len(local_cmds) == 1                     # only allow single cmd
        # TODO expand this to multiple commands/tools? 
        assert isinstance(local_cmds[0].elems[0], str)  # cmd must start with string
        assert local_cmds[0].elems[0] != ''             # cmd must not start with ''
        
        # cmd can't have weird chars
        for elem in local_cmds[0].elems: 
            if isinstance(elem, str):
                assert not re.match(r'[)(}{:]', elem)
        
        self.cmdline = local_cmds[0]

    def update_base_command(self) -> None:
        text = self.cmdline.elems[0]
        assert isinstance(text, str)
        words = text.split(' ')
        i = 0
        while i < len(words):
            if re.match(FLAG, words[i]) or re.match(OPTION, words[i]):
                break
            i += 1
        
        # set base command
        base_cmd = words[:i]
        self.internal.base_command = base_cmd # type: ignore
        
        # update cmdline elems
        rest_line = ' '.join(words[i:])
        if rest_line == '':
            if len(self.cmdline.elems) > 1:
                self.cmdline.elems = self.cmdline.elems[1:]
            else:
                self.cmdline.elems = []
        else:
            self.cmdline.elems[0] = rest_line

    def greedy_parse(self) -> None:
        tokens = self.get_tokens()
        parsers = [
            PrefixExprOptionParser,
            ExprOptionParser,
            FlagParser,
            PositionalParser,
            ArgumentParser,
        ]

        ptr = 0
        while ptr < len(tokens):
            # check not about to go out of bounds
            if ptr == len(tokens) - 1:
                ctoken, ntoken = tokens[ptr], None
            else:
                ctoken, ntoken = tokens[ptr], tokens[ptr+1]
            
            parsed = False
            for parser_c in parsers:
                parser = parser_c(self.task, ctoken, ntoken)
                if parser.passes_check():
                    ptr += parser.update_command_tool(self.internal, ptr)
                    parsed = True
                    break
            assert parsed

    def get_tokens(self) -> list[str | WDL.Expr.Placeholder]:
        tokens = []
        for elem in self.cmdline.elems:
            if isinstance(elem, str):
                tokens += elem.split(' ')
            else:
                tokens += [elem]
        return tokens

    # def detect_options(self) -> None:
    #     raise NotImplementedError
    
    # def detect_flags(self) -> None:
    #     raise NotImplementedError
    
    # def detect_positionals(self) -> None:
    #     raise NotImplementedError
    
    # def detect_arguments(self) -> None:
    #     raise NotImplementedError

    # position
    # prefix
    # delim
    # item separator
    # unwrapping expression to janis
    # - wrapped in function(s)?
    # - should these become janis arguments?


    
class ShellCommandParserDep:
    """
    parses WDL command into a CommandToolBuilder object.
    uses ShellCommandRequirement. 
    does not create shell script. 
    """
    def __init__(self, internal: CommandToolBuilder, task: WDL.Tree.Task, cmds: list[CmdLine]):
        self.internal = internal
        self.task = task
        self.cmds = cmds
        
    @property
    def can_parse(self) -> bool:
        pass

    def parse(self) -> None:
        self.update_inputs()
        self.update_arguments()
        self.update_base_command()

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
    
    def update_base_command(self) -> None:
        raise NotImplementedError

