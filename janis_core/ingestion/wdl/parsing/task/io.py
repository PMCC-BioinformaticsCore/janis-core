

import WDL 

from janis_core import ToolInput, ToolOutput, CommandToolBuilder
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory

from ..types import parse_type
from ..expressions import parse_expr
from ..EntityParser import TaskParser



class TaskInputParser(TaskParser):

    def __init__(self, task: WDL.Tree.Task, cmdtool: CommandToolBuilder, wdl_inp: WDL.Tree.Decl) -> None:
        super().__init__(task, cmdtool)
        self.wdl_inp = wdl_inp

    def do_parse(self) -> ToolInput:
        default = None
        if self.wdl_inp.expr:
            default = parse_expr(self.wdl_inp.expr, self.task, self.cmdtool)
        tinput = ToolInput(self.wdl_inp.name, parse_type(self.wdl_inp.type, self.task, uuid=self.wdl_inp.name), default=default)
        return tinput
    
    def fallback(self) -> ToolInput:
        msg = f'Error parsing tool input: {self.wdl_inp.name}'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return ToolInput(self.wdl_inp.name, parse_type(self.wdl_inp.type, self.task, uuid=self.wdl_inp.name))


class TaskOutputParser(TaskParser):

    def __init__(self, task: WDL.Tree.Task, cmdtool: CommandToolBuilder, wdl_out: WDL.Tree.Decl) -> None:
        super().__init__(task, cmdtool)
        self.wdl_out = wdl_out

    def do_parse(self) -> ToolOutput:
        if self.wdl_out.expr is None:
            raise Exception(f"Output {self.wdl_out.name} has no expression")
        sel = parse_expr(self.wdl_out.expr, self.task, self.cmdtool)
        tout = ToolOutput(self.wdl_out.name, parse_type(self.wdl_out.type, self.task, uuid=self.wdl_out.name), selector=sel)
        return tout
    
    def fallback(self) -> ToolOutput:
        msg = 'Error parsing tool output'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return ToolOutput(self.wdl_out.name, parse_type(self.wdl_out.type, self.task, uuid=self.wdl_out.name))

