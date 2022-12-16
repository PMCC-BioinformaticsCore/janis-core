
from typing import Any

from janis_core import (
    ToolOutput, 
    TOutput, 
    CommandTool, 
    PythonTool,
)

from .model import ProcessOutput
from .factory_cmdtool import CmdtoolProcessOutputFactory
from .factory_pythontool import PythonToolProcessOutputFactory


def create_output(out: ToolOutput | TOutput, tool: CommandTool | PythonTool, sources: dict[str, Any]) -> ProcessOutput:
    if isinstance(out, ToolOutput) and isinstance(tool, CommandTool):
        factory = CmdtoolProcessOutputFactory(out, tool, sources)
    if isinstance(out, TOutput) and isinstance(tool, PythonTool):
        factory = PythonToolProcessOutputFactory(out, tool, sources)
    return factory.create()
