
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


def create_nextflow_process_outputs(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> list[ProcessOutput]:
    process_outputs: list[ProcessOutput] = []
    for out in tool.outputs():
        process_outputs.append(create_output(out, tool, sources))
    return process_outputs

def create_output(out: ToolOutput | TOutput, tool: CommandTool | PythonTool, sources: dict[str, Any]) -> ProcessOutput:
    if isinstance(out, ToolOutput) and isinstance(tool, CommandTool):
        factory = CmdtoolProcessOutputFactory(out, tool, sources)
    if isinstance(out, TOutput) and isinstance(tool, PythonTool):
        factory = PythonToolProcessOutputFactory(out, tool, sources)
    return factory.create()
