
from typing import Any

from janis_core import (
    ToolOutput, 
    TOutput, 
    CommandTool, 
    PythonTool,
)

from ....scope import Scope
from ....model.process.outputs import NFProcessOutput

from .factory_cmdtool import CmdtoolProcessOutputFactory
from .factory_pythontool import PythonToolProcessOutputFactory

from ....variables import VariableManager

def gen_nf_process_outputs(
    tool: CommandTool | PythonTool,
    variable_manager: VariableManager, 
    ) -> list[NFProcessOutput]:
    
    process_outputs: list[NFProcessOutput] = []
    
    for out in tool.outputs():
        if isinstance(out, ToolOutput) and isinstance(tool, CommandTool):
            factory = CmdtoolProcessOutputFactory(
                out=out, 
                tool=tool, 
                variable_manager=variable_manager, 
            )
        if isinstance(out, TOutput) and isinstance(tool, PythonTool):
            factory = PythonToolProcessOutputFactory(
                out=out, 
                tool=tool, 
                variable_manager=variable_manager, 
            )
    
        new_output = factory.create()
        process_outputs.append(new_output)
    
    return process_outputs

