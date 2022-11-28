

from abc import ABC, abstractmethod
from janis_core import CommandTool, ToolArgument, ToolInput
from janis_core.workflow.workflow import InputNode
from janis_core.types import Boolean, File
from . import nfgen_utils


class CmdtoolInsArgsStrategy(ABC):
    @abstractmethod
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        ...

class PositionStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda x: x.position if x.position else 0)
        

class AlphabeticalStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda x: x.prefix or 'zzz')

class ComponentTypeStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        positionals: list[ToolInput | ToolArgument] = []
        flags: list[ToolInput | ToolArgument] = []
        options: list[ToolInput | ToolArgument] = []

        for x in ins_args:
            # positionals
            if not x.prefix:
                positionals.append(x)
            
            # flag or opt tool inputs
            elif isinstance(x, ToolInput):
                if isinstance(x.input_type, Boolean):
                    flags.append(x)
                else:
                    options.append(x)
            
            # flag or opt tool arguments
            elif isinstance(x, ToolArgument):
                if x.value is None:
                    flags.append(x)
                else:
                    options.append(x)

        return positionals + options + flags

class InsPriorityStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda x: isinstance(x, ToolInput), reverse=True)

class FilePriorityStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        top: list[ToolInput | ToolArgument] = []
        bottom: list[ToolInput | ToolArgument] = []
        for elem in ins_args:
            if isinstance(elem, ToolInput):
                dtype = nfgen_utils.get_base_type_task_input(elem)
                if isinstance(dtype, File):
                    top.append(elem)
                else:
                    bottom.append(elem)
            else:
                bottom.append(elem)
        return top + bottom

class ExposedPriorityStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        top: list[ToolInput | ToolArgument] = []
        bottom: list[ToolInput | ToolArgument] = []
        for x in ins_args:
            if isinstance(x, ToolInput) and x.id() in tool.connections:
                top.append(x)
            else:
                bottom.append(x)
        return top + bottom

ins_args_strategies = [
    AlphabeticalStrategy,
    FilePriorityStrategy,
    ComponentTypeStrategy,
    ExposedPriorityStrategy,
    InsPriorityStrategy,
    PositionStrategy
]

def cmdtool_inputs_arguments(tool: CommandTool) -> list[ToolInput | ToolArgument]:
    ins_args: list[ToolInput | ToolArgument] = []
    ins_args += tool.inputs()
    if tool.arguments():
        ins_args += tool.arguments()
    for strategy in ins_args_strategies:
        ins_args = strategy().order(ins_args, tool)
    return ins_args


class WinpStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        ...

class AlphabeticalWinpStrategy(WinpStrategy):
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        return sorted(inputs, key=lambda x: x.id())

class FileWinpStrategy(WinpStrategy):
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        return sorted(inputs, key=lambda x: isinstance(x, File), reverse=True)

class MandatoryWinpStrategy(WinpStrategy):
    def order(self, inputs: list[InputNode]) -> list[InputNode]:
        return sorted(inputs, key=lambda x: x.datatype.optional == True)

workflow_input_strategies = [
    #AlphabeticalWinpStrategy, 
    FileWinpStrategy,
    MandatoryWinpStrategy,
]

def workflow_inputs(inputs: list[InputNode]) -> list[InputNode]:
    for strategy in workflow_input_strategies:
        inputs = strategy().order(inputs)
    return inputs



# essentially the same as above, but has to be different because 
# no shared interface for Workflow and CommandTool (with datatype etc)
class ToolStrategy(ABC):
    @abstractmethod
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        ...

class AlphabeticalToolStrategy(ToolStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        return sorted(inputs, key=lambda x: x.id())

class FileToolStrategy(ToolStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        return sorted(inputs, key=lambda x: isinstance(x, File), reverse=True)

class MandatoryToolStrategy(ToolStrategy):
    def order(self, inputs: list[ToolInput]) -> list[ToolInput]:
        return sorted(inputs, key=lambda x: x.input_type.optional == True)

tool_input_strategies = [
    #AlphabeticalToolStrategy, 
    FileToolStrategy,
    MandatoryToolStrategy,
]

def tool_inputs(inputs: list[ToolInput]) -> list[ToolInput]:
    for strategy in tool_input_strategies:
        inputs = strategy().order(inputs)
    return inputs