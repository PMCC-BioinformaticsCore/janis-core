

from abc import ABC, abstractmethod
from janis_core import CommandTool, ToolArgument, ToolInput
from janis_core.types import Boolean

class CmdtoolInsArgsStrategy(ABC):
    @abstractmethod
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        ...

class PositionStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda x: x.position or 0)

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
            else:
                if any([delim in x.value for delim in [' ', ':', '=']]):
                    flags.append(x)
                else:
                    options.append(x)

        return positionals + options + flags

class InsPriorityStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda x: isinstance(x, ToolInput), reverse=True)

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
    InsPriorityStrategy,
    ComponentTypeStrategy,
    ExposedPriorityStrategy,
    PositionStrategy
]

def cmdtool_inputs_arguments(tool: CommandTool) -> list[ToolInput | ToolArgument]:
    ins_args: list[ToolInput | ToolArgument] = tool.arguments() or [] + tool.inputs()
    for strategy in ins_args_strategies:
        ins_args = strategy().order(ins_args, tool)
    return ins_args