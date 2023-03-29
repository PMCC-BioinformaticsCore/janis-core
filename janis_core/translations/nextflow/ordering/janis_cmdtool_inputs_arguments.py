

from abc import ABC, abstractmethod
from janis_core import CommandTool, ToolArgument, ToolInput
from janis_core.types import Boolean, File
from .. import nfgen_utils

class CmdtoolInsArgsStrategy(ABC):
    @abstractmethod
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        ...
        
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

# class ExposedPriorityStrategy(CmdtoolInsArgsStrategy):
#     def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
#         top: list[ToolInput | ToolArgument] = []
#         bottom: list[ToolInput | ToolArgument] = []
#         for x in ins_args:
#             if isinstance(x, ToolInput) and x.id() in tool.connections:
#                 top.append(x)
#             else:
#                 bottom.append(x)
#         return top + bottom


class PositionStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda a: (a.position or 0))

class NoPrefixPriorityStrategy(CmdtoolInsArgsStrategy):
    def order(self, ins_args: list[ToolInput | ToolArgument], tool: CommandTool) -> list[ToolInput | ToolArgument]:
        return sorted(ins_args, key=lambda a: (a.prefix is None))


ins_args_strategies = [
    # aesthetic
    AlphabeticalStrategy,
    FilePriorityStrategy,
    ComponentTypeStrategy,
    # ExposedPriorityStrategy,
    InsPriorityStrategy,
    # correctness
    NoPrefixPriorityStrategy,
    PositionStrategy
]

def order_cmdtool_inputs_arguments(tool: CommandTool) -> list[ToolInput | ToolArgument]:
    ins_args: list[ToolInput | ToolArgument] = []

    if tool.arguments():
        ins_args += [a for a in tool.arguments() if a.position is not None or a.prefix is not None]
    
    if tool.inputs():
        ins_args += [a for a in tool.inputs() if a.position is not None or a.prefix is not None]
    
    for strategy in ins_args_strategies:
        ins_args = strategy().order(ins_args, tool)
    
    return ins_args

