

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any, Optional
from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.gx.command.components.inputs.InputComponent import InputComponent

from ..gxtool.param.Param import Param 

from .cmdstr.CommandString import CommandString
from .components.CommandComponent import CommandComponent
from .components.inputs.InputComponent import InputComponent
from .components.inputs.Positional import Positional
from .components.inputs.Flag import Flag
from .components.inputs.Option import Option
from .components.outputs.RedirectOutput import RedirectOutput
from .components import factory


class Updater(ABC):
    command: Command
    incoming: Positional | Flag | Option | RedirectOutput
    
    @abstractmethod
    def update(self, command: Command, incoming: Any) -> None:
        """updates the command's store of CommandComponents with the incoming component"""
        ...

    @abstractmethod
    def should_merge(self) -> bool:
        """determines whether to merge the incoming CommandComponent with an existing one, or add new"""
        ...
    
    @abstractmethod
    def merge(self) -> None:
        """updates via merging the incoming CommandComponent with a known CommandComponent"""
        ...
    
    @abstractmethod
    def add(self) -> None:
        """updates via adding a new CommandComponent"""
        ...


class PositionalUpdater(Updater):
    command: Command
    incoming: Positional

    def update(self, command: Command, incoming: Positional) -> None:
        self.command = command
        self.incoming = incoming
        if self.should_merge():
            self.merge()
        else:
            self.add()

    def should_merge(self) -> bool:
        cmd_pos = self.incoming.cmd_pos
        existing_comp = self.command.get_positional(cmd_pos)
        if existing_comp:
            return True
        return False
    
    def merge(self) -> None:
        cmd_pos = self.incoming.cmd_pos
        existing_comp = self.command.get_positional(cmd_pos)
        if existing_comp:
            existing_comp.update(self.incoming)
    
    def add(self) -> None:
        cmd_pos = self.incoming.cmd_pos
        self.command.positionals[cmd_pos] = self.incoming


class FlagUpdater(Updater):
    command: Command
    incoming: Flag

    def update(self, command: Command, incoming: Flag) -> None:
        self.command = command
        self.incoming = incoming
        if self.should_merge():
            self.merge()
        else:
            self.add()

    def should_merge(self) -> bool:
        query_prefix = self.incoming.prefix
        existing_comp = self.command.get_flag(query_prefix)
        if existing_comp:
            return True
        return False
    
    def merge(self) -> None:
        query_prefix = self.incoming.prefix
        existing_comp = self.command.get_flag(query_prefix)
        if existing_comp:
            existing_comp.update(self.incoming)
    
    def add(self) -> None:
        prefix = self.incoming.prefix
        self.command.flags[prefix] = self.incoming


class OptionUpdater(Updater):
    command: Command
    incoming: Option

    def update(self, command: Command, incoming: Option) -> None:
        self.command = command
        self.incoming = incoming
        if self.should_merge():
            self.merge()
        else:
            self.add()

    def should_merge(self) -> bool:
        query_prefix = self.incoming.prefix
        existing_comp = self.command.get_option(query_prefix)
        if existing_comp:
            return True
        return False
    
    def merge(self) -> None:
        query_prefix = self.incoming.prefix
        existing_comp = self.command.get_option(query_prefix)
        if existing_comp:
            existing_comp.update(self.incoming)
    
    def add(self) -> None:
        prefix = self.incoming.prefix
        self.command.options[prefix] = self.incoming


class RedirectOutputUpdater(Updater):
    command: Command
    incoming: RedirectOutput

    def update(self, command: Command, incoming: RedirectOutput) -> None:
        self.command = command
        self.incoming = incoming
        if self.should_merge():
            pass
        else:
            self.add()

    def should_merge(self) -> bool:
        if self.command.redirect:
            return True
        return False
    
    def merge(self) -> None:
        raise RuntimeError("can't update an output")
    
    def add(self) -> None:
        self.command.redirect = self.incoming


class Command:
    def __init__(self, xmlcmdstr: CommandString):
        self.xmlcmdstr = xmlcmdstr
        self.positionals: dict[int, Positional] = {}
        self.flags: dict[str, Flag] = {}
        self.options: dict[str, Option] = {}
        self.redirect: Optional[RedirectOutput] = None

    def gxparam_is_attached(self, gxparam: Param) -> bool:
        components = self.list_inputs()
        components += self.list_outputs()
        for component in components:
            if component.gxparam and component.gxparam.name == gxparam.name:
                return True
        return False

    def list_inputs(self, include_base_cmd: bool=True) -> list[CommandComponent]:
        components: list[CommandComponent] = []
        if include_base_cmd:
            components += self.get_positionals()
        else:
            components += self.get_non_base_positionals()
        components += self.get_flags()
        components += self.get_options()
        return components

    def list_outputs(self) -> list[RedirectOutput]:
        # just returns redirect component if present.
        # other outputs are handled by ToolFactory
        components: list[RedirectOutput] = []
        if self.redirect:
            components.append(self.redirect)
        return components
    
    def get_positional(self, cmd_pos: int) -> Optional[Positional]:
        if cmd_pos in self.positionals:
            return self.positionals[cmd_pos]
        return None

    def get_flag(self, query_prefix: str) -> Optional[Flag]:
        if query_prefix in self.flags:
            return self.flags[query_prefix]
        return None

    def get_option(self, query_prefix: str) -> Optional[Option]:
        if query_prefix in self.options:
            return self.options[query_prefix]
        return None
    
    def get_positionals(self) -> list[Positional]:
        """returns positionals in sorted order"""
        positions_components = list(self.positionals.items())
        positions_components.sort(key=lambda x: x[0])
        return [p[1] for p in positions_components]
    
    def get_flags(self) -> list[Flag]:
        return list(self.flags.values())
    
    def get_options(self) -> list[Option]:
        return list(self.options.values())

    def get_base_positionals(self) -> list[InputComponent]:
        positionals = self.get_positionals()
        out: list[InputComponent] = []
        for p in positionals:
            if self.positional_is_based(p):
                out.append(p)
            else:
                break
        return out

    def positional_is_based(self, p: Positional) -> bool:
        if p.before_opts:
            if not p.gxparam:
                if len(p.values.unique) == 1:
                    if not expressions.is_var(p.values.unique[0]):
                        return True
        return False

    def get_non_base_positionals(self) -> list[Positional]:
        base_positionals = self.get_base_positionals()
        all_positionals = self.get_positionals()
        return [p for p in all_positionals if p not in base_positionals]

    def update(self, incoming: CommandComponent) -> None:
        component = self.refine_component(incoming)
        updater = self.select_updater(component)
        updater.update(self, component)

    def refine_component(self, incoming: CommandComponent) -> CommandComponent:
        # migrate incorrect option to flag
        if isinstance(incoming, Option):
            if incoming.prefix in self.flags:
                return factory.flag(incoming.prefix, incoming.gxparam)
        # migrate incorrect flag to option
        # if isinstance(incoming, Flag):
        #     if incoming.prefix in self.options:
        #         return factory.option(prefix=incoming.prefix, gxparam=incoming.gxparam)
        return incoming

    def select_updater(self, incoming: CommandComponent) -> Updater:
        match incoming:
            case Positional():
                return PositionalUpdater()
            case Flag():
                return FlagUpdater()
            case Option():
                return OptionUpdater()
            case RedirectOutput():
                return RedirectOutputUpdater()
            case _:
                raise RuntimeError(f'must pass CommandComponent to Command.update(). received {type(incoming)}')

    def set_cmd_positions(self) -> None:
        options_pos: int = self.get_options_position()
        for flag in self.get_flags():
            flag.cmd_pos = options_pos
        for option in self.get_options():
            option.cmd_pos = options_pos
        
        pos_ptr: int = 1
        for positional in self.get_non_base_positionals():
            if pos_ptr == options_pos:
                pos_ptr += 1
            positional.cmd_pos = pos_ptr
            pos_ptr += 1

    def get_options_position(self) -> int:
        """
        returns cmd_pos for options and flags. 
        base command will always occupy cmd_pos == 0
        """
        i: int = 0
        positionals = self.get_non_base_positionals()
        while i < len(positionals) and positionals[i].before_opts:
            i += 1
        return i + 1

    # string representations
    def __str__(self) -> str:
        return f""" \
##### Command #####

positionals: ------------
{'main value':20}{'optional':>10}
{self._get_components_list_as_str('positional')}

flags: ------------
{'prefix':20}{'optional':>10}
{self._get_components_list_as_str('flag')}

options: ------------
{'prefix':30}{'main value':20}{'optional':>10}
{self._get_components_list_as_str('option')}"""

    def _get_components_list_as_str(self, ctype: str='positional') -> str:
        funcmap = {
            'positional': self.get_positionals,
            'flag': self.get_flags,
            'option': self.get_options,
        }
        outstr: str = ''
        for comp in funcmap[ctype]():
            outstr += comp.__str__() + '\n'
        return outstr
    