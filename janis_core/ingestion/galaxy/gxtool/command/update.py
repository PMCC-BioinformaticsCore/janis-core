


from __future__ import annotations
from typing import Any, TYPE_CHECKING

if TYPE_CHECKING:
    from .Command import Command

from abc import ABC, abstractmethod

from .components import CommandComponent
from .components.inputs.Positional import Positional
from .components.inputs.Flag import Flag
from .components.inputs.Option import Option
from .components.outputs.RedirectOutput import RedirectOutput
# from .components import factory


def update_command(command: Command, incoming: CommandComponent) -> None:
    updater = select_updater(incoming)
    updater.update(command, incoming)

def select_updater(incoming: CommandComponent) -> Updater:
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

# def update(self, incoming: CommandComponent) -> None:
#     component = self.refine_component(incoming)
#     updater = self.select_updater(component)
#     updater.update(self, component)

# def refine_component(self, incoming: CommandComponent) -> CommandComponent:
#     # migrate incorrect option to flag
#     if isinstance(incoming, Option):
#         if incoming.prefix in self.flags:
#             return factory.flag(incoming.prefix, incoming.gxparam)
#     return incoming





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

