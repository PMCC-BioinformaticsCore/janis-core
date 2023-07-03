


from __future__ import annotations
from typing import Any, TYPE_CHECKING
from abc import ABC, abstractmethod

if TYPE_CHECKING:
    from .Command import Command

from janis_core.ingestion.galaxy.gxtool.model import XMLSelectParam
from .components import CommandComponent
from .components.inputs.Positional import Positional
from .components.inputs.Flag import Flag
from .components.inputs.Option import Option
from .components.outputs.RedirectOutput import RedirectOutput


def update_command(command: Command, incoming: CommandComponent, final_pass: bool=False) -> None:
    if should_update(command, incoming, final_pass):
        updater = select_updater(incoming)
        updater.update(command, incoming)


### CHECKING WHETHER TO UPDATE COMMAND ###

def should_update(command: Command, incoming: CommandComponent, final_pass: bool=False) -> bool:
    """
    exists so we know whether to update the command. 
    this function ensures we are updating the command with GOOD information, not BAD. 
    this function assumes components are already refined (context aware) in the epath annotation phase.
    
    in some cases its easy:
        - if there is no existing component similar to the existing component, we always update
        - if there is an existing component but it's essentially the same, we always update (merge)
        - if the incoming component is a positional, we always update as positionals are only extracted once. 
          (there will never be a duplicate already in the command)
    in other cases its hard:
        - when there is a similar existing component, but they're not actually the same (eg flag vs option)
        - we will need to decide which to keep.
    """
    # 
    # only add positionals on final cmdstr pass
    if isinstance(incoming, Positional):
        if final_pass:
            return True
        return False

    # select params which are list of simple flags
    if isinstance(incoming, Flag) and incoming.gxparam and isinstance(incoming.gxparam, XMLSelectParam):
        return True
    
    # get existing components similar to incoming component
    similar_components = get_similar_components(command, incoming)
    
    # first observation of this component
    if not similar_components:
        return True
    
    # single similar component exists
    elif len(similar_components) == 1:
        # same component
        existing = similar_components[0]
        if components_are_the_same(existing, incoming):
            return True
        
        # different components
        best_component = select_best_component_from_mismatch(existing, incoming)
        if best_component.uuid == existing.uuid:
            return False
        elif best_component.uuid == incoming.uuid:
            command.delete_component(existing)
            return True
        
    # multiple similar components exist
    elif len(similar_components) == 1:
        # ...shit hope this doesn't happen
        raise NotImplementedError

    return False

def get_similar_components(command: Command, incoming: CommandComponent) -> list[CommandComponent]:
    all_components = command.list_components(include_base_cmd=True)
    similar_components: list[CommandComponent] = []
    similar_components_ids: set[str] = set()
    
    if incoming.gxparam:
        for comp in all_components:
            if comp.uuid in similar_components_ids:
                continue
            if comp.gxparam and comp.gxparam.name == incoming.gxparam.name:
                similar_components.append(comp)
                similar_components_ids.add(comp.uuid)
                
    if isinstance(incoming, Flag | Option):
        for comp in all_components:
            if comp.uuid in similar_components_ids:
                continue
            if isinstance(comp, Flag | Option) and comp.prefix == incoming.prefix:
                similar_components.append(comp)
                similar_components_ids.add(comp.uuid)

    return similar_components

def components_are_the_same(existing: CommandComponent, incoming: CommandComponent) -> bool:
    if isinstance(existing, Flag) and isinstance(incoming, Flag):
        if existing.prefix == incoming.prefix:
            return True
        elif have_same_gxparam(existing, incoming):
            return True
    
    elif isinstance(existing, Option) and isinstance(incoming, Option):
        if existing.prefix == incoming.prefix:
            return True
        elif have_same_gxparam(existing, incoming):
            return True
    
    elif isinstance(existing, RedirectOutput) and isinstance(incoming, RedirectOutput):
        if have_same_gxparam(existing, incoming):
            return True
    
    return False

def have_same_gxparam(comp1: CommandComponent, comp2: CommandComponent) -> bool:
    if comp1.gxparam and comp2.gxparam:
        if comp1.gxparam.name == comp2.gxparam.name:
            return True
    return False

def select_best_component_from_mismatch(existing: CommandComponent, incoming: CommandComponent) -> CommandComponent:
    """components are not the same. there is some sort of mismatch. select the best one."""

    # incoming is option
    if isinstance(incoming, Option) and isinstance(existing, Flag):
        # will either have the same prefix or the same gxparam (or both)
        # preference the option
        if existing.confidence.value > incoming.confidence.value:
            return existing
        return incoming

    elif isinstance(incoming, Option) and isinstance(existing, RedirectOutput):
        # WTF?
        raise NotImplementedError
    
    # incoming is flag
    elif isinstance(incoming, Flag) and isinstance(existing, Option):
        # preference the option
        if incoming.confidence.value > incoming.confidence.value:
            return incoming
        return existing

    elif isinstance(incoming, Flag) and isinstance(existing, RedirectOutput):
        # WTF?
        raise NotImplementedError
    
    # incoming is redirect
    elif isinstance(incoming, RedirectOutput) and isinstance(existing, RedirectOutput):
        if incoming.gxparam and not existing.gxparam:
            return incoming
        return existing
    
    elif isinstance(incoming, RedirectOutput) and isinstance(existing, Option):
        # WTF?
        raise NotImplementedError

    elif isinstance(incoming, RedirectOutput) and isinstance(existing, Flag):
        # WTF?
        raise NotImplementedError

    raise NotImplementedError(f'cant compare {type(existing)} and {type(incoming)}')




### DOING UPDATES ###

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


class Updater(ABC):
    command: Command
    incoming: Positional | Flag | Option | RedirectOutput
    
    @abstractmethod
    def update(self, command: Command, incoming: Any) -> None:
        """updates the command's store of CommandComponents with the incoming component"""
        ...


class PositionalUpdater(Updater):
    command: Command
    incoming: Positional

    def update(self, command: Command, incoming: Positional) -> None:
        command.positionals[incoming.cmd_pos] = incoming


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
        existing_flag = self.command.get_flag(query_prefix)
        if existing_flag:
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

