

from abc import ABC, abstractmethod
from typing import Optional

from ..components import InputComponent
from ..components import Positional
from ..components import Flag
from ..components import Option


class ComponentOrderingStrategy(ABC):
    @abstractmethod
    def order(self, components: list[InputComponent], disallow: Optional[list[str]]=None) -> list[InputComponent]:
        """annotates components with command line positions"""
        ...

# TODO WTF epath cmd_pos orderings for components vs Command() orderings for components
class SimplifiedComponentOrderingStrategy(ComponentOrderingStrategy):
    def order(self, components: list[InputComponent], disallow: Optional[list[str]]=None) -> list[InputComponent]:
        start_pos_count = self.count_starting_positionals(components)
        self.assign_positional_positions(components, start_pos_count)
        self.assign_flag_option_positions(components, start_pos_count)
        return components

    def count_starting_positionals(self, components: list[InputComponent]) -> int:
        i: int = 0
        while i < len(components) and isinstance(components[i], Positional):
            i += 1
        return i 

    def assign_positional_positions(self, components: list[InputComponent], start_pos_count: int) -> None:
        positionals = [c for c in components if isinstance(c, Positional)]
        for i, comp in enumerate(positionals):
            if i >= start_pos_count:
                comp.cmd_pos = i + 1
            else:
                comp.cmd_pos = i
                comp.before_opts = True

    def assign_flag_option_positions(self, components: list[InputComponent], start_pos_count: int) -> None:
        options_flags = [c for c in components if isinstance(c, Option) or isinstance(c, Flag)]
        for comp in options_flags:
            comp.cmd_pos = start_pos_count


class RealisticComponentOrderingStrategy(ComponentOrderingStrategy):
    def order(self, components: list[InputComponent], disallow: Optional[list[str]]=None) -> list[InputComponent]:
        # wayyy more complicated
        # 
        pos: int = 0
        for component in components:
            component.cmd_pos = pos
            if isinstance(component, Positional):
                pos += 1
        return components


