

from __future__ import annotations
from typing import Optional
from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.gxtool.command.components.inputs.InputComponent import InputComponent

from ..model import XMLParam
from .components import CommandComponent
from .components.inputs.InputComponent import InputComponent
from .components.inputs.Positional import Positional
from .components.inputs.Flag import Flag
from .components.inputs.Option import Option
from .components.outputs.RedirectOutput import RedirectOutput


class Command:
    def __init__(self):
        # this structure sucks
        self.positionals: dict[int, Positional] = {}
        self.flags: dict[str, Flag] = {}
        self.options: dict[str, Option] = {}
        self.redirect: Optional[RedirectOutput] = None

    def delete_component(self, component: CommandComponent) -> None:
        if isinstance(component, Flag):
            del self.flags[component.prefix]
        elif isinstance(component, Option):
            del self.options[component.prefix]
        elif isinstance(component, RedirectOutput):
            self.redirect = None
        else:
            raise RuntimeError(f"cant delete component type: {type(component)}")

    def gxparam_is_attached(self, gxparam: XMLParam) -> bool:
        components = self.list_components(include_base_cmd=True)
        for component in components:
            if component.gxparam and component.gxparam.name == gxparam.name:
                return True
        return False
    
    def list_components(self, include_base_cmd: bool=True) -> list[CommandComponent]:
        components: list[CommandComponent] = []
        components += self.list_inputs(include_base_cmd=include_base_cmd)
        components += self.list_outputs()
        return components

    def list_inputs(self, include_base_cmd: bool=True) -> list[InputComponent]:
        components: list[InputComponent] = []
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
    