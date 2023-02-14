

from dataclasses import dataclass
from typing import Optional

from ..tokens import Token
from ..tokens import spawn_end_sentinel
from ..components import CommandComponent
from ..components import StreamMerge
from ..components import Tee

from .ComponentOrderingStrategy import SimplifiedComponentOrderingStrategy


@dataclass
class EPathPosition:
    ptr: int
    token: Token
    ignore: bool = False
    component: Optional[CommandComponent] = None


class ExecutionPath:
    id: int

    def __init__(self, tokens: list[Token]):
        self.positions: list[EPathPosition] = self._init_positions(tokens)
        self.tokens_to_excise: list[int] = []

    def get_components(self) -> list[CommandComponent]:
        components = self._get_component_list()
        ordering_strategy = SimplifiedComponentOrderingStrategy()
        return ordering_strategy.order(components)
    
    def _get_component_list(self) -> list[CommandComponent]:
        """
        gets the CommandComponents in this EPath (unique)
        preserves ordering
        """
        ignore = [Tee, StreamMerge]
        out: list[CommandComponent] = []
        for position in self.positions:
            if position.component and type(position.component) not in ignore and position.component not in out:
                out.append(position.component)
        return out
    
    def _init_positions(self, tokens: list[Token]) -> list[EPathPosition]:
        positions = [EPathPosition(i, token) for i, token in enumerate(tokens)]
        end_sentinel = spawn_end_sentinel()
        positions.append(EPathPosition(len(positions), end_sentinel))
        return positions

    def __str__(self) -> str:
        out: str = '\n'
        for position in self.positions:
            out += f'{str(position.ptr):<3}{position.token.text[:39]:<40}{position.token.ttype:<35}\n'
        return out


