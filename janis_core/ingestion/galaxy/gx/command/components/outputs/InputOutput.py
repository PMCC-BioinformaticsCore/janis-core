

from __future__ import annotations
from .OutputComponent import OutputComponent

from typing import Any, Optional

from ..CommandComponent import CommandComponent


class InputOutput(OutputComponent):
    def __init__(self, input_component: CommandComponent):
        super().__init__()
        self.input_component = input_component
        self.gxparam = self.input_component.gxparam
        assert(self.gxparam)

    @property
    def name(self) -> str:
        name = ''
        if self.gxparam:
            name = self.gxparam.name
        if not name.startswith('out'):
            name = f'out_{name}'
        return name

    @property
    def default_value(self) -> Any:
        raise NotImplementedError()

    @property
    def optional(self) -> bool:
        if self.forced_optionality is not None:
            return self.forced_optionality
        elif self.gxparam: 
            return self.gxparam.optional
        elif self.input_component.optional:
            return True
        return False

    @property
    def docstring(self) -> Optional[str]:
        if self.gxparam:
            return self.gxparam.docstring
        return f'output created during runtime. file relates to the {self.input_component.name} input'

    def update(self, incoming: Any) -> None:
        raise NotImplementedError()

