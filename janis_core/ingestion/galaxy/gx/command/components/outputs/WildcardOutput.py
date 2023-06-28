



from __future__ import annotations
from .OutputComponent import OutputComponent
from typing import Optional


class WildcardOutput(OutputComponent):
    def __init__(self):
        super().__init__()
        self.verified: bool = False

    @property
    def name(self) -> str:
        name = ''
        if self.gxparam:
            name = self.gxparam.name
        if not name.startswith('out'):
            name = f'out_{name}'
        return name

    @property
    def optional(self) -> bool:
        if self.forced_optionality is not None:
            return self.forced_optionality
        elif self.gxparam: 
            return self.gxparam.optional
        return False
   
    @property
    def docstring(self) -> Optional[str]:
        if self.gxparam:
            return self.gxparam.docstring
        return 'output created during runtime. file is collected from working directory'



