

from __future__ import annotations
from typing import Optional

from .OutputComponent import OutputComponent
from ..linux.streams import Stream
from ...components.ValueRecord import ValueRecord


class RedirectOutput(OutputComponent):
    def __init__(self, redirect: str, filepath: str):
        super().__init__()
        self.redirect = redirect
        self.stream: Stream = self.extract_stream()
        self.values: ValueRecord = ValueRecord()
        self.values.add(filepath)

    @property
    def name(self) -> str:
        name: str = 'redirect'
        # get name from galaxy param if available
        if self.gxparam:
            name = self.gxparam.name
        # otherwise, most commonly witnessed option value as name
        elif self.values.most_common_value:
            name = self.values.most_common_value
        if not name.startswith('out'):
            name = f'out_{name}'
        return name
    
    @property
    def text(self) -> str:
        filepath = self.values.most_common_value
        assert(filepath)
        return f'{self.redirect} {filepath}'

    @property
    def optional(self) -> bool:
        return False

    @property
    def docstring(self) -> Optional[str]:
        if self.gxparam:
            return self.gxparam.docstring
        return ''
    
    def is_append(self) -> bool:
        if '>>' in self.redirect:
            return True
        return False

    def extract_stream(self) -> Stream:
        if self.redirect[0] == '2':
            return Stream.STDERR
        return Stream.STDOUT


