
from abc import ABC, abstractmethod
from dataclasses import dataclass



@dataclass
class ProcessOutput(ABC):
    name: str
    is_optional: bool

    @abstractmethod
    def get_string(self) -> str:
        ...
    
    def emit(self, alias: str) -> str:
        return f', emit: {alias}'
    
    @property
    def optional(self) -> str:
        if self.is_optional:
            return ', optional: true'
        return ''


@dataclass
class StdoutProcessOutput(ProcessOutput):

    def get_string(self) -> str:
        return f'stdout{self.optional}{self.emit(self.name)}'


@dataclass
class ValProcessOutput(ProcessOutput):
    expression: str

    def get_string(self) -> str:
        return f'val {self.expression}{self.optional}{self.emit(self.name)}'


@dataclass
class PathProcessOutput(ProcessOutput):
    expression: str

    def get_string(self) -> str:
        return f'path {self.expression}{self.optional}{self.emit(self.name)}'


@dataclass
class TupleProcessOutput(ProcessOutput):
    qualifiers: list[str]
    expressions: list[str]

    @property
    def fields(self) -> str:
        out: str = ''
        for qual, expr in zip(self.qualifiers, self.expressions):
            out += f'{qual}({expr}), '
        out = out.rstrip(', ') # strip off the last comma & space
        return out
    
    def get_string(self) -> str:
        return f'tuple {self.fields}{self.optional}{self.emit(self.name)}'


@dataclass
class SecondariesArrayProcessOutput(ProcessOutput):
    qualifiers: list[str]
    expressions: list[str]
    names: list[str]

    def get_string(self) -> str:
        lines: list[str] = []
        for qual, expr, name in zip(self.qualifiers, self.expressions, self.names):
            line = f'{qual} {expr}{self.optional}{self.emit(name)}'
            lines.append(line)
        return '\n'.join(lines)