

from typing import Optional
from abc import ABC, abstractmethod
from dataclasses import dataclass



"""
tuple, path, val


PATH
path x, stageAs: 'data.txt'
"""


class ProcessInput(ABC):
    @abstractmethod
    def get_string(self) -> str:
        ...


@dataclass
class ValProcessInput(ProcessInput):
    name: str

    def get_string(self) -> str:
        return f'val {self.name}'


@dataclass
class PathProcessInput(ProcessInput):
    name: str
    stage_as: Optional[str]=None

    @property
    def attributes(self) -> list[str]:
        out: list[str] = []
        if self.stage_as:
            attr = f"stageAs: '{self.stage_as}'"
            out.append(attr)
        return out

    def get_string(self) -> str:
        base = f'path {self.name}'
        if self.attributes:
            attrs = ', '.join(self.attributes)
            declaration = f'{base}, {attrs}'
        else:
            declaration = base
        return declaration


@dataclass
class TupleProcessInput(ProcessInput):
    names: list[str]
    qualifiers: list[str]

    def get_string(self) -> str:
        declaration = 'tuple '
        for name, qual in zip(self.names, self.qualifiers):
            declaration += f'{qual}({name}), '
        declaration = declaration.rstrip(', ')
        return declaration










# path


    # val = "val"
    # env = "env"
    # file = "file"
    # path = "path"
    # tuple = "tuple"
    # stdout = "stdout"


# class ProcessInput:
#     def __init__(
#         self,
#         qualifier: Qualifier,
#         name: str,
#         attributes: Optional[str | list[str]] = None,
#     ):
#         self.qualifier = qualifier
#         self.name = name
#         self.attributes = attributes

#     def get_string(self) -> str:
#         els = [self.qualifier.value, self.name]

#         if self.attributes:
#             if isinstance(self.attributes, list):
#                 els.extend(str(a) for a in self.attributes)
#             else:
#                 els.append(str(self.attributes))

#         return " ".join(str(e) for e in els).strip()







