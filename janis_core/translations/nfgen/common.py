
import os

from abc import ABC, abstractmethod
from typing import Optional, List
from textwrap import indent
from . import settings
from .casefmt import to_case


def filter_null(iterable):
    if iterable is None:
        return None
    if not isinstance(iterable, list):
        raise Exception(
            f"filter_null: Expected iterable ('{str(iterable)}') to be of type list, received '{type(iterable)}'"
        )

    return [el for el in iterable if el is not None]


class NFBase(ABC):
    @abstractmethod
    def get_string(self) -> str:
        raise Exception("Subclass must override .get_string() method")

class ImportItem(NFBase):
    def __init__(self, name: str, alias: Optional[str] = None):
        self.name = name
        self.alias = alias

    def get_string(self) -> str:
        name = to_case(self.name, settings.NF_PROCESS_CASE)
        if self.alias:
            return f"{name} as {self.alias}"
        return name


class Import(NFBase):
    def __init__(self, items: List[ImportItem], source: str):
        self.items = items
        self.source = source

    def get_string(self) -> str:
        if len(self.items) == 0:
            raise Exception(
                f"NF workflow import from '{self.source}' did not contain any items"
            )

        items = "; ".join(i.get_string() for i in self.items)
        return f"include {{ {items} }} from '{self.source}'"


class NFFile(NFBase):
    def __init__(self, subtype: str, imports: List[Import], items: List[NFBase], name: Optional[str]=None):
        self.subtype = subtype
        self.imports = imports
        self.items = items
        self._name = name

    @property
    def name(self) -> str:
        return self._name  # type: ignore
    
    @property
    def path(self) -> str:
        if self.subtype == 'process':
            directory = 'modules'
        elif self.subtype == 'sub_workflow':
            directory = 'subworkflows'
        elif self.subtype == 'main_workflow':
            directory = ''
        return os.path.join(directory, self.name)

    def get_string(self) -> str:
        components = [f"nextflow.enable.dsl=2"]

        if self.imports:
            components.append("\n".join(i.get_string() for i in self.imports))
        if self.items:
            components.extend(i.get_string() for i in self.items)

        return "\n\n".join(components)
    
class Function(NFBase):
    def __init__(self, name: str, parameters: List[str], definition: str):
        self.name = name
        self.definition = definition
        self.parameters = parameters

    def get_string(self) -> str:
        return f"""
def {self.name}({", ".join(self.parameters)}) {{
  {indent(self.definition, settings.NF_INDENT)}
}}
"""
