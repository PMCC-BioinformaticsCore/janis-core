from enum import Enum
from abc import ABC, abstractmethod
from typing import Optional, List
from textwrap import indent


def filter_null(iterable):
    if iterable is None:
        return None
    if not isinstance(iterable, list):
        raise Exception(
            f"filter_null: Expected iterable ('{str(iterable)}') to be of type list, received '{type(iterable)}'"
        )

    return [el for el in iterable if el is not None]


def get_nf_string(value):
    if isinstance(value, list):
        return ", ".join(get_nf_string(v for v in value))

    elif isinstance(value, Enum):
        return value.value

    elif isinstance(value, NFBase):
        return value.get_string()

    return str(value)


class NFBase(ABC):
    @abstractmethod
    def get_string(self):
        raise Exception("Subclass must override .get_string() method")


class ImportItem(NFBase):
    def __init__(self, name: str, alias: Optional[str] = None):
        self.name = name
        self.alias = alias

    def get_string(self):
        if self.alias:
            return f"{self.name} as {self.alias}"
        return self.name


class Import(NFBase):
    def __init__(self, items: List[ImportItem], source: str):
        self.items = items
        self.source = source

    def get_string(self):
        if len(self.items) == 0:
            raise Exception(
                f"NF workflow import from '{self.source}' did not contain any items"
            )

        items = "; ".join(i.get_string() for i in self.items)
        return f"include {{ {items} }} from '{self.source}'"


class NFFile(NFBase):
    def __init__(self, imports: List[Import], items: List[NFBase]):
        self.imports = imports
        self.items = items

    def get_string(self):

        components = [f"nextflow.enable.dsl=2"]

        if self.imports:
            components.append("\n".join(i.get_string() for i in self.imports))
        if self.items:
            components.extend(i.get_string() for i in self.items)

        return "\n\n".join(components)


class Channel(NFBase):
    def __init__(self, method, value):
        self.method = method
        self.value = value

    def get_string(self):
        value = get_nf_string(self.value)
        return f"Channel.{self.method}( {value} )"

    @classmethod
    def from_(cls, value):
        return cls("from", value)


class Function(NFBase):
    def __init__(self, name: str, parameters: List[str], definition: str):
        self.name = name
        self.definition = definition
        self.parameters = parameters

    def get_string(self):
        return f"""
def {self.name}({", ".join(self.parameters)}) {{
  {indent(self.definition, '  ')}
}}
"""
