
from __future__ import annotations

import os
from typing import Optional, Any
from textwrap import indent 

from janis_core import settings

from ... import naming
from ...channels import Channel


class NFFile:
    def __init__(self, subtype: str, items: list[Any], name: str):
        self.subtype = subtype
        self.items = items
        self.name = name

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
        if self.items:
            components.extend(i.get_string() for i in self.items)
        return "\n\n".join(components)


class NFImportItem:
    def __init__(self, name: str, alias: Optional[str] = None):
        self.name = name
        self.alias = alias

    def get_string(self) -> str:
        name = naming.constructs.gen_varname_process(self.name)
        if self.alias:
            return f"{name} as {self.alias}"
        return name


class NFImport:
    def __init__(self, method: str, items: list[NFImportItem], source: str):
        self.method = method
        self.items = items
        self.source = source

    def get_string(self) -> str:
        if len(self.items) == 0:
            raise Exception(
                f"NF workflow import from '{self.source}' did not contain any items"
            )

        items = "; ".join(i.get_string() for i in self.items)
        return f"{self.method} {{ {items} }} from '{self.source}'"


class NFFunction:
    """its actually a groovy function but whatever"""
    def __init__(self, name: str, parameters: list[str], definition: str):
        self.name = name
        self.definition = definition
        self.parameters = parameters

    def get_string(self) -> str:
        return f"""
def {self.name}({", ".join(self.parameters)}) {{
  {indent(self.definition, settings.translate.nextflow.NF_INDENT)}
}}
"""

class NFFunctionsBlock:
    def __init__(self, functions: list[str]) -> None:
        self.functions = functions

    def get_string(self) -> str:
        return '\n\n'.join(self.functions)

class NFImportsBlock:
    def __init__(self, imports: list[str], declarations: list[str]) -> None:
        # assert(len(methods) == len(imports))
        # self.methods = methods
        self.imports = imports
        self.declarations = declarations

    def get_string(self) -> str:
        # imports = [f'{meth} {imp}' for meth, imp in zip(self.methods, self.imports)]
        imports = [f'{imp}' for imp in self.imports]
        declarations = [f'def {dec}' for dec in self.declarations]
        imports = '\n'.join(imports)
        declarations = '\n'.join(declarations)
        return f'{imports}\n{declarations}'


class NFChannelsBlock: 
    def __init__(self, channels: list[Channel]) -> None:
        self.channels = channels
    
    def get_string(self) -> str:
        outstr = ''
        channels = self.channels
        channels = [ch for ch in channels if ch.define]
        width_col_1 = max([ch.width for ch in channels])
        for ch in channels:
            outstr += f'{ch.name:<{width_col_1}} = {ch.get_string()}\n'
        return outstr