
from enum import Enum
from textwrap import indent

from typing import Optional

from ..directives import ProcessDirective
from ..casefmt import to_case
from .. import settings

from .inputs.factory import ProcessInput
from .outputs import ProcessOutput
from .. import ordering


def filter_null(iterable):
    if iterable is None:
        return None
    if not isinstance(iterable, list):
        raise Exception(
            f"filter_null: Expected iterable ('{str(iterable)}') to be of type list, received '{type(iterable)}'"
        )

    return [el for el in iterable if el is not None]


class ProcessScriptType(Enum):
    script = "script"
    shell = "shell"
    exec = "exec"


class Process:
    def __init__(
        self,
        name: str,
        script: str,
        script_type: Optional[ProcessScriptType] = None,
        script_quote: Optional[str] = '"',
        inputs: Optional[list[ProcessInput]] = None,
        outputs: Optional[list[ProcessOutput]] = None,
        when: Optional[str] = None,  # TODO unimplemented?
        directives: Optional[list[ProcessDirective]] = None,
        pre_script: Optional[str] = None,
    ):
        self.name = to_case(name, settings.NF_PROCESS_CASE)

        self.script = script
        self.script_type = script_type
        self.script_quote = script_quote

        self.inputs: list[ProcessInput] = inputs or []
        self.outputs: list[ProcessOutput] = outputs or []
        self.directives: list[ProcessDirective] = directives or []
        self.pre_script = pre_script

    def prepare_script(self):
        script_body = str(self.script).strip()
        
        if self.script_type:
            script = ''
            script += f'{self.script_type.value}:\n'
            script += f'{self.pre_script}\n' if self.pre_script else ''
            script += f'{3 * self.script_quote}\n' if self.script_quote else ''
            script += f'{script_body}\n'
            script += f'{3 * self.script_quote}\n' if self.script_quote else ''
            script = indent(script, settings.NF_INDENT)
        else:
            script = indent(script_body, settings.NF_INDENT)
        
        return script

    def prepare_inputs(self):
        if not self.inputs:
            return None
        return indent(
            "input:\n" + "\n".join(i.get_string() for i in self.inputs), 
            settings.NF_INDENT
        )

    def prepare_outputs(self):
        if not self.outputs:
            return None
        return indent(
            "output:\n" + "\n".join(o.get_string() for o in self.outputs),
            settings.NF_INDENT,
        )

    def prepare_directives(self):
        if not self.directives:
            return None
        directives = ordering.order_nf_directives(self.directives)
        return "\n".join(settings.NF_INDENT + d.get_string() for d in directives)

    def get_string(self) -> str:
        components = filter_null(
            [
                self.prepare_directives(),
                self.prepare_inputs(),
                self.prepare_outputs(),
                self.prepare_script(),
            ]
        )
        tool_definition = "\n\n".join(components)

        return f"""\
process {self.name} {{
{tool_definition}
}}
"""
