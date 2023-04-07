
from enum import Enum
from textwrap import indent

from typing import Optional
from dataclasses import dataclass, field

from janis_core import settings
from ... import ordering

from .directives import NFProcessDirective
from .inputs import NFProcessInput
from .outputs import NFProcessOutput

INDENT = settings.translate.nextflow.NF_INDENT


class NFProcessScriptType(Enum):
    SCRIPT = "script"
    SHELL = "shell"
    EXEC = "exec"


@dataclass
class NFProcess:
    name: str
    script: str
    script_type: NFProcessScriptType = NFProcessScriptType.SCRIPT
    script_quote: Optional[str] = '"'
    directives: list[NFProcessDirective] = field(default_factory=list)
    inputs: list[NFProcessInput] = field(default_factory=list)
    outputs: list[NFProcessOutput] = field(default_factory=list)
    
    main_exec: Optional[str] = None
    when: Optional[str] = None  # TODO unimplemented?
    pre_script: Optional[str] = None
    
    @property
    def formatted_directives(self):
        if not self.directives:
            return None
        directives = ordering.order_nf_directives(self.directives)
        return "\n".join(INDENT + d.get_string() for d in directives)

    @property
    def formatted_inputs(self):
        if not self.inputs:
            return None
        return indent(
            "input:\n" + "\n".join(i.get_string() for i in self.inputs), 
            INDENT
        )

    @property
    def formatted_outputs(self):
        if not self.outputs:
            return None
        return indent(
            "output:\n" + "\n".join(o.get_string() for o in self.outputs),
            INDENT,
        )
    
    @property
    def formatted_exec(self) -> Optional[str]:
        if self.main_exec is not None:
            outstr = ''
            outstr += 'exec:'
            if self.main_exec != '':
                outstr += f'\n{self.main_exec}'
            outstr = indent(outstr, INDENT)
            return outstr
        return None
    
    @property
    def formatted_script(self):
        script_body = str(self.script).strip()
        
        if self.script_type == NFProcessScriptType.SCRIPT:
            script = ''
            script += f'{self.script_type.value}:\n'
            script += f'{self.pre_script}\n' if self.pre_script else ''
            script += f'{3 * self.script_quote}\n' if self.script_quote else ''
            script += f'{script_body}\n'
            script += f'{3 * self.script_quote}\n' if self.script_quote else ''
            script = indent(script, INDENT)
        else:
            script = indent(script_body, INDENT)
        
        return script

    def get_string(self) -> str:
        possible_items = [
            self.formatted_directives,
            self.formatted_inputs,
            self.formatted_outputs,
            self.formatted_exec,
            self.formatted_script,
        ]
        final_items = [x for x in possible_items if x is not None]
        tool_definition = "\n\n".join(final_items)
        return f"""\
process {self.name} {{
{tool_definition}
}}
"""
