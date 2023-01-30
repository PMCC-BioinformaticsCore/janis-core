
from typing import Any, Tuple, Optional

from janis_core import (
    CommandTool,
    ToolArgument,
    ToolInput,
)

from ...unwrap import unwrap_expression
from ...scope import Scope
from ... import ordering
from ... import settings
from ... import nfgen_utils as nfgen_utils
from .. import inputs
from .ScriptFormatter import ScriptFormatter


def gen_script_for_cmdtool(
    tool: CommandTool,
    stdout_filename: str,
    scope: Scope,
    referenced_variables: set[str],
    sources: dict[str, Any],
) -> Tuple[Optional[str], str]:
    return ProcessScriptGenerator(
        tool=tool,
        stdout_filename=stdout_filename,
        scope=scope,
        referenced_variables=referenced_variables,
        sources=sources,
    ).generate()



class ProcessScriptGenerator:
    def __init__(
        self,
        tool: CommandTool, 
        stdout_filename: str,
        scope: Scope,
        referenced_variables: set[str],
        sources: Optional[dict[str, Any]]=None,
    ):
        self.tool = tool
        self.scope = scope
        self.referenced_variables = referenced_variables
        self.process_name = scope.labels[-1]
        self.stdout_filename = stdout_filename

        self.sources = sources if sources is not None else {}
        self.process_inputs = inputs.get_process_inputs(self.sources)
        self.param_inputs = inputs.get_param_inputs(self.sources)
        self.internal_inputs = inputs.get_internal_inputs(tool, self.sources)

        self.prescript: list[str] = []
        self.script: list[str] = []

    def generate(self) -> Tuple[Optional[str], str]:
        """Generate the script content of a Nextflow process for Janis command line tool"""
        self.handle_cmdtool_preprocessing()
        self.handle_cmdtool_base_command()
        self.handle_cmdtool_inputs()
        prescript = self.finalise_prescript()
        script = self.finalise_script()
        return prescript, script
    
    def handle_cmdtool_inputs(self) -> None:
        for inp in ordering.order_cmdtool_inputs_arguments(self.tool):
            if isinstance(inp, ToolInput):
                prescript, script = ScriptFormatter(
                    tinput=inp, 
                    tool=self.tool, 
                    process_inputs=self.process_inputs, 
                    param_inputs=self.param_inputs, 
                    internal_inputs=self.internal_inputs, 
                    sources=self.sources
                ).format()
                if prescript:
                    self.prescript.append(prescript)
                if script:
                    self.script.append(script)
            else:
                self.handle_tool_argument(inp)

    def handle_tool_argument(self, arg: ToolArgument) -> None:
        prefix = ''
        space = ''
        if arg.prefix:
            prefix = arg.prefix
            if arg.separate_value_from_prefix != False:
                space = ' '

        expr = unwrap_expression(
            val=arg.value,
            tool=self.tool,
            sources=self.sources,
            process_inputs=self.process_inputs,
            param_inputs=self.param_inputs,
            internal_inputs=self.internal_inputs,
            in_shell_script=True,
        )
        line = f'{prefix}{space}{expr}'
        self.script.append(line)

    def handle_cmdtool_preprocessing(self) -> None:
        for dirpath in self.tool.directories_to_create() or []:
            unwrapped_dir = unwrap_expression(
                val=dirpath, 
                tool=self.tool, 
                sources=self.sources,
                process_inputs=self.process_inputs,
                param_inputs=self.param_inputs,
                internal_inputs=self.internal_inputs,
                in_shell_script=True
            ) 
            line = f"mkdir -p '{unwrapped_dir}'"
            self.script.append(line)

    def handle_cmdtool_base_command(self) -> None:
        bc = self.tool.base_command()
        if bc is not None:
            if isinstance(bc, list):
                self.script += [' '.join([str(cmd) for cmd in bc])]
            else:
                self.script += [str(bc)]

    def finalise_prescript(self) -> Optional[str]:
        if self.prescript:
            return '\n'.join(self.prescript)
        return None

    def finalise_script(self) -> str:
        script = self.script
        if settings.JANIS_ASSISTANT:
            script = script + [f'| tee {self.stdout_filename}_{self.process_name}']
        script = [f'{ln} \\' for ln in script]
        return '\n'.join(script)


