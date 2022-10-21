

from typing import Any
from janis_core import CommandTool, ToolArgument, ToolInput
from .unwrap import unwrap_expression



def gen_script_for_cmdtool(
    tool: CommandTool, 
    inputs: list[ToolInput], 
    input_in_selectors: dict[str, Any],
    process_name: str,
    stdout_filename: str
) -> str:
    return ProcessScriptGenerator(
        tool,
        inputs,
        input_in_selectors,
        process_name,
        stdout_filename
    ).generate()


class ProcessScriptGenerator:
    def __init__(
        self,
        tool: CommandTool, 
        inputs: list[ToolInput], 
        input_in_selectors: dict[str, Any],
        process_name: str, 
        stdout_filename: str
    ):
        assert(tool)
        self.tool = tool
        self.inputs = inputs        # exposed_inputs vs internal_inputs???
        self.input_in_selectors = input_in_selectors
        self.process_name = process_name
        self.stdout_filename = stdout_filename

    def generate(self) -> str:
        """Generate the script content of a Nextflow process for Janis command line tool"""
        lines = []
        lines += self.gen_cmdtool_preprocessing()
        lines += self.gen_cmdtool_base_command()
        lines += self.gen_cmdtool_args()
        lines = [f'{ln} \\' for ln in lines]
        lines += [f'| tee {self.stdout_filename}_{self.process_name}']
        return '\n'.join(lines)

    def gen_cmdtool_preprocessing(self) -> list[str]:
        lines: list[str] = []
        for dirpath in self.tool.directories_to_create() or []:
            unwrapped_dir = unwrap_expression(
                value=dirpath, 
                input_in_selectors=self.input_in_selectors,
                inputs_dict=self.tool.inputs_map(),
                tool=self.tool, 
                in_shell_script=True
            ) 
            line = f"mkdir -p '{unwrapped_dir}'"
            lines.append(line)
        return lines

    def gen_cmdtool_base_command(self) -> list[str]:
        bc = self.tool.base_command()
        if bc is None:
            lines = []
        elif bc and isinstance(bc, list):
            lines = [' '.join([str(cmd) for cmd in bc])]
        else:
            lines = [str(bc)]
        return lines

    def gen_cmdtool_args(self) -> list[str]:
        lines: list[str] = []
        
        for inp in self.get_ordered_cmdtool_arguments():
            if isinstance(inp, ToolInput):
                lines.append(f"${inp.id()}WithPrefix")
            
            elif isinstance(inp, ToolArgument):
                expression = unwrap_expression(
                    value=inp.value,
                    input_in_selectors=self.input_in_selectors,
                    tool=self.tool,
                    inputs_dict=self.tool.inputs_map(),
                    skip_inputs_lookup=True,
                    quote_string=False,
                    in_shell_script=True,
                )
                
                if inp.prefix is not None:
                    space = ""
                    if inp.separate_value_from_prefix is not False:
                        space = " "
                    line = f'{inp.prefix}{space}"{expression}"'
                else:
                    line = expression

                lines.append(line)
            else:
                raise Exception("unknown input type")
        return lines

    def get_ordered_cmdtool_arguments(self) -> list[ToolInput | ToolArgument]:
        arguments = self.tool.arguments() or []
        args = [a for a in arguments if a.position is not None or a.prefix is not None]
        args += [a for a in self.inputs if a.position is not None or a.prefix is not None]
        args = sorted(args, key=lambda a: (a.prefix is None))
        args = sorted(args, key=lambda a: (a.position or 0))
        return args

