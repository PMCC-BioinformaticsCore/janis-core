
from typing import Any, Tuple, Optional

from janis_core import (
    CommandTool,
    ToolArgument,
    ToolInput,
)

from .... import naming
from .... import ordering

from ....unwrap import unwrap_expression
from .... import trace

from ....variables import VariableManager
from ....variables import VariableType
from .ScriptFormatter import ScriptFormatter


def gen_nf_process_script(
    tool: CommandTool,
    variable_manager: VariableManager,
    sources: dict[str, Any],
    stdout_filename: str,
) -> Tuple[Optional[str], str]:
    return ProcessScriptGenerator(
        tool=tool,
        variable_manager=variable_manager,
        stdout_filename=stdout_filename,
        sources=sources,
    ).generate()


class ProcessScriptGenerator:
    def __init__(
        self,
        tool: CommandTool, 
        variable_manager: VariableManager,
        stdout_filename: str,
        sources: Optional[dict[str, Any]]=None,
    ):
        self.tool = tool
        self.variable_manager = variable_manager
        self.stdout_filename = stdout_filename
        self.sources = sources if sources is not None else {}

        self.prescript: list[str] = []
        self.script: list[str] = []

    def generate(self) -> Tuple[Optional[str], str]:
        """Generate the script content of a Nextflow process for Janis command line tool"""
        self.handle_undefined_variable_references()
        self.handle_cmdtool_directories()
        self.handle_cmdtool_base_command()
        self.handle_cmdtool_inputs()
        prescript = self.finalise_prescript()
        script = self.finalise_script()
        return prescript, script
    
    def handle_cmdtool_inputs(self) -> None:
        tool_input_formatter = ScriptFormatter(
            tool=self.tool, 
            variable_manager=self.variable_manager,
            sources=self.sources
        )
        for inp in ordering.order_cmdtool_inputs_arguments(self.tool):
            if isinstance(inp, ToolInput):
                prescript, script = tool_input_formatter.format(inp)
                self.prescript += prescript
                self.script += script
            else:
                self.handle_tool_argument(inp)

    def handle_tool_argument(self, arg: ToolArgument) -> None:
        prefix = ''
        space = ''
        if arg.prefix:
            prefix = arg.prefix
            if arg.separate_value_from_prefix != False:
                space = ' '

        # unwrap toolargument value
        expr = unwrap_expression(
            val=arg.value,
            context='process_script',
            variable_manager=self.variable_manager,
            tool=self.tool,
            in_shell_script=True,
            quote_strings=arg.shell_quote
        )
        line = f'{prefix}{space}{expr}'
        self.script.append(line)

    def handle_undefined_variable_references(self) -> None:
        """
        create definitions for referenced tool inputs in pre-script section
        ensures all referenced variables in script are defined
        """
        undef_variables = self.get_undefined_variable_references()
        for varname in undef_variables:
            line = f'def {varname} = null'
            self.prescript.append(line)
                
    def get_undefined_variable_references(self) -> set[str]:
        """
        ensure all referenced variables are defined
        references to process / param tool inputs are ok as always have variable defined.
        references to internal inputs are not defined. 

        for tool inputs / arguments / outputs referencing internal input variables:
          - can autofill their value if a default is present, or
          - must define as null ('def [var] = null') if no default.
        """
        undef_variables: set[str] = set()

        all_entities = self.tool.inputs() + self.tool.outputs()
        if self.tool.arguments():
            arguments: list[ToolArgument] = self.tool.arguments()
            all_entities += arguments
        
        for entity in all_entities:
            # get all referenced tool inputs
            referenced_ids = trace.trace_referenced_variables(entity, self.tool)
            
            # check if any are internal inputs with no default value
            for ref in referenced_ids:
                var = self.variable_manager.get(ref).original
                if var.vtype == VariableType.IGNORED:
                    tinput = [x for x in self.tool.inputs() if x.id() == ref][0]
                    if tinput.default is None:
                        varname = naming.process.generic(tinput)
                        undef_variables.add(varname)
                        self.variable_manager.update(
                            tinput_id=tinput.id(), 
                            vtype_str='local',
                            value=varname
                        )
        
        return undef_variables

    def handle_cmdtool_directories(self) -> None:
        for dirpath in self.tool.directories_to_create() or []:
            unwrapped_dir = unwrap_expression(
                val=dirpath, 
                in_shell_script=True
            ) 
            line = f"mkdir -p '{unwrapped_dir}';"
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
        # if settings.JANIS_ASSISTANT:
        #     script = script + [f'| tee {self.stdout_filename}_{self.process_name}']
        
        script = self.script
        script = [f'{ln} \\' for ln in script]
        return '\n'.join(script)


