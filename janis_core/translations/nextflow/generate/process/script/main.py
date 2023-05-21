
from typing import Tuple, Optional

from janis_core import CommandTool, ToolArgument
from janis_core import settings

from .... import naming
from .....common import trace

from ....unwrap import unwrap_expression
from ....variables import VariableManager
from ....variables import VariableType

from .prescript import gen_prescript_lines
from .script import gen_script_lines


def gen_nf_process_script(
    tool: CommandTool,
    vmanager: VariableManager,
) -> Tuple[Optional[str], str]:
    return ProcessScriptGenerator(
        tool=tool,
        vmanager=vmanager,
    ).generate()


class ProcessScriptGenerator:
    def __init__(
        self,
        tool: CommandTool, 
        vmanager: VariableManager
    ):
        self.tool = tool
        self.vmanager = vmanager

        self.prescript: list[str] = []
        self.script: list[str] = []

    def generate(self) -> Tuple[Optional[str], str]:
        """Generate the script content of a Nextflow process for Janis command line tool"""
        self.handle_undefined_variable_references()
        self.handle_cmdtool_directories()
        self.handle_cmdtool_base_command()
        if settings.translate.MODE != 'skeleton':
            self.handle_cmdtool_inputs_arguments()
        prescript = self.finalise_prescript()
        script = self.finalise_script()
        return prescript, script
    
    def handle_cmdtool_inputs_arguments(self) -> None:
        self.prescript += gen_prescript_lines(self.tool, self.vmanager)
        self.script += gen_script_lines(self.tool, self.vmanager)
        
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
                var = self.vmanager.get(ref).original
                if var.vtype == VariableType.IGNORED:
                    tinput = [x for x in self.tool.inputs() if x.id() == ref][0]
                    if tinput.default is None:
                        varname = naming.process.generic(tinput)
                        undef_variables.add(varname)
                        self.vmanager.update(
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
        script = self.script
        script = [f'{ln} \\' for ln in script]
        return '\n'.join(script)

