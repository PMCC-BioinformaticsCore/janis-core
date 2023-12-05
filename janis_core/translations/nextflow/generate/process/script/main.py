from typing import Tuple, Optional
import regex as re

from janis_core import CommandToolBuilder, Stdout, StringFormatter
from ....unwrap import unwrap_expression
from ....variables import VariableManager
from .prescript import gen_prescript_lines
from .script import gen_script_lines


def gen_nf_process_script(
    tool: CommandToolBuilder,
    vmanager: VariableManager,
) -> Tuple[Optional[str], str]:
    generator = ProcessScriptGenerator(tool, vmanager)
    return generator.generate()


class ProcessScriptGenerator:
    def __init__(
        self,
        tool: CommandToolBuilder, 
        vmanager: VariableManager
    ):
        self.tool = tool
        self.vmanager = vmanager
        self.prescript_lines: list[str] = []
        self.script_lines: list[str] = []

    def generate(self) -> Tuple[Optional[str], str]:
        """Generate the script content of a Nextflow process for Janis command line tool"""
        prescript = self.generate_prescript()
        
        if self.tool.is_shell_script:
            script = self.generate_script_shell()
        else:
            script = self.generate_script_default()
        
        return prescript, script
    
    def generate_prescript(self) -> Optional[str]:
        self.prescript_lines += gen_prescript_lines(self.tool, self.vmanager)
        if self.prescript_lines:
            return '\n'.join(self.prescript_lines)
        return None
    
    def generate_script_shell(self) -> str:
        assert 'script.sh' in self.tool._files_to_create
        contents = self.tool._files_to_create['script.sh']
        if not isinstance(contents, StringFormatter):
            raise NotImplementedError
        script = unwrap_expression(
            val=contents, 
            context='process_script',
            variable_manager=self.vmanager,
            tool=self.tool,
        ) 
        assert isinstance(script, str)
        script = self.deindent(script)
        return script
    
    def deindent(self, text: str) -> str:
        #wrapped in try except cuz seems dodgy
        try:
            # remove common indent
            PATTERN = r'(\n[ \t]+)'
            match = list(re.finditer(PATTERN, text))[0]
            indent = match.group(1)
            text = '\n' + text
            text = text.replace(indent, '\n')

            # remove beinning newlines
            text = re.sub(r'^\n*\n', '\n', text)
        except Exception as e:
            pass
        return text
    
    def generate_script_default(self) -> str:
        self.handle_cmdtool_directories()
        self.handle_cmdtool_base_command()
        self.handle_cmdtool_inputs_arguments()
        self.handle_stdout_redirect()
        return self.finalise_script()

    def handle_cmdtool_directories(self) -> None:
        """generate cmdline 'mkdir' statement for each directory to create"""
        for dirpath in self.tool.directories_to_create() or []:
            unwrapped_dir = unwrap_expression(
                val=dirpath, 
                context='process_script',
                variable_manager=self.vmanager,
                tool=self.tool,
            ) 
            line = f"mkdir -p '{unwrapped_dir}';"
            self.script_lines.append(line)

    def handle_cmdtool_base_command(self) -> None:
        bc = self.tool.base_command()
        if bc is not None:
            if isinstance(bc, list):
                self.script_lines += [' '.join([str(cmd) for cmd in bc])]
            else:
                self.script_lines += [str(bc)]

    def handle_cmdtool_inputs_arguments(self) -> None:
        self.script_lines += gen_script_lines(self.tool, self.vmanager)

    def handle_stdout_redirect(self) -> None:
        assert(len([out for out in self.tool.outputs() if isinstance(out.output_type, Stdout)]) <= 1)
        for out in self.tool.outputs():
            if isinstance(out.output_type, Stdout):
                if hasattr(out.output_type.subtype, 'extension') and out.output_type.subtype.extension is not None:
                    suffix = out.output_type.subtype.extension
                else:
                    suffix = ''
                self.script_lines.append(f'> {out.id()}{suffix}')

    def finalise_script(self) -> str: 
        script = self.script_lines
        if len(script) == 0:
            return ''
        elif len(script) == 1:
            return script[0]
        else:
            script = [f"{ln} \\" for ln in script[:-1]] + [script[-1]]  # type: ignore
        return '\n'.join(script) # type: ignore


        
                
    # def get_undefined_variable_references(self) -> set[str]:
    #     """
    #     ensure all referenced variables are defined
    #     references to process / param tool inputs are ok as always have variable defined.
    #     references to internal inputs are not defined. 

    #     for tool inputs / arguments / outputs referencing internal input variables:
    #       - can autofill their value if a default is present, or
    #       - must define as null ('def [var] = null') if no default.
    #     """
    #     undef_variables: set[str] = set()

    #     all_entities = self.tool.inputs() + self.tool.outputs()
    #     if self.tool.arguments():
    #         arguments: list[ToolArgument] = self.tool.arguments()
    #         all_entities += arguments
        
    #     for entity in all_entities:
    #         # get all referenced tool inputs
    #         referenced_ids = trace.trace_referenced_variables(entity, self.tool)
            
    #         # check if any are internal inputs with no default value
    #         for ref in referenced_ids:
    #             var = self.vmanager.get(ref).original
    #             if var.vtype == VariableType.IGNORED:
    #                 tinput = [x for x in self.tool.inputs() if x.id() == ref][0]
    #                 if tinput.default is None:
    #                     varname = naming.process.generic(tinput)
    #                     undef_variables.add(varname)
    #                     self.vmanager.update(
    #                         tinput_id=tinput.id(), 
    #                         vtype_str='local',
    #                         value=varname
    #                     )
        
    #     return undef_variables


    # def handle_undefined_variable_references(self) -> None:
    #     """
    #     create definitions for referenced tool inputs in pre-script section
    #     ensures all referenced variables in script are defined
    #     """
    #     undef_variables = self.get_undefined_variable_references()
    #     for varname in undef_variables:
    #         line = f'def {varname} = null'
    #         self.prescript_lines.append(line)


# handle_cmdtool_directories
# handle_cmdtool_base_command
# handle_cmdtool_inputs_arguments
# handle_stdout_redirect
# finalise_script