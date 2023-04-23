
from typing import Any

from janis_core import settings
from janis_core import CommandTool, PythonTool, Workflow
from janis_core.types import DataType, Array, Int, Float, Double, Boolean
NoneType = type(None)

from ... import naming
from ...model.process import NFProcess

from .directives import gen_nf_process_directives
from .inputs import gen_nf_process_inputs
from .script import gen_nf_process_script
from .outputs import gen_nf_process_outputs

from ...variables import init_variable_manager_for_task
from ...variables import VariableType


def generate_processes(wf: Workflow) -> dict[str, NFProcess]:
    """for each CommandTool | PythonTool in workflow, generate a nextflow process"""
    process_dict = {}
    return do_generate_processes(wf, process_dict)

def do_generate_processes(wf: Workflow, process_dict: dict[str, NFProcess]) -> dict[str, NFProcess]:
    for step in wf.step_nodes.values():
        
        # create process for each cmdtool / pytool if not already done
        if isinstance(step.tool, CommandTool) or isinstance(step.tool, PythonTool):
            tool_id = step.tool.id()
            if tool_id not in process_dict:
                process = generate_process(step.tool)
                process_dict[tool_id] = process
        
        # recursively do for subworkflows 
        elif isinstance(step.tool, Workflow):
            process_dict = do_generate_processes(step.tool, process_dict)
        
        else:
            raise RuntimeError
    
    return process_dict

def generate_process(tool: CommandTool | PythonTool) -> NFProcess:
    if isinstance(tool, CommandTool):
        return generate_process_cmdtool(tool)
    elif isinstance(tool, PythonTool):  # type: ignore
        return generate_process_pythontool(tool)
    else:
        raise RuntimeError

def generate_process_cmdtool(tool: CommandTool) -> NFProcess:
    generator = CmdToolProcessGenerator(tool)
    return generator.generate()

def generate_process_pythontool(tool: PythonTool) -> NFProcess:
    generator = PythonToolProcessGenerator(tool)
    return generator.generate()



# helper class
class CmdToolProcessGenerator:
    def __init__(self, tool: CommandTool) -> None:
        self.tool = tool
        self.vmanager = init_variable_manager_for_task(self.tool)

    @property
    def name(self) -> str:
        return naming.constructs.gen_varname_process(self.tool.id())

    def generate(self) -> NFProcess:
        # directives
        resources = {}
        process_directives = gen_nf_process_directives(self.tool, resources)

        # inputs
        process_inputs = gen_nf_process_inputs(self.tool)

        # script
        pre_script, main_script = gen_nf_process_script(
            tool=self.tool,
            variable_manager=self.vmanager,
            stdout_filename=settings.translate.nextflow.TOOL_STDOUT_FILENAME,
        )
        
        # outputs
        process_outputs = gen_nf_process_outputs(
            tool=self.tool, 
            variable_manager=self.vmanager,
        )

        # process
        process = NFProcess(
            name=self.name,
            pre_script=pre_script,
            script=main_script,
            inputs=process_inputs,
            outputs=process_outputs,
            directives=process_directives
        )

        return process


PYINDENT = '    '

# helper class
class PythonToolProcessGenerator:
    def __init__(self, tool: PythonTool) -> None:
        self.tool = tool
        self.vmanager = init_variable_manager_for_task(self.tool)

    @property
    def name(self) -> str:
        return naming.constructs.gen_varname_process(self.tool.id())
    
    @property
    def args(self) -> list[str]:
        # TODO: handle args of type list of string (need to quote them)
        args: list[str] = []
        
        # supply args to pythontool call
        for tinput in self.tool.tool_inputs():
            tag: str = tinput.tag
            value: Any = None
            dtype: DataType = tinput.intype

            # get final arg for this tinput 
            cvar = self.vmanager.get(tag).current
            if cvar.vtype in [VariableType.TASK_INPUT, VariableType.PARAM]:
                
                # get varname
                # workaround for secondaries
                if isinstance(cvar.value, list):
                    varname = cvar.value[0]
                else:
                    varname = cvar.value
                
                # format varname
                if isinstance(dtype, Array):
                    value = f'"${{{varname}}}".split(" ")'  # this is python not groovy
                else:
                    value = f'${{{varname}}}'

            elif cvar.vtype == VariableType.STATIC:
                value = cvar.value
            
            elif cvar.vtype == VariableType.IGNORED:
                continue
            
            elif cvar.vtype == VariableType.LOCAL:
                value = cvar.value
            
            else:
                raise RuntimeError(f"Unknown variable type: {cvar.vtype}")

            # wrap in quotes unless numeric or bool
            if not isinstance(dtype, (Array, Int, Float, Double, Boolean, NoneType)):
                value = f'"{value}"'

            arg = f"{tag}={value}"
            args.append(arg) 
        
        return args 
    
    @property
    def call_block(self) -> str:
        # one-liner if only one arg
        if len(self.args) <= 1:
            args_str = ", ".join(a for a in self.args)
            args_str = f'result = code_block({args_str})'
            return args_str
        
        # multi-liner if more than one arg
        else:
            args_str = ''
            args_str += f'result = code_block(\n'
            
            # middle lines get commas on the end
            for arg in self.args[:-1]:
                args_str += f'{PYINDENT}{arg},\n'
            
            # last line
            args_str += f'{PYINDENT}{self.args[-1]}\n'
            args_str += f')'
            return args_str


    @property
    def script(self) -> str:
        """
        Generate the content of the script section in a Nextflow process for Janis python code tool
        """

        script = f"""\
{settings.translate.nextflow.PYTHON_SHEBANG}

from ${{code_file}} import code_block
import os
import json

{self.call_block}

work_dir = os.getcwd()
for key in result:
    with open(os.path.join(work_dir, f"{settings.translate.nextflow.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{{key}}"), "w") as fp:
        fp.write(json.dumps(result[key]))
"""
        return script

    def generate(self) -> NFProcess:
        # directives
        resources = {}
        process_directives = gen_nf_process_directives(self.tool, resources)
        
        # inputs
        process_inputs = gen_nf_process_inputs(self.tool)

        # outputs
        process_outputs = gen_nf_process_outputs(
            tool=self.tool, 
            variable_manager=self.vmanager,
        )

        # process
        process = NFProcess(
            name=self.name,
            script=self.script,
            inputs=process_inputs,
            outputs=process_outputs,
            main_exec='',
            directives=process_directives
        )

        return process
