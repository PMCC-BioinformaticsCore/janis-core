
from typing import Any

from janis_core import settings
from janis_core import CommandTool, PythonTool, Workflow
from janis_core.types import DataType, Array, Int, Float, Double, Boolean
NoneType = type(None)

from ... import naming
from ...model.process import NFProcess
from ...model.process import NFProcessInput
from ...model.process import NFPathProcessInput

from .directives import gen_nf_process_directives
from .inputs import gen_nf_process_inputs
from .script import gen_nf_process_script
from .outputs import gen_nf_process_outputs

from ...variables import init_variable_manager_for_task


def generate_processes(wf: Workflow) -> dict[str, NFProcess]:
    process_dict = {}
    return do_generate_processes(wf, process_dict)

def do_generate_processes(wf: Workflow, process_dict: dict[str, NFProcess]) -> dict[str, NFProcess]:
    for step in wf.step_nodes.values():
        
        # create process for each cmdtool / pytool if not already done
        if isinstance(step.tool, CommandTool) or isinstance(step.tool, PythonTool):
            tool_id = step.tool.id()
            if tool_id not in process_dict:
                process = generate_process(step.tool, step.sources, wf)
                process_dict[tool_id] = process
        
        # recursively do for subworkflows 
        elif isinstance(step.tool, Workflow):
            process_dict = do_generate_processes(step.tool, process_dict)
        
        else:
            raise RuntimeError
    
    return process_dict

def generate_process(tool: CommandTool | PythonTool, sources: dict[str, Any]) -> NFProcess:
    if isinstance(tool, CommandTool):
        return generate_process_cmdtool(tool, sources)
    elif isinstance(tool, PythonTool):  # type: ignore
        return generate_process_pythontool(tool, sources)
    else:
        raise RuntimeError

def generate_process_cmdtool(tool: CommandTool, sources: dict[str, Any]) -> NFProcess:
    generator = CmdToolProcessGenerator(tool, sources)
    return generator.generate()

def generate_process_pythontool(tool: PythonTool, sources: dict[str, Any]) -> NFProcess:
    generator = PythonToolProcessGenerator(tool, sources)
    return generator.generate()



# helper class
class CmdToolProcessGenerator:
    def __init__(self, tool: CommandTool, sources: dict[str, Any]) -> None:
        self.tool = tool
        self.sources = sources

    @property
    def name(self) -> str:
        return naming.constructs.gen_varname_process(self.tool.id())

    def generate(self) -> NFProcess:
        # variable mappings 
        self.vmanager = init_variable_manager_for_task(self.tool)

        # directives
        resources = {}
        process_directives = gen_nf_process_directives(self.tool, resources)

        # inputs
        process_inputs = gen_nf_process_inputs(self.tool)

        # script
        pre_script, main_script = gen_nf_process_script(
            tool=self.tool,
            variable_manager=self.vmanager,
            sources=self.sources,
            stdout_filename=settings.translate.nextflow.TOOL_STDOUT_FILENAME,
        )
        
        # outputs
        process_outputs = gen_nf_process_outputs(
            tool=self.tool, 
            variable_manager=self.vmanager,
            sources=self.sources
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


# helper class
class PythonToolProcessGenerator:
    def __init__(self, tool: PythonTool, sources: dict[str, Any]) -> None:
        self.tool = tool
        self.sources = sources

    @property
    def name(self) -> str:
        return naming.constructs.gen_varname_process(self.tool.id())
    
    @property
    def script(self) -> str:
        """
        Generate the content of the script section in a Nextflow process for Janis python code tool

        :param tool:
        :type tool:
        :param inputs:
        :type inputs:
        :return:
        :rtype:
        """
        # TODO: handle args of type list of string (need to quote them)
        args: list[str] = []
        
        for tinput in self.tool.inputs():
            tag: str = tinput.tag
            value: Any = None
            dtype: DataType = tinput.intype

            if tinput.id() in data_sources.task_inputs(scope) or tinput.id() in data_sources.param_inputs(scope):
                varname = data_sources.get(scope, tinput).value
                if isinstance(varname, list):
                    varname = varname[0]
                value = f'${{{varname}}}'
                if isinstance(dtype, Array):
                    value = f'"{value}".split(" ")'

            elif tinput.default is not None:
                value = tinput.default

            elif tinput.intype.optional == True:
                value = None

            else:
                raise NotImplementedError

            # wrap in quotes unless numeric or bool
            if not isinstance(dtype, (Array, Int, Float, Double, Boolean, NoneType)):
                value = f'"{value}"'

            arg = f"{tag}={value}"
            args.append(arg)  

        args_str = ", ".join(a for a in args)
        script = f"""\
    {settings.translate.nextflow.PYTHON_SHEBANG}

    from ${{code_file.simpleName}} import code_block
    import os
    import json

    result = code_block({args_str})

    work_dir = os.getcwd()
    for key in result:
        with open(os.path.join(work_dir, f"{settings.translate.nextflow.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{{key}}"), "w") as fp:
            fp.write(json.dumps(result[key]))
    """
        return script

    def generate(self) -> NFProcess:
        # variable mappings 
        self.vmanager = init_variable_manager_for_task(self.tool)
        
        # directives
        resources = {}
        process_directives = gen_nf_process_directives(tool, resources)
        
        # inputs
        process_inputs: list[NFProcessInput] = []
        
        # inputs: python script
        python_file_input = NFPathProcessInput(name=settings.translate.nextflow.PYTHON_CODE_FILE_SYMBOL)
        process_inputs.append(python_file_input)

        # inputs: tool inputs
        process_inputs += gen_nf_process_inputs(scope, tool)

        # outputs
        process_outputs = gen_nf_process_outputs(
            scope=scope, 
            tool=tool, 
            variable_manager=variable_manager, 
            sources=sources
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
