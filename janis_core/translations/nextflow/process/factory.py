
from typing import Any, Optional

from janis_core import CommandTool, PythonTool
from janis_core.types import DataType, Array, Int, Float, Double, Boolean
NoneType = type(None)

from ..scope import Scope
from janis_core import settings
from .. import naming
from janis_core import translation_utils as utils

from . import model
from . import script
from . import directives
from . import inputs
from . import outputs
from ..plumbing import trace_entity_counts



get_primary_files_code = """\
def get_primary_files(var) {
    def primary_files = []
    var.eachWithIndex {item, index -> 
        if (index % 2 == 0) {
            primary_files.add(item)
        }
    }
    return primary_files
}"""


def gen_imports_for_process(tool: CommandTool) -> Optional[model.ImportsBlock]:
    imports: list[str] = []
    declarations: list[str] = []

    if should_add_json_slurper(tool):
        imports.append('groovy.json.JsonSlurper')
        declarations.append('jsonSlurper = new JsonSlurper()')
    
    if imports:
        return model.ImportsBlock(imports, declarations)
    else:
        return None

def should_add_json_slurper(tool: CommandTool) -> bool:
    for toutput in tool.outputs():
        entity_counts = trace_entity_counts(toutput.selector, tool)
        if 'ReadJsonOperator' in entity_counts:
            return True
    return False

def gen_functions_for_process(tool: CommandTool) -> Optional[model.FunctionsBlock]:
    funcs: list[str] = []

    if should_add_get_primary_files(tool):
        funcs.append(get_primary_files_code)
    
    if funcs:
        return model.FunctionsBlock(funcs)
    else:
        return None
    
def should_add_get_primary_files(tool: CommandTool) -> bool:
    for tinput in tool.inputs():
        if utils.is_array_secondary_type(tinput.input_type):
            return True
    return False


def gen_process_from_cmdtool(tool: CommandTool, sources: dict[str, Any], scope: Scope,
) -> model.Process:
    """
    Generate a Nextflow Process object for a Janis Command line tool

    :param tool:
    :type tool:
    :param name: Generally, this is a workflow step id, so that we can prefix variables or process names
    :type name:
    :param provided_inputs:
    :type provided_inputs:
    :return:
    :rtype:
    """
    # name
    process_name = scope.labels[-1]

    # directives
    resources = {}
    process_directives = directives.gen_directives_for_process(tool, resources, scope)

    # inputs
    process_inputs = inputs.create_nextflow_process_inputs(tool, sources)

    # outputs
    process_outputs = outputs.create_nextflow_process_outputs(tool, sources)

    # script
    pre_script, main_script = script.gen_script_for_cmdtool(
        tool=tool,
        scope=scope,
        sources=sources,
        stdout_filename=settings.translate.nextflow.TOOL_STDOUT_FILENAME,
    )
    
    # process
    process = model.Process(
        name=process_name,
        pre_script=pre_script,
        main_script=main_script,
        script_type=model.ProcessScriptType.script,
        inputs=process_inputs,
        outputs=process_outputs,
        directives=process_directives
    )

    return process



def gen_process_from_codetool(
    tool: PythonTool,
    sources: dict[str, Any],   # values fed to tool inputs (step translation)
    scope: Scope,
) -> model.Process:
    """
    Generate a Nextflow Process object for Janis python code tool

    :param tool:
    :type tool:
    :param name:
    :type name:
    :param provided_inputs:
    :type provided_inputs:
    :return:
    :rtype:
    """   

    # name
    process_name = scope.labels[-1] if scope.labels else tool.id()

    # directives
    resources = {}
    process_directives = directives.gen_directives_for_process(tool, resources, scope)
    
    # inputs
    process_inputs: list[inputs.ProcessInput] = []
    
    # inputs: python script
    python_file_input = inputs.PathProcessInput(name=settings.translate.nextflow.PYTHON_CODE_FILE_SYMBOL)
    process_inputs.append(python_file_input)

    # inputs: tool inputs
    process_inputs += inputs.create_nextflow_process_inputs(tool, sources)

    # outputs
    process_outputs = outputs.create_nextflow_process_outputs(tool, sources)

    # script
    main_script = prepare_script_for_python_code_tool(tool, sources=sources)

    # process
    process = model.Process(
        name=process_name,
        main_script=main_script,
        script_type=model.ProcessScriptType.script,
        inputs=process_inputs,
        outputs=process_outputs,
        main_exec='',
        directives=process_directives
    )

    return process



def prepare_script_for_python_code_tool(tool: PythonTool, sources: dict[str, Any]) -> str:
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
    process_inputs = inputs.get_process_inputs(sources)
    param_inputs = inputs.get_param_inputs(sources)
    
    for inp in tool.inputs():
        tag: str = inp.tag
        value: Any = None
        dtype: DataType = inp.intype

        if inp.id() in process_inputs or inp.id() in param_inputs:
            src = naming.get_varname_toolinput(inp, process_inputs, param_inputs, sources)
            value = f'${{{src}}}'
            if isinstance(dtype, Array):
                value = f'"{value}".split(" ")'

        elif inp.default is not None:
            value = inp.default

        elif inp.intype.optional == True:
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
