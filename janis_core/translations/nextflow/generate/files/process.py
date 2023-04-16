
from typing import Optional

from janis_core import CommandTool, PythonTool
from janis_core import translation_utils as utils

from ...model.files import NFFile
from ...model.files import NFImportsBlock
from ...model.files import NFFunctionsBlock
from ...model.process import NFProcess

from ...trace import trace_entity_counts


get_primary_files_code = """\
def get_primary_files(var, element_count) {
    def primary_files = []
    var.eachWithIndex {item, index -> 
        if (index % element_count == 0) {
            primary_files.add(item)
        }
    }
    return primary_files
}"""


def generate_file_pythontool(process: NFProcess, tool: PythonTool) -> NFFile:
    """generates nextflow file for nextflow process derived from PythonTool"""
    nf_file = NFFile(subtype='process', name=process.name)
    nf_file.items.append(process)
    return nf_file


def generate_file_cmdtool(process: NFProcess, tool: CommandTool) -> NFFile:
    """generates nextflow file for nextflow process derived from CommandTool"""
    nf_file = NFFile(subtype='process', name=process.name)

    # groovy library imports & groovy functions used in process
    # item: imports
    imports_item = gen_imports_for_process_file(tool)
    if imports_item:
        nf_file.items.append(imports_item)
    
    # item: functions
    functions_item = gen_functions_for_process_file(tool)
    if functions_item:
        nf_file.items.append(functions_item)

    # item: process
    nf_file.items.append(process)

    # process_item = cls.handle_container(scope, tool, process_item)
    # cls.item_register.add(scope, process_item)
    return nf_file


def gen_imports_for_process_file(tool: CommandTool) -> Optional[NFImportsBlock]:
    # methods: list[str] = []
    imports: list[str] = []
    declarations: list[str] = []

    if _should_add_json_slurper(tool):
        # methods.append('include')
        imports.append('import groovy.json.JsonSlurper')
        declarations.append('jsonSlurper = new JsonSlurper()')
    
    if imports:
        # return NFImportsBlock(methods, imports, declarations)
        return NFImportsBlock(imports, declarations)
    else:
        return None

def _should_add_json_slurper(tool: CommandTool) -> bool:
    for toutput in tool.outputs():
        entity_counts = trace_entity_counts(toutput.selector, tool)
        if 'ReadJsonOperator' in entity_counts:
            return True
    return False

def gen_functions_for_process_file(tool: CommandTool) -> Optional[NFFunctionsBlock]:
    funcs: list[str] = []

    if _should_add_get_primary_files(tool):
        funcs.append(get_primary_files_code)
    
    if funcs:
        return NFFunctionsBlock(funcs)
    else:
        return None
    
def _should_add_get_primary_files(tool: CommandTool) -> bool:
    for tinput in tool.inputs():
        if utils.is_array_secondary_type(tinput.input_type):
            return True
    return False