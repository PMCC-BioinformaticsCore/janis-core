
from typing import Optional

from janis_core import CommandTool

from janis_core import translation_utils as utils

from ...model.files import NFImportsBlock
from ...model.files import NFFunctionsBlock

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