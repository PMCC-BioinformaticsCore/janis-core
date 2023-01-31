
from typing import Any
from janis_core import DataType, CommandTool, PythonTool, StepOutputSelector, String

from .. import nfgen_utils



def satisfies_edge_case(
    src: Any, 
    desttype: DataType, 
    tool: CommandTool | PythonTool
    ) -> bool:
        if is_pythontool_array_string_output(src):
            return True
        # add more plumbing edge cases here as they arise
        return False

def is_pythontool_array_string_output(src: Any) -> bool:
    source = src.source_map[0].source
    # source is from step output
    if isinstance(source, StepOutputSelector):
        tool = source.node.tool
        output = tool.outputs_map()[source.tag]
        dtype = output.outtype
        basetype = nfgen_utils.get_base_type(dtype)

        # source tool is pythontool
        if isinstance(tool, PythonTool):
            # source datatype is Array(String())
            if dtype.is_array() and isinstance(basetype, String):
                return True
    return False

def handle_edge_case(src: Any, desttype: DataType, tool: CommandTool | PythonTool) -> str:
    if is_pythontool_array_string_output(src):
        operation = ".filter{ it != '' }.map{ it -> it.split(', ') }.ifEmpty( null )"
    else:
        raise RuntimeError('DEV: there should be an edge case handled here, but apparently not')
    return operation