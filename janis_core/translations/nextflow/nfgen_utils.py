

import ast
from typing import Any, Optional

from janis_core.types import (
    Boolean,
    Int,
    Float,
)
from janis_core import (
    DataType,
    ToolInput,
    TInput,
)

from janis_core.workflow.workflow import InputNode
from janis_core import translation_utils as utils


"""
FOR THIS SECTION

In the current approach, only Files are supplied using channels. 
Process inputs can therefore only be File types. 
All other inputs are supplied a value using params. 

Additionally, when doing MINIMAL_PROCESS tool translation:
    - inputs with a default are autofilled
    - inputs which are optional are ignored

If doing workflow translation, we get extra info about process inputs
using step inputs. 

For example in MINIMAL_PROCESS mode, some optional files may be being passed 
values in the workflow. these would usually be ignored due to their optionality, but here 
shoud should be included.
"""


### GENERAL

def items_with_id(the_list: list[Any], ids: set[str]) -> list[Any]:
    return [x for x in the_list if x.id() in ids]

def to_groovy(
        val: Any, 
        dtype: Optional[DataType]=None, 
        quote_override: Optional[bool]=None,
        delim: Optional[str]=None,
    ) -> Any:
    # must work with str version. 
    if dtype is not None:
        dtype = utils.get_base_type(dtype)
    val = str(val)
    
    # # secondary files
    # if val != 'None' and isinstance(dtype, File) and dtype.has_secondary_files():
    #     primary_file = wrap(val)
    #     secondary_files: list[str] = []
    #     for suffix in dtype.secondary_files():
    #         sec_file = apply_secondary_file_format_to_filename(primary_file, suffix)
    #         sec_file = wrap(sec_file)
    #         secondary_files.append(sec_file)
    #     # Note: we want primary file to always be the first item in the array
    #     val = [primary_file] + secondary_files

    # wrap in quotes if necessary (for actual string values, file paths etc)
    if _should_wrap(val, dtype, quote_override):
        val = _wrap(val)

    # remove dollar variable references (unsure if needed)
    if '$' in val:
        val = val.replace('$', '')

    # 'None' -> 'null' etc
    val = _cast_keywords(val)
    return val

def _should_wrap(val: str, dtype: Optional[DataType], quote_override: Optional[bool]) -> bool:
    if quote_override is not None:
        return True if quote_override else False
    
    # don't quote lists
    try:
        literal_val = ast.literal_eval(val)
        if isinstance(literal_val, list):
            return False
    except Exception as e: 
        pass

    # don't quote None
    if val == 'None':
        return False

    # don't quote outer array, boolean, numeric types
    no_quote_types: list[type[DataType]] = [Boolean, Int, Float]
    if dtype:
        if type(dtype) in no_quote_types:
            return False

    # don't quote nextflow input channel
    if val.startswith('ch_'):
        return False
    
    # don't quote nextflow referenced variable
    if val.startswith('$'):
        return False
    
    # quote everything else
    return True

def _wrap(val: Any) -> Any:
    return f"'{val}'"

def _cast_keywords(val: str) -> str:
    # this is done in string world - need a better way of handling lists!
    keyword_map: dict[str, str] = {
        'None': 'null',
        'False': 'false',
        'True': 'true',
    }
    for python_val, groovy_val in keyword_map.items():
        if python_val in val:
            val = val.replace(python_val, groovy_val)
    return val


### TYPES 

def get_base_type_task_input(task_input: ToolInput | InputNode | TInput) -> DataType:
    match task_input:
        case ToolInput():
            dtype = task_input.input_type
        case InputNode():
            dtype = task_input.datatype
        case TInput():
            dtype = task_input.intype
        case _:
            raise NotImplementedError
    return utils.get_base_type(dtype)


