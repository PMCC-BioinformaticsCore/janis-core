

from janis_core import translation_utils as utils
from janis_core.types import DataType, File
from janis_core import settings

NULL = settings.translate.nextflow.NULL


def get_null_value(dtype: DataType, as_param: bool=False) -> str:
    if utils.is_secondary_array_type(dtype) and dtype.optional:
        expr = secondary_array_optional_null(dtype)
    
    elif utils.is_secondary_array_type(dtype):
        expr = secondary_array_null()
    
    elif utils.is_secondary_type(dtype) and dtype.optional:
        expr = secondary_optional_null(dtype)
    
    elif utils.is_secondary_type(dtype):
        expr = secondary_null()

    elif utils.is_file_pair_array_type(dtype) and dtype.optional:
        expr = file_pair_array_optional_null()

    elif utils.is_file_pair_array_type(dtype):
        expr = file_pair_array_null()

    elif utils.is_file_pair_type(dtype) and dtype.optional:
        expr = file_pair_optional_null()

    elif utils.is_file_pair_type(dtype):
        expr = file_pair_null()

    elif utils.is_file_array_type(dtype) and dtype.optional:
        expr = file_array_optional_null()

    elif utils.is_file_array_type(dtype):
        expr = file_array_null()

    elif utils.is_file_type(dtype) and dtype.optional:
        expr = file_optional_null()

    elif utils.is_file_type(dtype):
        expr = file_null()

    elif utils.is_array_type(dtype):
        expr = array_null()

    else:
        expr = nonfile_null()
    
    # present as reference to NULL param if requested
    if as_param:
        expr = expr.replace(NULL, f'params.{NULL}')
    
    return expr


# NON-OPTIONAL
def secondary_array_null() -> str:
    return '[[]]'

def secondary_null() -> str:
    return '[]'

def file_pair_array_null() -> str:
    return '[[]]'

def file_pair_null() -> str:
    return '[]'

def file_array_null() -> str:
    return '[]'

def file_null() -> str:
    return NULL

def array_null() -> str:
    return '[]'

def nonfile_null() -> str:
    return NULL


# OPTIONAL
def secondary_array_optional_null(dtype: DataType) -> str:
    basetype = utils.get_base_type(dtype)
    assert(isinstance(basetype, File))
    exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
    val = [f"{NULL}"] * len(exts)
    val = ', '.join(val)
    val = f"[[{val}]]"
    return val

def secondary_optional_null(dtype: DataType) -> str:
    assert(isinstance(dtype, File))
    exts = utils.get_extensions(dtype, remove_prefix_symbols=True)
    val = [f"{NULL}"] * len(exts)
    val = ', '.join(val)
    val = f"[{val}]"
    return val

def file_pair_array_optional_null() -> str:
    return f"[[{NULL}, {NULL}]]"

def file_pair_optional_null() -> str:
    return f"[{NULL}, {NULL}]"

def file_array_optional_null() -> str:
    return f"[{NULL}]"

def file_optional_null() -> str:
    return f"'{NULL}'"

