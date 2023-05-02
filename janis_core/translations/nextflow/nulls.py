

from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType
from janis_core.types import DataType, File
from janis_core import settings

NULL = settings.translate.nextflow.NULL


def get_null_value(dtype: DataType, as_param: bool=False, should_add_file_cast: bool=False) -> str:
    dtt = utils.get_dtt(dtype)

    if dtt == DTypeType.SECONDARY_ARRAY and dtype.optional:
        expr = secondary_array_optional_null(dtype)
    
    elif dtt == DTypeType.SECONDARY_ARRAY:
        expr = secondary_array_mandatory_null()
    
    elif dtt == DTypeType.SECONDARY and dtype.optional:
        expr = secondary_optional_null(dtype)
    
    elif dtt == DTypeType.SECONDARY:
        expr = secondary_mandatory_null()

    elif dtt == DTypeType.FILE_PAIR_ARRAY and dtype.optional:
        expr = file_pair_array_optional_null()

    elif dtt == DTypeType.FILE_PAIR_ARRAY:
        expr = file_pair_array_mandatory_null()

    elif dtt == DTypeType.FILE_PAIR and dtype.optional:
        expr = file_pair_optional_null()

    elif dtt == DTypeType.FILE_PAIR:
        expr = file_pair_mandatory_null()

    elif dtt == DTypeType.FILE_ARRAY and dtype.optional:
        expr = file_array_optional_null()

    elif dtt == DTypeType.FILE_ARRAY:
        expr = file_array_mandatory_null()

    elif dtt == DTypeType.FILE and dtype.optional:
        expr = file_optional_null()

    elif dtt == DTypeType.FILE:
        expr = file_mandatory_null()

    elif dtt == DTypeType.GENERIC_ARRAY and dtype.optional:
        expr = generic_array_optional_null()

    elif dtt == DTypeType.GENERIC_ARRAY:
        expr = generic_array_mandatory_null()
    
    elif dtt == DTypeType.GENERIC and dtype.optional:
        expr = generic_optional_null()

    elif dtt == DTypeType.GENERIC:
        expr = generic_mandatory_null()
    
    elif dtt == DTypeType.FLAG_ARRAY and dtype.optional:
        expr = generic_optional_null()

    elif dtt == DTypeType.FLAG_ARRAY:
        expr = generic_optional_null()
    
    elif dtt == DTypeType.FLAG and dtype.optional:
        expr = generic_optional_null()

    elif dtt == DTypeType.FLAG:
        expr = generic_optional_null()

    else:
        raise RuntimeError(f"Unknown datatype: {dtt}")
    
    # present as reference to NULL param if requested
    if as_param:
        expr = expr.replace(NULL, f'params.{NULL}')

    if should_add_file_cast:
        expr = add_file_cast(dtype, expr)
    
    return expr


def add_file_cast(dtype: DataType, expr: str) -> str:
    dtt = utils.get_dtt(dtype)
    if dtt == DTypeType.SECONDARY_ARRAY:
        expr = f'{expr}.collect{{ it.collect{{ file(it) }} }}'
    
    elif dtt == DTypeType.SECONDARY:
        expr = f'{expr}.collect{{ file(it) }}'
    
    elif dtt == DTypeType.FILE_PAIR_ARRAY:
        expr = f'{expr}.collect{{ it.collect{{ file(it) }} }}'
    
    elif dtt == DTypeType.FILE_PAIR:
        expr = f'{expr}.collect{{ file(it) }}'
    
    elif dtt == DTypeType.FILE_ARRAY:
        expr = f'{expr}.collect{{ file(it) }}'
    
    elif dtt == DTypeType.FILE:
        expr = f'file( {expr} )'
    
    return expr


# NON-OPTIONAL
def secondary_array_mandatory_null() -> str:
    return '[[]]'

def secondary_mandatory_null() -> str:
    return '[]'

def file_pair_array_mandatory_null() -> str:
    return '[[]]'

def file_pair_mandatory_null() -> str:
    return '[]'

def file_array_mandatory_null() -> str:
    return '[]'

def file_mandatory_null() -> str:
    return NULL

def generic_array_mandatory_null() -> str:
    return '[]'

def generic_mandatory_null() -> str:
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
    return NULL

def generic_array_optional_null() -> str:
    return NULL

def generic_optional_null() -> str:
    return NULL

