

from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType
from janis_core.types import DataType, File
from janis_core import settings

NULL_VAL = settings.translate.nextflow.NULL_VALUE


def get_null_value(dtype: DataType, as_param: bool=False, should_add_file_cast: bool=False) -> str:
    dtt = utils.get_dtt(dtype)

    if dtt == DTypeType.SECONDARY_ARRAY:
        expr = secondary_array_null()
    
    elif dtt == DTypeType.SECONDARY:
        expr = secondary_null()

    elif dtt == DTypeType.FILE_PAIR_ARRAY:
        expr = file_pair_array_null()

    elif dtt == DTypeType.FILE_PAIR:
        expr = file_pair_null()

    elif dtt == DTypeType.FILE_ARRAY:
        expr = file_array_null()

    elif dtt == DTypeType.FILE:
        expr = file_null()

    elif dtt == DTypeType.GENERIC_ARRAY:
        expr = generic_array_null()
    
    elif dtt == DTypeType.GENERIC:
        expr = generic_null()
    
    elif dtt == DTypeType.FLAG_ARRAY:
        expr = generic_null()
    
    elif dtt == DTypeType.FLAG:
        expr = generic_null()

    else:
        raise RuntimeError(f"Unknown datatype: {dtt}")
    
    # present as reference to NULL param if requested
    if as_param:
        expr = expr.replace(NULL_VAL, f'params.{NULL_VAL}')

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
    return NULL_VAL

def generic_array_null() -> str:
    return NULL_VAL

def generic_null() -> str:
    return NULL_VAL


# # OPTIONAL
# def secondary_array_optional_null(dtype: DataType) -> str:
#     basetype = utils.get_base_type(dtype)
#     assert(isinstance(basetype, File))
#     exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
#     val = [f"{NULL}"] * len(exts)
#     val = ', '.join(val)
#     val = f"[[{val}]]"
#     return val

# def secondary_optional_null(dtype: DataType) -> str:
#     assert(isinstance(dtype, File))
#     exts = utils.get_extensions(dtype, remove_prefix_symbols=True)
#     val = [f"{NULL}"] * len(exts)
#     val = ', '.join(val)
#     val = f"[{val}]"
#     return val

# def file_pair_array_optional_null() -> str:
#     return f"[[{NULL}, {NULL}]]"

# def file_pair_optional_null() -> str:
#     return f"[{NULL}, {NULL}]"

# def file_array_optional_null() -> str:
#     return f"[{NULL}]"

# def file_optional_null() -> str:
#     return NULL

# def generic_array_optional_null() -> str:
#     return NULL

# def generic_optional_null() -> str:
#     return NULL

