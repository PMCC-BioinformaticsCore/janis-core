

from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType
from janis_core.types import DataType, File
from janis_core import settings

NULL_VAL = settings.translate.nextflow.NULL_VALUE


def get_null_value(dtype: DataType, as_param: bool=False, should_add_file_cast: bool=False) -> str:
    dtt = utils.get_dtt(dtype)

    if dtt == DTypeType.SECONDARY_ARRAY:
        expr = '[[]]'
    
    elif dtt == DTypeType.SECONDARY:
        expr = '[]'

    elif dtt == DTypeType.FILE_PAIR_ARRAY:
        expr = '[[]]'

    elif dtt == DTypeType.FILE_PAIR:
        expr = '[]'

    elif dtt == DTypeType.FILE_ARRAY:
        expr = '[]'

    elif dtt == DTypeType.FILE:
        expr = NULL_VAL
    
    elif dtt == DTypeType.FILENAME:
        expr = NULL_VAL

    elif dtt == DTypeType.GENERIC_ARRAY:
        expr = NULL_VAL
    
    elif dtt == DTypeType.GENERIC:
        expr = NULL_VAL
    
    elif dtt == DTypeType.FLAG_ARRAY:
        expr = NULL_VAL
    
    elif dtt == DTypeType.FLAG:
        expr = 'false'

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

