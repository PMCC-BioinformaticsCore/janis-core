

from janis_core.types import DataType

from .. import nfgen_utils


# helper functions to identify types involved in plumbing situations

def single_type(dtype: DataType) -> bool:
    if not dtype.is_array():
        if not secondary_type(dtype):
            return True
    return False

def array_type(dtype: DataType) -> bool:
    if not secondary_type(dtype) and not secondary_array_type(dtype):
        if dtype.is_array():
            return True
    return False

def secondary_type(dtype: DataType) -> bool:
    if nfgen_utils.is_secondary_type(dtype) or nfgen_utils.is_file_pair_type(dtype):
        return True
    return False

def secondary_array_type(dtype: DataType) -> bool:
    if nfgen_utils.is_array_secondary_type(dtype) or nfgen_utils.is_array_file_pair_type(dtype):
        return True
    return False