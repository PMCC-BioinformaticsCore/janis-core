

from janis_core.types import DataType, Array

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

def get_collate_size(dtype: DataType) -> int:
    basetype = nfgen_utils.get_base_type(dtype)
    assert(basetype)
    if basetype.name() in ['FastqPair', 'FastqGzPair']:
        size = 2
    else:
        exts = nfgen_utils.get_extensions(basetype)
        size = len(exts)
    return size


# from janis_bioinformatics.data_types import FastqGzPair