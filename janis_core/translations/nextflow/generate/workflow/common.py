

from typing import Optional, Type

from janis_core.types import DataType, UnionType

from janis_core import translation_utils as utils


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
    if utils.is_secondary_type(dtype) or utils.is_file_pair_type(dtype):
        return True
    return False

def secondary_array_type(dtype: DataType) -> bool:
    if utils.is_array_secondary_type(dtype) or utils.is_array_file_pair_type(dtype):
        return True
    return False

def get_collate_size(dtype: DataType) -> int:
    basetype = utils.get_base_type(dtype)
    assert(basetype)
    if basetype.name() in ['FastqPair', 'FastqGzPair']:
        size = 2
    else:
        exts = utils.get_extensions(basetype)
        size = len(exts)
    return size


# UnionType stuff
def union_type(dtype: DataType) -> bool:
    basetype = utils.get_base_type(dtype)
    if isinstance(basetype, UnionType):
        return True
    return False

def get_common_type(srctype: DataType, desttype: DataType) -> Optional[DataType]:
    srctype_base = utils.get_base_type(srctype)
    desttype_base = utils.get_base_type(desttype)

    # both are non-union types, & types match
    if not union_type(srctype_base) and not union_type(desttype_base):
        if srctype_base.name() == desttype_base.name():
            return srctype_base

    # at least one is union type, & common type exists
    elif union_type(srctype_base) or union_type(desttype_base):
        common_types = _get_type_intersection(srctype_base, desttype_base)
        if common_types:
            return common_types[0]

    # no common types
    return None

def _get_type_intersection(srctype: DataType, desttype: DataType) -> list[DataType]:
    # get types as list for each 
    srctype_list: list[DataType] = srctype.subtypes if isinstance(srctype, UnionType) else [srctype] # type: ignore
    desttype_list: list[DataType] = desttype.subtypes if isinstance(desttype, UnionType) else [desttype] # type: ignore

    # turn lists of instantiated DataTypes into type DataType
    srctypes: set[Type[DataType]] = set([type(x) for x in srctype_list])
    desttypes: set[Type[DataType]] = set([type(x) for x in desttype_list])
    
    # get intersection
    intersection = srctypes & desttypes

    # return original instantiated types
    out: list[DataType] = []
    for dtype in srctype_list:
        if type(dtype) in intersection:
            out.append(dtype)
    return out







