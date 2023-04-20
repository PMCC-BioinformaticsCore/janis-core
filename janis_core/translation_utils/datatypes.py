


from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from janis_core.types import DataType

from janis_core.types import (
    File, 
    Directory,
    UnionType,
    Filename
)


# GENERAL

def get_base_type(dtype: DataType) -> DataType:
    while dtype.name() == 'Array' and dtype.subtype():
        dtype = dtype.subtype()
    return dtype

def ensure_single_type(dtype: DataType) -> DataType:
    if isinstance(dtype, UnionType):
        return dtype.subtypes[0]
    return dtype    

def is_file_type(dtype: DataType) -> bool:
    basetype = get_base_type(dtype)
    if isinstance(basetype, (File, Filename, Directory)):
        return True
    elif is_file_pair_type(dtype):
        return True
    return False



# FILE PAIRS

known_file_pair_types = set([
    'FastqPair',
    'FastqGzPair',
])

def is_file_pair_type(dtype: DataType, recursive: bool=True) -> bool:
    if recursive:
        dtype = get_base_type(dtype)
    if dtype.name() in known_file_pair_types:
        return True
    return False

def is_array_file_pair_type(dtype: DataType) -> bool:
    if dtype.name() == 'Array':
        if is_file_pair_type(dtype):
            return True
    return False



### SECONDARIES 

def is_secondary_type(dtype: DataType) -> bool:
    if not is_array_secondary_type(dtype):
        basetype = get_base_type(dtype)
        if isinstance(basetype, File) and basetype.has_secondary_files():
            return True
    return False

def is_array_secondary_type(dtype: DataType) -> bool:
    basetype = get_base_type(dtype)
    if dtype.is_array() and isinstance(basetype, File) and basetype.has_secondary_files():
        return True
    return False

def get_extensions(dtype: File, remove_prefix_symbols: bool=False) -> list[str]:
    """returns extension of each file for File types with secondaries"""
    primary_ext: str = ''
    secondary_exts: list[str] = []

    # primary extension
    if isinstance(dtype, Directory):
        return []
    if len(dtype.get_extensions()) > 0:
        primary_ext = dtype.get_extensions()[0]
    else:
        primary_ext = 'primary'
    
    # secondary extensions
    if dtype.secondary_files() is not None:
        secondary_exts = dtype.secondary_files()
    else:
        secondary_exts = []

    exts = _sort_extensions(primary_ext, secondary_exts)
    if remove_prefix_symbols:
        exts = [x.lstrip('.^') for x in exts]
    return exts

def _sort_extensions(primary_ext: str, secondary_exts: list[str]) -> list[str]:
    out: list[str] = []
    out.append(primary_ext)
    secondary_exts = sorted(secondary_exts, key=lambda x: x.rsplit('.')[-1])
    out += secondary_exts
    return out
