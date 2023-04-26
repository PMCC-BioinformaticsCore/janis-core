


from typing import Optional
from dataclasses import dataclass
from abc import ABC, abstractmethod

"""
SECONDARY_ARRAY_GATHER = 'def {name} = get_primary_files({src})'
FILE_PAIR_ARRAY_GATHER = 'def {name} = {src}.collate(2, 1)'

FLAG_FALSE_FMT = 'def {name} = {src} ? "{prefix}" : ""'
FLAG_TRUE_FMT  = 'def {name} = {src} == false ? "" : "{prefix}"'

ARR_JOIN_BASIC                = "{src}.join('{delim}')"
ARR_JOIN_PREFIXEACH           = "{src}.collect{{ \"{prefix}\" + it }}" + ".join('{delim}')"
ARR_JOIN_FILE_PAIR            = "{pair1} + '{delim}' + {pair2}"
ARR_JOIN_FILE_PAIR_PREFIXEACH = '"{prefix}{spacer}${{{pair1}}} {prefix}{spacer}${{{pair2}}}"'

DECLARATION_FMT = 'def {name} = {value}'
DECLARATION_NAME_FMT1 = '{name}'
DECLARATION_NAME_FMT2 = '{name}_joined'
DECLARATION_VALUE_FMT = '{value}'

CONDITION_FMT = '{cond_check} ? {cond_true} : {cond_false}'

COND_CHECK_FMT1 = '{src} != params.NULL'
COND_CHECK_FMT2 = '{src}.simpleName != params.NULL'
COND_CHECK_FMT3 = '{src}[0].simpleName != params.NULL'
COND_TRUE_FMT1 = '{src}'
COND_TRUE_FMT2 = '"{prefix}{spacer}" + {src}'
COND_FALSE_FMT = '{default}'
"""


@dataclass
class Composer(ABC):
    
    @abstractmethod
    def compose(self) -> str:
        ...


@dataclass
class DeclarationComposer(Composer):
    flavour: str
    dest: str
    src: str
    prefix: Optional[str]=None

    SECONDARY_ARRAY_GATHER  = 'def {dest} = get_primary_files({src})'
    FILE_PAIR_ARRAY_GATHER  = 'def {dest} = {src}.collate(2, 1)'
    FLAG_FALSE_REDEFINE     = 'def {dest} = {src} ? "{prefix}" : ""'
    FLAG_TRUE_REDEFINE      = 'def {dest} = {src} == false ? "" : "{prefix}"'
    GENERIC_REDEFINE        = 'def {dest} = {src}'

    def compose(self) -> str:
        if self.flavour == 'secondary_array_gather':
            return self.SECONDARY_ARRAY_GATHER.format(dest=self.dest, src=self.src)
        
        elif self.flavour == 'file_pair_array_gather':
            return self.FILE_PAIR_ARRAY_GATHER.format(dest=self.dest, src=self.src)
        
        elif self.flavour == 'flag_false_redefine':
            return self.FLAG_FALSE_REDEFINE.format(dest=self.dest, src=self.src, prefix=self.prefix)
        
        elif self.flavour == 'flag_true_redefine':
            return self.FLAG_TRUE_REDEFINE.format(dest=self.dest, src=self.src, prefix=self.prefix)
        
        elif self.flavour == 'generic_redefine':
            return self.GENERIC_REDEFINE.format(dest=self.dest, src=self.src)
        
        else:
            raise ValueError(f"Unknown operation: {self.flavour}")
        

@dataclass
class ArrJoinComposer(Composer):
    flavour: str
    src: Optional[str]=None
    pair1: Optional[str]=None
    pair2: Optional[str]=None
    prefix: Optional[str]=None
    delim: Optional[str]=None
    spacer: Optional[str]=None

    BASIC                = "{src}.join('{delim}')"
    PREFIXEACH           = "{src}.collect{{ \"{prefix}\" + it }}" + ".join('{delim}')"
    FILE_PAIR            = "{pair1} + '{delim}' + {pair2}"
    FILE_PAIR_PREFIXEACH = '"{prefix}{spacer}${{{pair1}}} {prefix}{spacer}${{{pair2}}}"'
    
    def compose(self) -> str:
        if self.flavour == 'basic':
            return self.BASIC.format(src=self.src, delim=self.delim)
        
        elif self.flavour == 'prefixeach':
            return self.PREFIXEACH.format(src=self.src, prefix=self.prefix, delim=self.delim)
        
        elif self.flavour == 'file_pair':
            return self.FILE_PAIR.format(pair1=self.pair1, delim=self.delim, pair2=self.pair2)
        
        elif self.flavour == 'file_pair_prefixeach':
            return self.FILE_PAIR_PREFIXEACH.format(
                prefix=self.prefix, 
                spacer=self.spacer, 
                pair1=self.pair1, 
                pair2=self.pair2
            )
        
        else:
            raise ValueError(f"Unknown operation: {self.flavour}")


@dataclass
class ConditionComposer(Composer):
    vartype: str # generic, file, file_array
    src: str
    prefix: Optional[str]=None
    spacer: Optional[str]=None
    default: Optional[str]=None
    
    CONDITION_FMT = '{cond_check} ? {cond_true} : {cond_false}'
    
    CHECK_FMT1 = '{src} != params.NULL'
    CHECK_FMT2 = '{src}.simpleName != params.NULL'
    CHECK_FMT3 = '{src}[0].simpleName != params.NULL'
    
    COND_TRUE_FMT1 = '{src}'
    COND_TRUE_FMT2 = '"{prefix}{spacer}${{{src}}}"'
    
    COND_FALSE_FMT = '{default}'
    
    def compose(self) -> str:
        return self.CONDITION_FMT.format(
            self.cond_check,
            self.cond_true,
            self.cond_false
        )
    
    @property 
    def cond_check(self) -> str:
        if self.vartype == 'generic':
            return self.CHECK_FMT1.format(src=self.src)
        elif self.vartype == 'file':
            return self.CHECK_FMT2.format(src=self.src)
        elif self.vartype == 'file_array':
            return self.CHECK_FMT3.format(src=self.src)
        else:
            raise ValueError(f"Unknown vartype: {self.vartype}")
    
    @property 
    def cond_true(self) -> str:
        if self.prefix and self.spacer:
            return self.COND_TRUE_FMT2.format(prefix=self.prefix, spacer=self.spacer, src=self.src)
        else:
            return self.COND_TRUE_FMT1.format(src=self.src)
    
    @property 
    def cond_false(self) -> str:
        return self.default  # type: ignore


