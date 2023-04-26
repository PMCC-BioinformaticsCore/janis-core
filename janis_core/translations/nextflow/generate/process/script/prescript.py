

from typing import Optional
from abc import ABC, abstractmethod

from janis_core import ToolInput, CommandTool, DataType
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from ....variables import VariableManager
from ....variables import VariableHistory
from ....variables import VariableType
from .... import naming

# from .composer import DeclarationComposer
# from .composer import ArrJoinComposer
# from .composer import ConditionComposer

# from .ctype import CType, get_ctype
from .attributes import Attributes, get_attributes

from . import common


def gen_prescript_lines(
    tinput: ToolInput,
    tool: CommandTool,
    vmanager: VariableManager,
) -> list[str]:
    
    if _should_ignore(tinput, vmanager):
        return []
    
    dtt = utils.get_dtt(tinput.input_type)
    formatter_map = {
        DTypeType.SECONDARY_ARRAY: SecondaryArrayFormatter,
        DTypeType.SECONDARY: SecondaryFormatter,
        DTypeType.FILE_PAIR_ARRAY: FilePairArrayFormatter,
        DTypeType.FILE_PAIR: FilePairFormatter,
        DTypeType.FILE_ARRAY: FileArrayFormatter,
        DTypeType.FILE: FileFormatter,
        DTypeType.FLAG_ARRAY: FlagArrayFormatter,
        DTypeType.FLAG: FlagFormatter,
        DTypeType.GENERIC_ARRAY: GenericArrayFormatter,
        DTypeType.GENERIC: GenericFormatter,
    }

    formatter = formatter_map[dtt](tinput, tool, vmanager)
    formatter.format()
    return formatter.prescript


def _should_ignore(tinput: ToolInput, vmanager: VariableManager) -> bool:
    varhistory = vmanager.get(tinput.id())
    
    # tinput not supplied value in any process call, no default.
    # no prescript, no script
    if varhistory.original.vtype == VariableType.IGNORED:
        if tinput.default is None:
            return True
    
    # tinput has consistent static value for each process call
    # no prescript, script autofilled
    elif varhistory.original.vtype == VariableType.STATIC:
        return True
    
    return False



SECONDARY_ARRAY_GATHER = 'def {dest} = get_primary_files({src})'
FILE_PAIR_ARRAY_GATHER = 'def {dest} = {src}.collate(2, 2)'

FLAG_FALSE_FMT = 'def {dest} = {src} ? "{prefix}" : ""'
FLAG_TRUE_FMT  = 'def {dest} = {src} == false ? "" : "{prefix}"'

ARR_JOIN_BASIC                = "{src}.join('{delim}')"
ARR_JOIN_PREFIXEACH           = "{src}.collect{{ \"{prefix}\" + it }}" + ".join('{delim}')"
ARR_JOIN_FILE_PAIR            = "{pair1} + '{delim}' + {pair2}"
ARR_JOIN_FILE_PAIR_PREFIXEACH = '"{prefix}{spacer}${{{pair1}}} {prefix}{spacer}${{{pair2}}}"'

DECLARATION_FMT = 'def {dest} = {value}'
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


class PreScriptFormatter(ABC):
    
    def __init__(
        self, 
        tinput: ToolInput,
        tool: CommandTool,
        vmanager: VariableManager,
    ) -> None:
        self.tool = tool
        self.tinput = tinput
        self.vmanager = vmanager
        # self.ctype = get_ctype(tinput)
        self.attributes = get_attributes(tinput)
        self.prescript: list[str] = []

    @abstractmethod
    def format(self) -> None:
        ...

    @property
    def dtype(self) -> DataType:
        return self.tinput.input_type  # type: ignore
    
    @property
    def basetype(self) -> DataType:
        basetype = utils.get_base_type(self.dtype)
        basetype = utils.ensure_single_type(basetype)
        return basetype
    
    @property
    def varhistory(self) -> VariableHistory:
        return self.vmanager.get(self.tinput.id())
    
    @property
    def prefix_str(self) -> str:
        return common.prefix_str(self.tinput)
    
    @property
    def spacer_str(self) -> str:
        return common.spacer_str(self.tinput)
    
    @property
    def delim_str(self) -> str:
        return common.delim_str(self.tinput)
    
    @property
    def default_str(self) -> str:
        default = self.tinput.default
        if default is not None:
            # eval_cmdline uses value (parameter), self.prefix, self.spacer, self.delim etc
            default = common.eval_default_cmdline(
                default=default,
                tinput=self.tinput,
                tool=self.tool,
                vmanager=self.vmanager,
            ) 
        else:
            default = ""
        return default

    def update_variable(self, new_value: Optional[str]) -> None:
        self.vmanager.update(
            tinput_id=self.tinput.id(),
            vtype_str='local',
            value=new_value
        )



class SecondaryArrayFormatter(PreScriptFormatter):
    

    def format(self) -> None:
        line1 = self.gather_declaration()
        line2 = self.redefine_declaration()
        self.prescript.append(line1)
        self.prescript.append(line2)

    def gather_declaration(self) -> str:
        new_varname = naming.process.generic(self.tinput)
        self.update_variable(new_varname)
        dest = self.varhistory.current.value 
        src = self.varhistory.previous.value 
        return SECONDARY_ARRAY_GATHER.format(dest=dest, src=src)
    
    def redefine_declaration(self) -> str:
        new_varname = f'{self.varhistory.current.value}_joined'
        self.update_variable(new_varname)
        if self.attributes.optional:
            return self.redefine_declaration_optional()
        return self.redefine_declaration_mandatory()
    
    ### HELPERS ###
    def redefine_declaration_mandatory(self) -> str:
        dest = self.varhistory.current.value
        value = self.arr_join()
        return DECLARATION_FMT.format(dest=dest, value=value)
    
    def redefine_declaration_optional(self) -> str:
        dest = self.varhistory.current.value
        value = self.condition()
        return DECLARATION_FMT.format(dest=dest, value=value)

    def arr_join(self) -> str:

        if self.attributes.prefixeach:
            src = self.varhistory.previous.value
            prefix = self.prefix_str
            delim = self.delim_str
            return ARR_JOIN_PREFIXEACH.format(src=src, prefix=prefix, delim=delim)
        
        else:
            src = self.varhistory.previous.value
            delim = self.delim_str
            return ARR_JOIN_BASIC.format(src=src, delim=delim)
    
    def condition(self) -> str:
        return CONDITION_FMT.format(
            cond_check=self.cond_check(),
            cond_true=self.cond_true(),
            cond_false=self.cond_false(),
        )

    def cond_check(self) -> str:
        src = self.varhistory.previous.value
        return COND_CHECK_FMT3.format(src=src) 
    
    def cond_true(self) -> str:
        if self.attributes.prefix and not self.attributes.prefixeach:
            src = self.arr_join()
            prefix = self.prefix_str
            spacer = self.spacer_str
            return COND_TRUE_FMT2.format(src=src, prefix=prefix, spacer=spacer)
        else:
            src = self.arr_join()
            return COND_TRUE_FMT1.format(src=src)
    
    def cond_false(self) -> str:
        default = self.default_str
        return COND_FALSE_FMT.format(default=default)



class SecondaryFormatter(PreScriptFormatter):

    def format(self) -> None:
        # prescript only needed when optional
        pass


class FilePairArrayFormatter(PreScriptFormatter):

    def format(self) -> None:
        # prescript only needed when optional
        pass


class FilePairFormatter(PreScriptFormatter):

    def format(self) -> None:
        pass


class FileArrayFormatter(PreScriptFormatter):

    def format(self) -> None:
        pass


class FileFormatter(PreScriptFormatter):

    def format(self) -> None:
        pass


class FlagArrayFormatter(PreScriptFormatter):

    def format(self) -> None:
        # prescript mandatory: array join, possibly with optional check
        pass


class FlagFormatter(PreScriptFormatter):

    def format(self) -> None:
        # prescript only needed when optional
        pass


class GenericArrayFormatter(PreScriptFormatter):

    def format(self) -> None:
        # prescript mandatory: array join, possibly with optional check
        pass


class GenericFormatter(PreScriptFormatter):

    def format(self) -> None:
        # prescript only needed when optional
        pass


# ### MAIN FORMATTING CLASS ###

# class PreScriptFormatterOld:
#     def __init__(
#         self, 
#         tinput: ToolInput,
#         tool: CommandTool,
#         vmanager: VariableManager,
#     ) -> None:
#         self.tool = tool
#         self.tinput = tinput
#         self.vmanager = vmanager
#         self.ctype = get_ctype(tinput)
#         self.attributes = get_attributes(tinput)
#         self.prescript: list[str] = []

#         self.func_map = {
#             DTypeType.SECONDARY_ARRAY: self.secondary_array,
#             DTypeType.SECONDARY: self.secondary,
#             DTypeType.FILE_PAIR_ARRAY: self.file_pair_array,
#             DTypeType.FILE_PAIR: self.file_pair,
#             DTypeType.FILE_ARRAY: self.file_array,
#             DTypeType.FILE: self.file,
#             DTypeType.GENERIC_ARRAY: self.generic_array,
#             DTypeType.GENERIC: self.generic,
#         }

#     ### MAIN FORMATTING METHODS ###

#     def format(self) -> list[str]:
#         self.add_preprocessing_line()
#         formatting_func = self.func_map[self.dtt]
#         formatting_func()
#         return self.prescript
    
#     def secondary_array(self) -> None:
#         # declaration to gather the primary files
#         self.do_secondary_array_gather()
#         self.do_secondary_array_redefine()
    
#     def do_secondary_array_gather(self) -> None:
#         pass
    
#     def do_secondary_array_redefine(self) -> None:
#         pass

#     def get_secondary_array_join(self) -> str:
#         # arr join
#         composer = ArrJoinComposer(
#             flavour='prefixeach' if self.attributes.prefixeach else 'basic',
#             src=self.varhistory.current.value,  # type: ignore
#             prefix=self.prefix_str,
#             delim=self.delim_str
#         )
#         arr_join = composer.compose()

#         if self.optional:

#         # redefine
#         new_varname = f'{self.varhistory.current.value}_joined'
#         self.update_variable(new_varname)

#         composer = DeclarationComposer(
#             flavour='secondary_array_gather',
#             dest=self.varhistory.current.value,
#             src=self.varhistory.previous.value
#         )
#         line = composer.compose()
#         self.prescript.append(line)

#     def secondary(self) -> None:
#         pass
    
#     def file_pair_array(self) -> None:
#         pass
    
#     def file_pair(self) -> None:
#         # arr join
#         if self.attributes.prefixeach:
            
            
#         pass
    
#     def file_array(self) -> None:
#         pass
    
#     def file(self) -> None:
#         pass
    
#     def generic_array(self) -> None:
#         pass
    
#     def generic(self) -> None:
#         pass
    
#     @property 
#     def arr_join(self) -> Optional[str]:
#         pass

#     @property 
#     def arr_join_basic(self) -> Optional[str]:
#         pass
    
#     @property 
#     def arr_join_prefixeach(self) -> Optional[str]:
#         pass
    
#     @property 
#     def arr_join_file_pair(self) -> Optional[str]:
#         pass
    
#     @property 
#     def arr_join_file_pair_prefixeach(self) -> Optional[str]:
#         pass

#     def add_preprocessing_line(self) -> None:
#         # preprocessing step for array secondaries.
#         # need to add prescript line to gather primary files from flattened process input
#         if self.dtt == DTypeType.SECONDARY_ARRAY:
#             line = self.secondary_array_gather
#             self.prescript.append(line)
        
#         # preprocessing step for file pair secondaries.
#         # need to add prescript line to collect file pairs from flattened process input
#         if self.dtt == DTypeType.FILE_PAIR_ARRAY:
#             line = self.file_pair_array_gather
#             self.prescript.append(line)

#     def flag_false(self) -> None:
#         new_varname = naming.process.generic(self.tinput)
#         self.update_variable(new_varname)
#         declaration = FLAG_FALSE_FMT.format(
#             name=self.varhistory.previous.value, 
#             src=self.varhistory.current.value, 
#             prefix=self.prefix_str
#         )
#         self.prescript.append(declaration)
    
#     def flag_true(self) -> None:
#         new_varname = naming.process.generic(self.tinput)
#         self.update_variable(new_varname)
#         declaration = FLAG_TRUE_FMT.format(
#             name=self.varhistory.previous.value, 
#             src=self.varhistory.current.value, 
#             prefix=self.prefix_str
#         )
#         self.prescript.append(declaration)

    

#     ### HELPER PROPERTIES ###
    

    
#     @property
#     def prefixeach(self) -> bool:
#         if self.tinput.prefix_applies_to_all_elements:
#             return True
#         return False
    
#     def arr_join_str(self) -> str:
#         if self.attributes.prefixeach and self.dtt == DTypeType.FILE_PAIR:
#             return ARR_JOIN_FILE_PAIR_PREFIXEACH(
#                 prefix=self.prefix_str,
#                 spacer=self.spacer_str,
#                 pair1=
#             )

#         if self.attributes.prefixeach:
#             return ARR_JOIN_FMT3.format(
#                 src=self.varhistory.current.value, 
#                 prefix=self.prefix_str, 
#                 delim=self.delim_str
#             )
        
#         # file pair format
#         elif utils.is_file_pair_type(self.dtype):
#             input_var = self.varhistory.original
#             assert(input_var.value)
#             assert(len(input_var.value) == 2)
#             pair1 = input_var.value[0]
#             pair2 = input_var.value[1]
#             return ARR_JOIN_FMT2.format(
#                 pair1=pair1, 
#                 delim=self.delim_str, 
#                 pair2=pair2
#             )
        
#         # generic format
#         else:
#             return ARR_JOIN_FMT1.format(
#                 src=self.varhistory.current.value,
#                 delim=self.delim_str
#             )
    
#     def declaration_str(self) -> str:
#         new_varname = naming.process.generic(self.tinput)
#         self.update_variable(new_varname)
#         return DECLARATION_FMT.format(
#             name=self.declaration_name_str, 
#             value=self.declaration_value_str
#         )
    
#     def declaration_name_str(self) -> str:
#         if self.doing_join_operation:
#             return DECLARATION_NAME_FMT2.format(
#                 src=self.varhistory.current.value
#             )
#         else:
#             return DECLARATION_NAME_FMT1.format(
#                 src=self.varhistory.current.value
#             )
    
#     def declaration_value_str(self) -> str:
#         if self.doing_join_operation:
#             value = self.arr_join_str
#         else:
#             value = self.condition_str
#         return DECLARATION_VALUE_FMT.format(value=value)

#     def condition_str(self) -> str:
#         return CONDITION_FMT.format(
#             cond_check=self.cond_check_str, 
#             cond_true=self.cond_true_str, 
#             cond_false=self.cond_false_str
#         )
    
#     @property
#     def cond_check_str(self) -> str:
#         """
#         cond_check  | <src> != params.NULL
#                     | <src>.simpleName != params.NULL
#                     | <src>[0].simpleName != params.NULL
#         """
#         src = self.varhistory.current.value

#         if self.dtt in [DTypeType.SECONDARY_ARRAY, DTypeType.FILE_PAIR_ARRAY, DTypeType.FILE_ARRAY]:
#             # eg path secondary_array_flat
#             # eg path file_pair_array_flat
#             return COND_CHECK_FMT3.format(src=src)

#         elif self.dtt in [DTypeType.SECONDARY, DTypeType.FILE_PAIR, DTypeType.FILE]:
#             # eg file pair: Tuple path(pair1), path(pair2) (src=pair1)
#             # eg file:      path fastq (src=fastq)
#             return COND_CHECK_FMT2.format(src=src)

#         elif self.dtt in [DTypeType.GENERIC_ARRAY, DTypeType.GENERIC]:
#             # eg int:           val kmer_size (src=kmer_size)
#             # eg string:        val adapter (src=adapter)
#             # eg string_array:  val adapters (src=adapters) 
#             # (ignored optional array passed as params.NULL, not [params.NULL])
#             return COND_CHECK_FMT1.format(src=src)
        
#         else:
#             raise RuntimeError
           
#     @property
#     def cond_true_str(self) -> str:
#         """
#         cond_true   | <src>
#                     | <arr_join>
#                     | "<prefix><spacer>" + <src>
        
#         using "<prefix><spacer>" + <src>:
#         - PS_OPT_OPTIONAL
#         - PS_OPT_OPTIONAL_FILETYPES
#         - PS_OPT_OPTIONAL_ARR
#         - PS_OPT_OPTIONAL_ARR_FILETYPES
        
#         - prefix | optional | array? | !prefixeach

#         using <arr_join>:
#         - PS_POS_DEFAULT_ARR
#         - PS_POS_OPTIONAL_ARR
#         - PS_POS_OPTIONAL_ARR_FILETYPES
#         - PS_OPT_DEFAULT_ARR
#         - PS_OPT_DEFAULT_ARR_PREFIXEACH
#         - PS_OPT_OPTIONAL_ARR_PREFIXEACH
#         - PS_OPT_OPTIONAL_ARR_PREFIXEACH_FILETYPES
        
#         - prefix?  | optional | array | prefixeach?

#         using <src>:
#         - PS_POS_DEFAULT
#         - PS_POS_OPTIONAL
#         - PS_POS_OPTIONAL_FILETYPES
#         - PS_OPT_DEFAULT*
        
#         - prefix?  | optional | !array | !prefixeach
#         """

#         # using "<prefix><spacer>" + <src>:
#         # prefix | optional | array? | !prefixeach
#         if self.prefix_str and self.dtype.optional and not self.prefixeach:
#             return COND_TRUE_FMT2.format(
#                 prefix=self.prefix_str, 
#                 spacer=self.spacer_str, 
#                 src=self.varhistory.current.value
#             )
        
#         # using <arr_join>:
#         # prefix?  | optional | array | prefixeach?
#         elif self.dtype.optional and utils.is_array_type(self.dtype):
#             return COND_TRUE_FMT1.format(src=self.arr_join_str)

#         # using <src>:
#         # prefix?  | optional | !array | !prefixeach
#         elif self.dtype.optional and not utils.is_array_type(self.dtype) and not self.prefixeach:
#             return COND_TRUE_FMT1.format(src=self.varhistory.current.value)
        
#         else:
#             raise RuntimeError
        
#     @property
#     def cond_false_str(self) -> str:
#         return COND_FALSE_FMT.format(default=self.default_str)
    
#     @property
#     def secondary_array_gather(self) -> str:
#         new_varname = naming.process.generic(self.tinput)
#         self.update_variable(new_varname)
#         return SECONDARY_ARRAY_GATHER.format(
#             name=self.varhistory.previous.value,
#             src=self.varhistory.current.value
#         )
    
#     @property
#     def file_pair_array_gather(self) -> str:
#         new_varname = naming.process.generic(self.tinput)
#         self.update_variable(new_varname)
#         return FILE_PAIR_ARRAY_GATHER.format(
#             name=self.varhistory.previous.value,
#             src=self.varhistory.current.value
#         )
    
#     ### HELPER METHODS ###
    
#     def update_variable(self, new_value: Optional[str]) -> None:
#         self.vmanager.update(
#             tinput_id=self.tinput.id(),
#             vtype_str='local',
#             value=new_value
#         )