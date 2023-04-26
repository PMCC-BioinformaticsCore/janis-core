

# from typing import Optional, Any, Tuple
# from enum import Enum, auto

# from janis_core import ToolInput, CommandTool
# from janis_core.types import Boolean, Filename, DataType, File, Directory
# from janis_core import translation_utils as utils

# from ....variables import VariableManager
# from ....variables import VariableType
# from ....variables import VariableHistory
# from ....variables import Variable

# from .... import naming
# from .... import nfgen_utils
# from ....unwrap import unwrap_expression

# from .... import nulls
# from .ctype import CType, get_ctype


# from . import autofill


# """
# PRESCRIPT

# edge PS_POS_BASIC
# edge PS_OPT_BASIC
# edge PS_FLAG_TRUE
# edge PS_FLAG_FALSE

# general
# - (name) new varname
# - (src) current varname
# - (default) default value for tinput (including "" ie empty string)
# - (arr_join) how to format an array join for this tinput
# - (prefix) tinput prefix 
# - spacer (arg=value, arg:value, arg value) for tinputs with prefix 
# - delim (1,2,3 vs 1 2 3) delim for arrays
# - for each optional tinput, {src} != params.NULL ? {src or arr_join} : {default} (remember default can be "")
# - prefix each 
# - formats for all the above cases 

# ### PRE-SCRIPT ###

# declaration | def <name> = <value>

# value   | <arr_join>
#         | <condition>
         
# condition   | <check> ? <cond_true> : <cond_false>

# check   | <src> != params.NULL
#         | <src>.simpleName != params.NULL
#         | <src>[0].simpleName != params.NULL

# cond_true   | <src>
#             | <arr_join>
#             | "<prefix><spacer>" + <src>
#             | "<prefix><spacer>" + <arr_join>

# cond_false  | <default>
#             | ""

# ### SCRIPT ###

# value   | <prefix><spacer>${<src>}
#         | ${<src>}

# """

# FLAG_FALSE_FMT = 'def {name} = {src} ? "{prefix}" : ""'
# FLAG_TRUE_FMT  = 'def {name} = {src} == false ? "" : "{prefix}"'

# ARR_JOIN_FMT1 = "{src}.join('{delim}')"
# ARR_JOIN_FMT2 = "{pair1} + '{delim}' + {pair2}"
# ARR_JOIN_FMT3 = "{src}.collect{{ \"{prefix}\" + it }}" + ".join('{delim}')"

# DECLARATION_FMT = 'def {name} = {value}'
# DECLARATION_NAME_FMT1 = '{name}'
# DECLARATION_NAME_FMT2 = '{name}_joined'
# DECLARATION_VALUE_FMT = '{value}'

# CONDITION_FMT = '{cond_check} ? {cond_true} : {cond_false}'

# COND_CHECK_FMT1 = '{src} != params.NULL'
# COND_CHECK_FMT2 = '{src}.simpleName != params.NULL'
# COND_CHECK_FMT3 = '{src}[0].simpleName != params.NULL'
# COND_TRUE_FMT1 = '{src}'
# COND_TRUE_FMT2 = '"{prefix}{spacer}" + {src}'
# COND_FALSE_FMT = '{default}'







# class ScriptFormatter:
#     def __init__(self, 
#         tool: CommandTool,
#         variable_manager: VariableManager,
#     ) -> None:
#         self.tool = tool
#         self.variable_manager = variable_manager
        
#         self.dtt: DtypeType     # undefined at init
#         self.itype: CType       # undefined at init
#         self.tinput: ToolInput  # undefined at init

#     # PUBLIC METHOD
#     def format(self, tinput: ToolInput) -> Tuple[list[str], list[str]]:
#         self.tinput = tinput
#         self.itype = get_ctype(tinput)
#         self.dtt = get_dtype_type(self.dtype)

#         if self.should_ignore:
#             return self.format_ignored()
#         elif self.should_autofill:
#             return self.format_autofill()
#         elif self.is_special_case:
#             return self.format_special_case()
#         else:
#             return self.format_general()
    
#     @property
#     def should_ignore(self) -> bool:
#         # not supplied value in step call
#         if self.varhistory.original.vtype == VariableType.IGNORED:
#             if self.tinput.default is None:
#                 return True
#         return False
    
#     def format_ignored(self) -> Tuple[list[str], list[str]]:
#         prescript: list[str] = []
#         script: list[str] = []
#         return (prescript, script)
    
#     @property
#     def should_autofill(self) -> bool:
#         """
#         For ToolInputs which are not fed via a process input or a param,
#         should a static value be evaluated? 
#         """
#         if not self.should_ignore:
#             if self.varhistory.original.vtype in [VariableType.STATIC, VariableType.IGNORED]:
#                 return True
#         return False
    
#     @property 
#     def is_special_case(self) -> bool:
#         if self.itype in [
#             CType.FLAG_FALSE,
#             CType.FLAG_TRUE,
#             CType.POS_BASIC,
#             CType.OPT_BASIC
#         ]:
#             return True
#         return False
        
#     def format_special_case(self) -> Tuple[list[str], list[str]]:
#         prescript: list[str] = []
#         script: list[str] = []

#         if self.itype == CType.FLAG_FALSE:
#             src = self.varhistory.original.value
#             new_name = naming.process.generic(src)
#             ps_line = FLAG_FALSE_FMT.format(name=new_name, src=src, prefix=self.prefix_str)
#             prescript.append(ps_line)
#             sc_line = SCRIPT_FMT1.format(src=new_name)
#             script.append(sc_line)
        
#         elif self.itype == CType.FLAG_TRUE:
#             src = self.varhistory.original.value
#             new_name = naming.process.generic(src)
#             ps_line = FLAG_TRUE_FMT.format(name=new_name, src=src, prefix=self.prefix_str)
#             prescript.append(ps_line)
#             sc_line = SCRIPT_FMT1.format(src=new_name)
#             script.append(sc_line)
        
#         elif self.itype == CType.POS_BASIC:
#             src = self.varhistory.original.value
#             sc_line = SCRIPT_FMT1.format(src=src)
#             script.append(sc_line)
        
#         elif self.itype == CType.OPT_BASIC:
#             src = self.varhistory.original.value
#             sc_line = SCRIPT_FMT2.format(prefix=self.prefix_str, spacer=self.spacer_str, src=src)
#             script.append(sc_line)

#         return (prescript, script)
    
#     def format_general(self) -> Tuple[list[str], list[str]]:
#         prescript: list[str] = []
#         script: list[str] = []

#         # preprocessing step for array secondaries
#         # if utils.is_secondary_array_type(self.dtype) or utils.is_file_pair_array_type(self.dtype):
#         if self.dtt == DtypeType.SECONDARY_ARRAY:
#             line = self.gen_secondary_array_gather()
#             prescript.append(line)
        
#         # preprocessing step for file pair secondaries
#         # if utils.is_secondary_array_type(self.dtype) or utils.is_file_pair_array_type(self.dtype):
#         if self.dtt == DtypeType.FILE_PAIR_ARRAY:
#             line = self.gen_file_pair_array_gather()
#             prescript.append(line)

#         formatting_func = self.func_map[self.itype]
#         prescript_ln, script_ln = formatting_func()
#         if prescript_ln is not None:
#             prescript.append(prescript_ln)
#         if script_ln is not None:
#             script.append(script_ln)

#         return (prescript, script)

    
    

   


    
   
        



#     ### FORMATTING METHODS BY IType
#     def flag_true(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = self.generic_variable_name
#         prescript = PS_FLAG_TRUE.format(
#             name=new_varname, 
#             src=self.current_var_value, 
#             prefix=self.prefix_str
#         )
#         self.update_variable(new_varname)
#         script = SC_FLAG_TRUE.format(var=self.generic_variable_name)
#         return prescript, script

#     def flag_false(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = self.generic_variable_name
#         prescript = PS_FLAG_FALSE.format(
#             name=new_varname, 
#             src=self.current_var_value, 
#             prefix=self.prefix_str
#         )
#         self.update_variable(new_varname)
#         script = SC_FLAG_FALSE.format(var=self.current_var_value)
#         return prescript, script

#     def pos_basic(self) -> Tuple[Optional[str], Optional[str]]:
#         """
#         A basic positional. Has no prefix, and is mandatory.
#         Will have either a process input or param input, or will be fed a value via InputSelector.
#         """
#         prescript = None
#         script = SC_POS_BASIC.format(var=self.current_var_value)
#         return prescript, script

#     def pos_default(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = naming.process.generic(self.tinput)
#         prescript = PS_POS_DEFAULT.format(
#             name=new_varname,
#             src=self.current_var_value,
#             default=self.unwrapped_tinput_default
#         )
#         self.update_variable(new_varname)
#         script = SC_POS_DEFAULT.format(var=self.current_var_value)
#         return prescript, script

#     def pos_optional(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = self.generic_variable_name
#         if self.is_optional_filetype:
#             prescript = PS_POS_OPTIONAL_FILETYPES.format(
#                 name=new_varname, 
#                 src=self.current_var_value,
#                 default=self.optional_default
#             )
#         else:
#             prescript = PS_POS_OPTIONAL.format(
#                 name=new_varname, 
#                 src=self.current_var_value,
#             )
#         self.update_variable(new_varname)
#         script = SC_POS_OPTIONAL.format(var=self.current_var_value)
#         return prescript, script

#     def pos_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = f'{self.generic_variable_name}_joined'
#         prescript = PS_POS_BASIC_ARR.format(
#             name=new_varname, 
#             arr_join=self.arr_join
#         )
#         self.update_variable(new_varname)
#         script = SC_POS_BASIC_ARR.format(var=self.current_var_value)
#         return prescript, script

#     def pos_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = f'{self.generic_variable_name}_joined'
#         prescript = PS_POS_DEFAULT_ARR.format(
#             name=new_varname, 
#             src=self.current_var_value, 
#             arr_join=self.arr_join, 
#             default=self.unwrapped_tinput_default
#         )
#         self.update_variable(new_varname)
#         script = SC_POS_DEFAULT_ARR.format(var=self.current_var_value)
#         return prescript, script

#     def pos_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = f'{self.generic_variable_name}_joined'
#         if self.is_optional_filetype:
#             prescript = PS_POS_OPTIONAL_ARR_FILETYPES.format(
#                 name=new_varname, 
#                 src=self.current_var_value, 
#                 default=self.optional_default,
#                 arr_join=self.arr_join
#             )
#         else:
#             prescript = PS_POS_OPTIONAL_ARR.format(
#                 name=new_varname, 
#                 src=self.current_var_value, 
#                 arr_join=self.arr_join
#             )
#         self.update_variable(new_varname)
#         script = SC_POS_OPTIONAL_ARR.format(var=self.current_var_value)
#         return prescript, script

#     def opt_basic(self) -> Tuple[Optional[str], Optional[str]]:
#         prescript = None
#         script = SC_OPT_BASIC.format(prefix=self.prefix_str, var=self.current_var_value)
#         return prescript, script

#     def opt_default(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = self.generic_variable_name
#         prescript = PS_OPT_DEFAULT.format(
#             name=new_varname,
#             src=self.current_var_value,
#             default=self.unwrapped_tinput_default
#         )
#         self.update_variable(new_varname)
#         script = SC_OPT_DEFAULT.format(prefix=self.prefix_str, var=self.current_var_value)
#         return prescript, script

#     def opt_optional(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = self.generic_variable_name
#         if self.is_optional_filetype:
#             prescript = PS_OPT_OPTIONAL_FILETYPES.format(
#                 name=new_varname,
#                 src=self.current_var_value,
#                 default=self.optional_default,
#                 prefix=self.prefix_str
#             )
#         else:
#             prescript = PS_OPT_OPTIONAL.format(
#                 name=new_varname,
#                 src=self.current_var_value,
#                 prefix=self.prefix_str
#             )
#         self.update_variable(new_varname)
#         script = SC_OPT_OPTIONAL.format(var=self.current_var_value)
#         return prescript, script

#     def opt_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = f'{self.generic_variable_name}_joined'
#         prescript = PS_OPT_BASIC_ARR.format(
#             name=new_varname, 
#             arr_join=self.arr_join
#         )
#         self.update_variable(new_varname)
#         script = SC_OPT_BASIC_ARR.format(prefix=self.prefix_str, var=self.current_var_value)
#         return prescript, script

#     def opt_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = f'{self.generic_variable_name}_joined'
#         prescript = PS_OPT_DEFAULT_ARR.format(
#             name=new_varname, 
#             src=self.current_var_value, 
#             arr_join=self.arr_join, 
#             default=self.unwrapped_tinput_default
#         )
#         self.update_variable(new_varname)
#         script = SC_OPT_DEFAULT_ARR.format(prefix=self.prefix_str, var=self.current_var_value)
#         return prescript, script
    
#     def opt_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = f'{self.generic_variable_name}_joined'
#         if self.is_optional_filetype:
#             prescript = PS_OPT_OPTIONAL_ARR_FILETYPES.format(
#                 name=new_varname, 
#                 src=self.current_var_value, 
#                 default=self.optional_default,
#                 prefix=self.prefix_str,
#                 arr_join=self.arr_join
#             )
#         else:
#             prescript = PS_OPT_OPTIONAL_ARR.format(
#                 name=new_varname, 
#                 src=self.current_var_value, 
#                 prefix=self.prefix_str,
#                 arr_join=self.arr_join
#             )
#         self.update_variable(new_varname)
#         script = SC_OPT_OPTIONAL_ARR.format(var=self.current_var_value)
#         return prescript, script
    
#     def opt_basic_arr_prefixeach(self) -> Tuple[Optional[str], Optional[str]]:
#         # no prefix on front
#         new_varname = f'{self.generic_variable_name}_items'
#         prescript = PS_OPT_BASIC_ARR_PREFIXEACH.format(
#             name=new_varname, 
#             arr_join=self.arr_join
#         )
#         self.update_variable(new_varname)
#         script = SC_OPT_BASIC_ARR_PREFIXEACH.format(prefix=self.prefix_str, var=self.current_var_value)
#         return prescript, script

#     def opt_default_arr_prefixeach(self) -> Tuple[Optional[str], Optional[str]]:
#         # no prefix on front
#         new_varname = f'{self.generic_variable_name}_items'
#         prescript = PS_OPT_DEFAULT_ARR_PREFIXEACH.format(
#             name=new_varname, 
#             src=self.current_var_value, 
#             arr_join=self.arr_join, 
#             default=self.unwrapped_tinput_default
#         )
#         self.update_variable(new_varname)
#         script = SC_OPT_DEFAULT_ARR_PREFIXEACH.format(prefix=self.prefix_str, var=self.current_var_value)
#         return prescript, script

#     def opt_optional_arr_prefixeach(self) -> Tuple[Optional[str], Optional[str]]:
#         new_varname = f'{self.generic_variable_name}_items'
#         if self.is_optional_filetype:
#             prescript = PS_OPT_OPTIONAL_ARR_PREFIXEACH_FILETYPES.format(
#                 name=new_varname, 
#                 src=self.current_var_value, 
#                 default=self.optional_default,
#                 prefix=self.prefix_str,
#                 arr_join=self.arr_join
#             )
#         else:
#             prescript = PS_OPT_OPTIONAL_ARR_PREFIXEACH.format(
#                 name=new_varname, 
#                 src=self.current_var_value, 
#                 prefix=self.prefix_str,
#                 arr_join=self.arr_join
#             )
#         self.update_variable(new_varname)
#         script = SC_OPT_OPTIONAL_ARR_PREFIXEACH.format(var=self.current_var_value)
#         return prescript, script