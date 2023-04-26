


# from enum import Enum, auto

# from janis_core import ToolInput

# "prefix | optional | array? | !prefixeach"

# class MType(Enum):
#     PREFIX                           = auto()
#     PREFIX_OPTIONAL                  = auto()
#     PREFIX_OPTIONAL_ARRAY            = auto()

    
#     PREFIX_OPTIONAL_ARRAY_PREFIXEACH = auto()
#     PREFIX_OPTIONAL_ARRAY_xxxxxxxxxx = auto()
#     PREFIX_OPTIONAL_xxxxx_xxxxxxxxxx = auto()
#     PREFIX_xxxxxxxx_xxxxx_xxxxxxxxxx = auto()
    

#     OPTIONAL                = auto()
#     OPTIONAL_ARRAY          = auto()

#     ARRAY



#     ARRAY                   = auto()
#     PREFIXEACH              = auto()
    


# # public function
# def get_mtype(tinput: ToolInput) -> MType:
#     # flags
#     if is_flag_true(tinput):
#         return CType.FLAG_TRUE
#     elif is_flag_false(tinput):
#         return CType.FLAG_FALSE

#     # positional arrays
#     elif is_array_type(tinput) and is_positional(tinput) and is_basic(tinput):
#         return CType.POS_BASIC_ARR
#     elif is_array_type(tinput) and is_positional(tinput) and has_default(tinput):
#         return CType.POS_DEFAULT_ARR
#     elif is_array_type(tinput) and is_positional(tinput) and is_optional(tinput):
#         return CType.POS_OPTIONAL_ARR
    
#     # option arrays
#     elif is_array_type(tinput) and is_option(tinput) and is_basic(tinput) and tinput.prefix_applies_to_all_elements:
#         return CType.OPT_BASIC_ARR_PREFIXEACH
#     elif is_array_type(tinput) and is_option(tinput) and has_default(tinput) and tinput.prefix_applies_to_all_elements:
#         return CType.OPT_DEFAULT_ARR_PREFIXEACH
#     elif is_array_type(tinput) and is_option(tinput) and is_optional(tinput) and tinput.prefix_applies_to_all_elements:
#         return CType.OPT_OPTIONAL_ARR_PREFIXEACH
#     elif is_array_type(tinput) and is_option(tinput) and is_basic(tinput):
#         return CType.OPT_BASIC_ARR
#     elif is_array_type(tinput) and is_option(tinput) and has_default(tinput):
#         return CType.OPT_DEFAULT_ARR
#     elif is_array_type(tinput) and is_option(tinput) and is_optional(tinput):
#         return CType.OPT_OPTIONAL_ARR
    
#     # positionals
#     elif is_positional(tinput) and is_basic(tinput):
#         return CType.POS_BASIC
#     elif is_positional(tinput) and has_default(tinput):
#         return CType.POS_DEFAULT
#     elif is_positional(tinput) and is_optional(tinput):
#         return CType.POS_OPTIONAL
    
#     # options
#     elif is_option(tinput) and is_basic(tinput):
#         return CType.OPT_BASIC
#     elif is_option(tinput) and has_default(tinput):
#         return CType.OPT_DEFAULT
#     elif is_option(tinput) and is_optional(tinput):
#         return CType.OPT_OPTIONAL

#     else:
#         # some future IType
#         raise NotImplementedError