


# from typing import Tuple, Any

# from janis_core import ToolInput, CommandTool

# from ....variables import VariableHistory
# from ....variables import VariableManager
# from ....variables import VariableType

# from .ctype import get_ctype
# from .ctype import CType
# from .dtypetype import get_dtype_type

# from . import common


# SCRIPT_FMT1 = '{src}' 
# SCRIPT_FMT2 = '{prefix}{spacer}{src}' 



# def autofill_script(tinput: ToolInput, tool: CommandTool, vmanager: VariableManager) -> list[str]:
#     # when autofill is possible, returns the str expression
#     # which can be injected directly into the nf script block.
#     # has no pre-script lines. 

#     formatter = AutoFiller(tool, tinput, vmanager)
#     return formatter.format()


# class AutoFiller:
#     def __init__(self, tool: CommandTool, tinput: ToolInput, vmanager: VariableManager) -> None:
#         self.tool = tool
#         self.tinput = tinput
#         self.vmanager = vmanager
#         self.itype = get_ctype(tinput)
#         self.dtt = get_dtype_type(tinput)

#     @property 
#     def varhistory(self) -> VariableHistory:
#         return self.vmanager.get(self.tinput.id())
    
#     @property
#     def value(self) -> Any:
#         # ignored input but has default
#         if self.varhistory.original.vtype == VariableType.IGNORED:
#             return self.tinput.default
#         # static value supplied to input
#         elif self.varhistory.original.vtype == VariableType.STATIC:
#             return self.varhistory.original.value
#         else:
#             raise RuntimeError

#     def format(self) -> list[str]:

#         if self.itype == CType.FLAG_TRUE:
#             return self.format_flag_true()
        
#         elif self.itype == CType.FLAG_FALSE:
#             return self.format_flag_false()

#         elif self.itype in [
#             CType.POS_BASIC,
#             CType.POS_BASIC_ARR,
#             CType.POS_DEFAULT,
#             CType.POS_OPTIONAL,
#             CType.POS_DEFAULT_ARR,
#             CType.POS_OPTIONAL_ARR,
#             CType.OPT_BASIC_ARR_PREFIXEACH,
#             CType.OPT_DEFAULT_ARR_PREFIXEACH,
#             CType.OPT_OPTIONAL_ARR_PREFIXEACH
#         ]:
#             return self.format_basic()
        
#         elif self.itype in [
#             CType.OPT_BASIC,
#             CType.OPT_BASIC_ARR,
#             CType.OPT_DEFAULT,
#             CType.OPT_OPTIONAL,
#             CType.OPT_DEFAULT_ARR,
#             CType.OPT_OPTIONAL_ARR
#         ]:
#             return self.format_prefix()

#         else:
#             raise NotImplementedError(f'cannot autofill a {self.itype}')

#     def format_flag_true(self) -> Tuple[list[str], list[str]]:
#         prescript: list[str] = []
#         script: list[str] = []
#         raise NotImplementedError
    
#     def format_flag_false(self) -> Tuple[list[str], list[str]]:
#         prescript: list[str] = []
#         script: list[str] = []
#         raise NotImplementedError
        
#     def format_basic(self) -> Tuple[list[str], list[str]]:
#         prescript: list[str] = []
#         script: list[str] = []
        
#         src = common.eval_default_cmdline(self.value)
#         sc_line = SCRIPT_FMT1.format(src=src)
#         script.append(sc_line)
        
#         return (prescript, script)
            
#     def format_prefix(self) -> Tuple[list[str], list[str]]:
#         prescript: list[str] = []
#         script: list[str] = []
        
#         src = common.eval_default_cmdline(self.value)
#         sc_line = SCRIPT_FMT2.format(
#             prefix=common.prefix_str(self.tinput), 
#             spacer=common.spacer_str(self.tinput), 
#             src=src
#         )
#         script.append(sc_line)
        
#         return (prescript, script)
        