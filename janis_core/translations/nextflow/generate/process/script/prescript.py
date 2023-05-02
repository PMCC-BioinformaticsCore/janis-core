

from typing import Any
from abc import ABC, abstractmethod

from janis_core import ToolInput, CommandTool, DataType
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from ....variables import VariableManager
from ....variables import VariableHistory
from ....variables import VariableType
from .... import naming

from .ctype import CType, get_ctype
from .attributes import get_attributes

from . import common
from . import selection
from . import ordering
 

def gen_prescript_lines(tool: CommandTool, vmanager: VariableManager) -> list[str]:
    lines: list[str] = []

    tinputs = selection.prescript_inputs(tool, vmanager)
    tinputs = ordering.prescript_inputs(tinputs)
    for tinput in tinputs:
        lines += gen_prescript_lines_for_input(tinput, tool, vmanager)
    
    return lines
        

def gen_prescript_lines_for_input(
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
        return True
    
    # tinput has consistent static value for each process call
    # no prescript, script autofilled
    elif varhistory.original.vtype == VariableType.STATIC:
        return True
    
    return False



class PreScriptFormatter(ABC):

    DECLARATION_FMT = 'def {dest} = {value}'
    CONDITION_FMT = '{cond_check} ? {cond_true} : {cond_false}'
    
    def __init__(
        self, 
        tinput: ToolInput,
        tool: CommandTool,
        vmanager: VariableManager,
    ) -> None:
        self.tool = tool
        self.tinput = tinput
        self.vmanager = vmanager
        self.ctype = get_ctype(tinput)
        self.attributes = get_attributes(tinput)
        self.dtt = utils.get_dtt(tinput.input_type)
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
    
    def default_str(self, apply_prefix: bool=False) -> str:
        default = None
        
        if self.varhistory.original.vtype == VariableType.STATIC:
            default = self.varhistory.original.value
        elif self.attributes.default:
            default = self.tinput.default
        
        if default is not None:
            # eval_cmdline uses value (parameter), self.prefix, self.spacer, self.delim etc
            expr = common.eval_cmdline_tinput(
                default=default,
                tinput=self.tinput,
                tool=self.tool,
                vmanager=self.vmanager,
                apply_prefix=apply_prefix
            ) 
            return expr
        
        return '""'


    def update_variable(self, new_value: Any) -> None:
        self.vmanager.update(
            tinput_id=self.tinput.id(),
            vtype_str='local',
            value=new_value
        )



### ARRAYS ### 

class GenericArrayFormatter(PreScriptFormatter):

    COND_CHECK_FMT = '{src} != params.NULL'
    COND_TRUE_FMT1 = '{src}'
    COND_TRUE_FMT2 = '"{prefix}{spacer}" + {src}'
    # COND_TRUE_FMT2 = '"{prefix}{spacer}${{{src}}}"'
    COND_FALSE_FMT = '{default}'

    ARR_JOIN_BASIC = '{src}.join(\'{delim}\')'
    ARR_JOIN_PREFIX = '"{prefix}{spacer}${{{src}.join(\'{delim}\')}}"'
    ARR_JOIN_PREFIXEACH = '{src}.collect{{ \"{prefix}{spacer}${{it}}\" }}.join(\'{delim}\')'

    def format(self) -> None:
        line = self.join_declaration()
        self.prescript.append(line)

    def join_declaration(self) -> str:
        new_varname = f'{naming.process.generic(self.tinput)}_joined'
        self.update_variable(new_varname)
        if self.attributes.optional or self.attributes.default:
            return self.join_declaration_optional()
        return self.join_declaration_mandatory()
    
    def join_declaration_mandatory(self) -> str:
        # eg: def adapters = "--adapters ${adapters.join(' ')}"
        dest = self.varhistory.items[1].value  # references the new local varname (eg adapters)
        value = self.arr_join()
        return self.DECLARATION_FMT.format(dest=dest, value=value)
    
    def join_declaration_optional(self) -> str:
        dest = self.varhistory.items[1].value  # references the second local varname (eg bams_joined)
        value = self.condition()
        return self.DECLARATION_FMT.format(dest=dest, value=value)
    
    def arr_join(self) -> str:
        prefix = self.prefix_str
        spacer = self.spacer_str
        src = self.varhistory.previous.value
        delim = self.delim_str

        if self.ctype in [
            CType.OPT_BASIC_ARR_PREFIXEACH,
            CType.OPT_DEFAULT_ARR_PREFIXEACH,
            CType.OPT_OPTIONAL_ARR_PREFIXEACH
        ]:
            # eg: adapters.collect{ "--adapters ${it}" }.join(' ')
            return self.ARR_JOIN_PREFIXEACH.format(src=src, prefix=prefix, spacer=spacer, delim=delim)
        
        else:
            # eg: "adapters.join(' ')"
            return self.ARR_JOIN_BASIC.format(src=src, delim=delim)

    def condition(self) -> str:
        # '{cond_check} ? {cond_true} : {cond_false}'
        # eg: adapters[0] != params.NULL ? "--adapters ${adapters.join(' ')}" : ""
        return self.CONDITION_FMT.format(
            cond_check=self.cond_check(),
            cond_true=self.cond_true(),
            cond_false=self.cond_false(),
        )

    def cond_check(self) -> str:
        # eg: adapters[0] != params.NULL
        src = self.varhistory.previous.value
        return self.COND_CHECK_FMT.format(src=src) 
    
    def cond_true(self) -> str:
        # eg: "--adapters ${adapters.join(' ')}"
        src = self.arr_join()

        if self.ctype == CType.OPT_OPTIONAL_ARR:
            return self.COND_TRUE_FMT2.format(prefix=self.prefix_str, spacer=self.spacer_str, src=src)
        else:
            return self.COND_TRUE_FMT1.format(src=src)

    def cond_false(self) -> str:
        if self.ctype == CType.OPT_OPTIONAL_ARR:
            expr = self.default_str(apply_prefix=True)
        else:
            expr = self.default_str(apply_prefix=False)
        
        return self.COND_FALSE_FMT.format(default=expr)
    


class SecondaryArrayFormatter(GenericArrayFormatter):
    
    SECONDARY_ARRAY_GATHER = 'def {dest} = get_primary_files({src}, {collate_size})'
    COND_CHECK_FMT = '{src}[0].simpleName != params.NULL'
    
    def format(self) -> None:
        line1 = self.gather_declaration()
        line2 = self.join_declaration()
        self.prescript.append(line1)
        self.prescript.append(line2)

    def gather_declaration(self) -> str:
        new_varname = naming.process.generic(self.tinput)
        self.update_variable(new_varname)
        src = self.varhistory.original.value   # references the process input varname
        dest = self.varhistory.items[1].value  # references the first local varname
        basetype = utils.get_base_type(self.dtype)
        collate_size = len(utils.get_extensions(basetype))
        return self.SECONDARY_ARRAY_GATHER.format(dest=dest, src=src, collate_size=collate_size)
    
    def join_declaration(self) -> str:
        new_varname = f'{self.varhistory.items[1].value}_joined' # references the first local varname
        self.update_variable(new_varname)
        if self.attributes.optional or self.attributes.default:
            return self.join_declaration_optional()
        return self.join_declaration_mandatory()
    
    def join_declaration_mandatory(self) -> str:
        # eg: def bams_joined = "--bams ${bams.join(' ')}"
        dest = self.varhistory.items[2].value  # references the second local varname (eg bams_joined)
        value = self.arr_join()
        return self.DECLARATION_FMT.format(dest=dest, value=value)
    
    def join_declaration_optional(self) -> str:
        # eg: def bams_joined = bams[0].simpleName != params.NULL ? "--bams ${bams.join(' ')}" : ""
        dest = self.varhistory.items[2].value  # references the second local varname (eg bams_joined)
        value = self.condition()
        return self.DECLARATION_FMT.format(dest=dest, value=value)


class FilePairArrayFormatter(PreScriptFormatter):

    FILE_PAIR_ARRAY_GATHER    = "def {dest} = {src}.collate(2, 2)"
    COND_CHECK_FMT            = '{src}[0].simpleName != params.NULL'
    ARR_JOIN_MANDATORY        = "{src}.collect{{ it.join('{delim}') }}"
    ARR_JOIN_MANDATORY_PREFIX = "{src}.collect{{ \"{prefix}{spacer}${{it.join('{delim}')}}\" }}"
    ARR_JOIN_OPTIONAL         = "{src}.collect{{ it[0].simpleName != params.NULL ? it.join('{delim}') : \"\" }}"
    ARR_JOIN_OPTIONAL_PREFIX  = "{src}.collect{{ it[0].simpleName != params.NULL ? \"{prefix}{spacer}${{it.join('{delim}')}}\" : \"\" }}"

    def format(self) -> None:
        line1 = self.gather_declaration()
        line2 = self.join_declaration()
        self.prescript.append(line1)
        self.prescript.append(line2)

    def gather_declaration(self) -> str:
        new_varname = naming.process.generic(self.tinput)
        self.update_variable(new_varname)
        src = self.varhistory.original.value   # references the process input varname
        dest = self.varhistory.items[1].value  # references the first local varname
        return self.FILE_PAIR_ARRAY_GATHER.format(dest=dest, src=src)
    
    def join_declaration(self) -> str:
        new_varname = f'{self.varhistory.items[1].value}_joined' # references the first local varname
        self.update_variable(new_varname)
        if self.attributes.prefix and (self.attributes.optional or self.attributes.default):
            return self.join_declaration_optional_prefix()
        elif self.attributes.optional or self.attributes.default:
            return self.join_declaration_optional_basic()
        elif self.attributes.prefix:
            return self.join_declaration_madatory_prefix()
        else:
            return self.join_declaration_mandatory_basic()
    
    def join_declaration_optional_prefix(self) -> str:
        # eg: def read_pairs_joined = read_pairs.collect{ it[0].simpleName != params.NULL ? "--reads ${it.join(' ')}" : "" }
        dest = self.varhistory.items[2].value  # references the second local varname (eg reads_joined)

        arr_join = self.ARR_JOIN_OPTIONAL_PREFIX.format(
            src=self.varhistory.items[1].value,
            prefix=self.prefix_str,
            spacer=self.spacer_str,
            delim=self.delim_str
        )

        return self.DECLARATION_FMT.format(dest=dest, value=arr_join)
    
    def join_declaration_optional_basic(self) -> str:
        dest = self.varhistory.items[2].value  # references the second local varname (eg reads_joined)

        arr_join = self.ARR_JOIN_OPTIONAL.format(
            src=self.varhistory.items[1].value,
            delim=self.delim_str
        )
        return self.DECLARATION_FMT.format(dest=dest, value=arr_join)
    
    def join_declaration_madatory_prefix(self) -> str:
        dest = self.varhistory.items[2].value  # references the second local varname (eg reads_joined)

        arr_join = self.ARR_JOIN_MANDATORY_PREFIX.format(
            src=self.varhistory.items[1].value,
            prefix=self.prefix_str,
            spacer=self.spacer_str,
            delim=self.delim_str
        )

        return self.DECLARATION_FMT.format(dest=dest, value=arr_join)
    
    def join_declaration_mandatory_basic(self) -> str:
        dest = self.varhistory.items[2].value  # references the second local varname (eg reads_joined)

        arr_join = self.ARR_JOIN_MANDATORY.format(
            src=self.varhistory.items[1].value,
            delim=self.delim_str
        )
        return self.DECLARATION_FMT.format(dest=dest, value=arr_join)
    


class FileArrayFormatter(GenericArrayFormatter):

    COND_CHECK_FMT = '{src}[0].simpleName != params.NULL'



class FlagArrayFormatter(PreScriptFormatter):

    def format(self) -> None:
        # prescript only needed when optional
        raise NotImplementedError




    
    
### SINGLES ###

class GenericFormatter(PreScriptFormatter):

    COND_CHECK_FMT = '{src} != params.NULL'
    COND_TRUE_FMT1 = '{src}'
    COND_TRUE_FMT2 = '"{prefix}{spacer}${{{src}}}"'
    COND_FALSE_FMT = '{default}'

    def format(self) -> None:
        # prescript only needed when optional
        if self.attributes.optional or self.attributes.default:
            line = self.optional_declaration()
            self.prescript.append(line)

    def optional_declaration(self) -> str:
        dest = naming.process.generic(self.tinput)
        self.update_variable(dest)
        condition = self.condition()
        return self.DECLARATION_FMT.format(dest=dest, value=condition)

    def condition(self) -> str:
        # '{cond_check} ? {cond_true} : {cond_false}'
        return self.CONDITION_FMT.format(
            cond_check=self.cond_check(),
            cond_true=self.cond_true(),
            cond_false=self.cond_false(),
        )

    def pvar(self) -> str:
        return self.varhistory.previous.value
    
    def cond_check(self) -> str:
        return self.COND_CHECK_FMT.format(src=self.pvar()) 
    
    def cond_true(self) -> str:
        src = self.pvar()

        if self.ctype == CType.OPT_OPTIONAL:
            return self.COND_TRUE_FMT2.format(prefix=self.prefix_str, spacer=self.spacer_str, src=src)
        
        elif self.ctype in [
            CType.POS_BASIC,
            CType.POS_DEFAULT,
            CType.POS_OPTIONAL,
            CType.OPT_BASIC,
            CType.OPT_DEFAULT,
        ]:
            return self.COND_TRUE_FMT1.format(src=src)
        
        else:
            raise RuntimeError(f'Component type {self.ctype} not allowed for this Formatter')

    def cond_false(self) -> str:
        if self.ctype == CType.OPT_OPTIONAL:
            expr = self.default_str(apply_prefix=True)
            return self.COND_FALSE_FMT.format(default=expr)

        elif self.ctype in [
            CType.POS_BASIC,
            CType.POS_DEFAULT,
            CType.POS_OPTIONAL,
            CType.OPT_BASIC,
            CType.OPT_DEFAULT,
        ]:
            expr = self.default_str(apply_prefix=False)
            return self.COND_FALSE_FMT.format(default=expr)
        
        else:
            raise RuntimeError(f'Component type {self.ctype} not allowed for this Formatter')
    

class FileFormatter(GenericFormatter):
    
    COND_CHECK_FMT = '{src}.simpleName != params.NULL'


class SecondaryFormatter(GenericFormatter):

    COND_CHECK_FMT = '{src}.simpleName != params.NULL'

    def pvar(self) -> str:
        return self.varhistory.previous.value[0]



class FilePairFormatter(PreScriptFormatter):

    ARR_JOIN_BASIC      = "{pair1} + '{delim}' + {pair2}"
    ARR_JOIN_PREFIX     = '"{prefix}{spacer}${{{pair1}}}{delim}${{{pair2}}}"'
    ARR_JOIN_PREFIXEACH  = '"{prefix}{spacer}${{{pair1}}}{delim}{prefix}{spacer}${{{pair2}}}"'
    
    COND_CHECK_FMT  = '{src}.simpleName != params.NULL'
    COND_TRUE_POS   = '{src}'
    COND_TRUE_OPT   = '"{prefix}{spacer}" + {src}'
    COND_FALSE_POS = '{default}'
    COND_FALSE_OPT = '{default}'
    
    def format(self) -> None:
        # join declaration
        self.join_declaration()
        
        # optional declaration for each file in pair if optional
        if self.attributes.optional:
            self.optional_declarations()

    def join_declaration(self) -> None:
        dest = f'{naming.process.generic(self.tinput)}_joined'
        self.update_variable(dest)
        
        if self.attributes.optional:
            line = self.join_declaration_optional()
        else:
            line = self.join_declaration_mandatory()
        
        self.prescript.append(line)

    def join_declaration_optional(self) -> str:
        dest = self.varhistory.current.value
        pair1 = self.varhistory.original.value[0]

        cond_check = self.COND_CHECK_FMT.format(src=pair1)
        cond_true = self.arr_join()
        cond_false = '""'
        condition = self.CONDITION_FMT.format(cond_check=cond_check, cond_true=cond_true, cond_false=cond_false)
        declaration = self.DECLARATION_FMT.format(dest=dest, value=condition)
        
        return declaration
    
    def join_declaration_mandatory(self) -> str:
        dest = self.varhistory.current.value
        value = self.arr_join() 
        return self.DECLARATION_FMT.format(dest=dest, value=value)
    
    def arr_join(self) -> str:
        pair1 = self.varhistory.original.value[0]
        pair2 = self.varhistory.original.value[1]

        if self.attributes.prefixeach:
            return self.ARR_JOIN_PREFIXEACH.format(
                prefix=self.prefix_str,
                spacer=self.spacer_str,
                pair1=pair1,
                delim=self.delim_str,
                pair2=pair2
            )

        elif self.attributes.prefix:
            return self.ARR_JOIN_PREFIX.format(
                prefix=self.prefix_str,
                spacer=self.spacer_str,
                pair1=pair1,
                delim=self.delim_str,
                pair2=pair2
            )
        
        else:
            return self.ARR_JOIN_BASIC.format(
                pair1=pair1,
                delim=self.delim_str,
                pair2=pair2
            )
    
    def optional_declarations(self) -> None:
        dest = [
            self.varhistory.original.value[0],
            self.varhistory.original.value[1],
        ]
        self.update_variable(dest)
        line1 = self.optional_declaration_pair(pair=0)
        line2 = self.optional_declaration_pair(pair=1)
        self.prescript.append(line1)
        self.prescript.append(line2)

    def optional_declaration_pair(self, pair: int) -> str:
        src = self.varhistory.original.value[pair]
        dest = self.varhistory.original.value[pair]

        ### CONDITION CHECK
        cond_check = self.COND_CHECK_FMT.format(src=src)
        
        ### CONDITION TRUE
        if self.attributes.prefix:
            cond_true = self.COND_TRUE_OPT.format(prefix=self.prefix_str, spacer=self.spacer_str, src=src)
        else:
            cond_true = self.COND_TRUE_POS.format(src=src)
        
        ### CONDITION FALSE
        if self.attributes.prefix:
            cond_false = self.COND_FALSE_OPT.format(prefix=self.prefix_str, spacer=self.spacer_str, default=self.default_str())
        else:
            cond_false = self.COND_FALSE_POS.format(default=self.default_str())

        condition = self.CONDITION_FMT.format(cond_check=cond_check, cond_true=cond_true, cond_false=cond_false)
        declaration = self.DECLARATION_FMT.format(dest=dest, value=condition)

        return declaration



class FlagFormatter(PreScriptFormatter):

    FLAG_FALSE_FMT = 'def {dest} = {src} ? "{prefix}" : ""'
    FLAG_TRUE_FMT  = 'def {dest} = {src} == false ? "" : "{prefix}"'

    def format(self) -> None:
        # prescript only needed when optional
        if self.attributes.optional or self.attributes.default:
            line = self.optional_declaration()
            self.prescript.append(line)

    def optional_declaration(self) -> str:
        dest = naming.process.generic(self.tinput)
        self.update_variable(dest)
        src = self.varhistory.previous.value
        prefix = self.prefix_str
        if str(self.tinput.default) == 'False':
            return self.FLAG_FALSE_FMT.format(dest=dest, src=src, prefix=prefix)
        return self.FLAG_TRUE_FMT.format(dest=dest, src=src, prefix=prefix)




