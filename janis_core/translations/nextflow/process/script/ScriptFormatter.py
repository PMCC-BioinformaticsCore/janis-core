

from typing import Optional, Any, Tuple

from janis_core import ToolInput, CommandTool
from janis_core.types import Boolean, Filename, DataType
from janis_core import translation_utils as utils

from ... import naming
from ... import nfgen_utils

from ...unwrap import unwrap_expression
from ...scope import Scope
from .. import data_sources

from .itype import IType, get_itype
from ..VariableManager import VariableManager


# pre-script (PS) formatting
PS_FLAG_TRUE                    = 'def {name} = {src} == false ? "" : "{prefix}"'
PS_FLAG_FALSE                   = 'def {name} = {src} ? "{prefix}" : ""'

PS_POS_BASIC                    = ''
PS_POS_DEFAULT                  = 'def {name} = {src} ? {src} : {default}'
PS_POS_OPTIONAL                 = 'def {name} = {src} ? {src} : ""'

PS_POS_BASIC_ARR                = 'def {name} = {arr_join}'
PS_POS_DEFAULT_ARR              = 'def {name} = {src} ? {arr_join} : "{default}"'
PS_POS_OPTIONAL_ARR             = 'def {name} = {src} ? {arr_join} : ""'

PS_OPT_BASIC                    = ''
PS_OPT_DEFAULT                  = 'def {name} = {src} ? {src} : {default}'
PS_OPT_OPTIONAL                 = 'def {name} = {src} ? "{prefix}${{{src}}}" : ""'

PS_OPT_BASIC_ARR                = 'def {name} = {arr_join}'
PS_OPT_DEFAULT_ARR              = 'def {name} = {src} ? {arr_join} : "{default}"'
PS_OPT_OPTIONAL_ARR             = 'def {name} = {src} ? "{prefix}" + {arr_join} : ""'

PS_OPT_BASIC_ARR_PREFIXEACH     = 'def {name} = {arr_join}'
PS_OPT_DEFAULT_ARR_PREFIXEACH   = 'def {name} = {src} ? {arr_join} : "{default}"'
PS_OPT_OPTIONAL_ARR_PREFIXEACH  = 'def {name} = {src} ? {arr_join} : ""'

# script (SC) formatting
SC_FLAG_TRUE                    = '${{{var}}}'
SC_FLAG_FALSE                   = '${{{var}}}'

SC_POS_BASIC                    = '${{{var}}}'
SC_POS_DEFAULT                  = '${{{var}}}'
SC_POS_OPTIONAL                 = '${{{var}}}'

SC_POS_BASIC_ARR                = '${{{var}}}'
SC_POS_DEFAULT_ARR              = '${{{var}}}'
SC_POS_OPTIONAL_ARR             = '${{{var}}}'

SC_OPT_BASIC                    = '{prefix}${{{var}}}'
SC_OPT_DEFAULT                  = '{prefix}${{{var}}}'
SC_OPT_OPTIONAL                 = '${{{var}}}'

SC_OPT_BASIC_ARR                = '{prefix}${{{var}}}'
SC_OPT_DEFAULT_ARR              = '{prefix}${{{var}}}'
SC_OPT_OPTIONAL_ARR             = '${{{var}}}'

SC_OPT_BASIC_ARR_PREFIXEACH     = '${{{var}}}'
SC_OPT_DEFAULT_ARR_PREFIXEACH   = '${{{var}}}'
SC_OPT_OPTIONAL_ARR_PREFIXEACH  = '${{{var}}}'



class ScriptFormatter:
    def __init__(self, 
        scope: Scope, 
        tool: CommandTool,
        variable_manager: VariableManager,
        sources: dict[str, Any]
    ) -> None:
        self.scope = scope
        self.tool = tool
        self.variable_manager = variable_manager
        self.sources = sources
        
        self.tinput: ToolInput  # undefined at init
        self.itype: IType       # undefined at init

        self.func_map = {
            IType.FLAG_TRUE:    self.flag_true,
            IType.FLAG_FALSE:   self.flag_false,
            
            IType.POS_BASIC:    self.pos_basic,
            IType.POS_DEFAULT:  self.pos_default,
            IType.POS_OPTIONAL: self.pos_optional,

            IType.POS_BASIC_ARR:    self.pos_basic_arr,
            IType.POS_DEFAULT_ARR:  self.pos_default_arr,
            IType.POS_OPTIONAL_ARR: self.pos_optional_arr,

            IType.OPT_BASIC:    self.opt_basic,
            IType.OPT_DEFAULT:  self.opt_default,
            IType.OPT_OPTIONAL: self.opt_optional,

            IType.OPT_BASIC_ARR:                self.opt_basic_arr,
            IType.OPT_DEFAULT_ARR:              self.opt_default_arr,
            IType.OPT_OPTIONAL_ARR:             self.opt_optional_arr,
            IType.OPT_BASIC_ARR_PREFIXEACH:     self.opt_basic_arr_prefixeach,
            IType.OPT_DEFAULT_ARR_PREFIXEACH:   self.opt_default_arr_prefixeach,
            IType.OPT_OPTIONAL_ARR_PREFIXEACH:  self.opt_optional_arr_prefixeach,
        }


    # PUBLIC METHOD
    def format(self, tinput: ToolInput) -> Tuple[list[str], list[str]]:
        self.tinput = tinput
        self.itype = get_itype(tinput)

        prescript: list[str] = []
        script: list[str] = []
        
        # preprocessing step for array secondaries
        if utils.is_array_secondary_type(self.dtype):
            func_call = self.gen_function_call('get_primary_files')
            prescript.append(func_call)
        
        # handling this ToolInput
        if self.should_ignore:
            return ([], [])

        elif self.should_autofill:
            expr = self.autofill_script_expr()
            if expr is not None:
                script.append(expr)
        
        else:
            formatting_func = self.func_map[self.itype]
            prescript_ln, script_ln = formatting_func()
            if prescript_ln is not None:
                prescript.append(prescript_ln)
            if script_ln is not None:
                script.append(script_ln)

        return prescript, script

    def gen_function_call(self, func_name: str) -> str:
        new_varname = self.generic_varname
        func_call = f'def {new_varname} = {func_name}({self.current_varname})'
        self.update_varname(new_varname)
        return f'{func_call}'

    # HELPER PROPERTIES / METHODS
    @property
    def is_process_input(self) -> bool:
        if self.tinput.id() in data_sources.process_inputs(self.scope):
            return True
        return False
    
    @property
    def is_param_input(self) -> bool:
        if self.tinput.id() in data_sources.param_inputs(self.scope):
            return True
        return False
    
    @property
    def is_internal_input(self) -> bool:
        # ToolInput which references another ToolInput using InputSelector
        if self.tinput.id() in data_sources.internal_inputs(self.scope):
            return True
        return False

    @property
    def should_ignore(self) -> bool:
        if self.is_internal_input and not self.should_autofill:
            return True
        return False
    
    @property
    def should_autofill(self) -> bool:
        """
        For ToolInputs which are not fed via a process input or a param,
        should a static value be evaluated? 
        eg. 
            ToolInput('myinp', String, default='hello!')
            Additional info: 
                In the step call for the Tool with ToolInput above, if no
                value is given for this ToolInput, we should use its default value. 
                This default value can be directly injected into the script. 
                If there was no default, we can't autofill anything.  
        
        Other cases for autofill exist including where values are driven using Filename. 
        """
        if self.is_internal_input:
            if self.tinput.default is not None:
                return True
            elif isinstance(self.tinput.input_type, Filename):
                return True
        return False

    @property
    def current_varname(self) -> Optional[str]:
        return self.variable_manager.current(self.tinput.id())
    
    @property
    def generic_varname(self) -> Optional[str]:
        return naming.process.generic(self.tinput)
    
    def update_varname(self, new_varname: Optional[str]) -> None:
        self.variable_manager.update(self.tinput.id(), new_varname)
    
    @property
    def dtype(self) -> DataType:
        return self.tinput.input_type  # type: ignore
    
    @property
    def basetype(self) -> DataType:
        basetype = utils.get_base_type(self.dtype)
        basetype = utils.ensure_single_type(basetype)
        return basetype

    @property
    def prefix(self) -> Optional[str]:
        if self.tinput.prefix is None:
            return None
        if isinstance(self.tinput.input_type, Boolean):
            return self.tinput.prefix
        elif self.tinput.separate_value_from_prefix == False:
            return self.tinput.prefix
        else:
            return f'{self.tinput.prefix} '
        
    @property
    def delim(self) -> str:
        return self.tinput.separator if self.tinput.separator else ' '
    
    @property
    def default(self) -> Any:
        default = self.tinput.default
        if default is None and isinstance(self.basetype, Filename):
            default = self.unwrap(self.basetype)
        elif self.itype in [IType.POS_DEFAULT_ARR, IType.OPT_DEFAULT_ARR, IType.OPT_DEFAULT_ARR_PREFIXEACH]:
            default = self.eval_cmdline(default)
        else:
            default = self.unwrap(default)
        return default

    def unwrap(self, val: Any) -> Any:
        return unwrap_expression(
            val=val,
            scope=self.scope,
            context='process_script',
            variable_manager=self.variable_manager,
            tool=self.tool,
            sources=self.sources,
            in_shell_script=True,
        )
    
    @property
    def arr_join(self) -> Optional[str]:
        ARR_JOIN_BASIC      = "{src}.join('{delim}')"
        # ARR_JOIN_PREFIX     = "\"{prefix}\" + {src}.join('{delim}')"
        ARR_JOIN_PREFIXEACH = "{src}.collect{{ \"{prefix}\" + it }}" + ".join('{delim}')"

        arr_join = None

        if self.itype.name.endswith('ARR') or self.itype.name.endswith('ARR_PREFIXEACH'):
            if self.itype in [IType.OPT_BASIC_ARR_PREFIXEACH, IType.OPT_DEFAULT_ARR_PREFIXEACH, IType.OPT_OPTIONAL_ARR_PREFIXEACH]:
                arr_join = ARR_JOIN_PREFIXEACH.format(src=self.current_varname, prefix=self.prefix, delim=self.delim)
            else:
                arr_join = ARR_JOIN_BASIC.format(src=self.current_varname, delim=self.delim)
        
        return arr_join
    
    def eval_cmdline(self, val: Any) -> str:
        if isinstance(val, list):
            if self.tinput.prefix_applies_to_all_elements:
                vals_groovy = [nfgen_utils.to_groovy(elem, self.basetype) for elem in val]
                elems = [f'{self.prefix}{elem}' for elem in vals_groovy]
                cmdline = ' '.join(elems)
                return cmdline
            
            else:
                prefix = self.prefix if self.prefix else ''
                vals_groovy = [nfgen_utils.to_groovy(elem, self.basetype) for elem in val]
                value = self.delim.join(vals_groovy)
                cmdline = f'{prefix}{value}'
                return cmdline
        else:
            prefix = self.prefix if self.prefix else ''
            value = nfgen_utils.to_groovy(val, self.basetype)
            cmdline = f'{prefix}{value}'
            return cmdline

    ### AUTOFILL METHODS BY IType
    def autofill_script_expr(self) -> Optional[str]: 
        # when autofill is possible, returns the str expression
        # which can be injected directly into the nf script block.
        if self.itype == IType.FLAG_TRUE:
            expr = self.prefix
        
        elif self.itype == IType.FLAG_FALSE:
            expr = None
        
        elif self.itype == IType.POS_BASIC:
            # TODO TEST IMPORTANT
            expr = self.default
        
        elif self.itype == IType.POS_BASIC_ARR:
            # TODO TEST IMPORTANT
            expr = self.default

        elif self.itype == IType.POS_DEFAULT:
            expr = self.default
        
        elif self.itype == IType.POS_DEFAULT_ARR:
            expr = self.default
        
        elif self.itype == IType.OPT_BASIC:
            expr = '{prefix}{default}'.format(
                prefix=self.prefix, 
                default=self.default
            )
        
        elif self.itype == IType.OPT_BASIC_ARR:
            # TODO TEST IMPORTANT
            expr = self.eval_cmdline(val=self.default)

        elif self.itype == IType.OPT_DEFAULT:
            expr = '{prefix}{default}'.format(
                prefix=self.prefix, 
                default=self.default
            )
        
        elif self.itype == IType.OPT_DEFAULT_ARR:
            # TODO TEST IMPORTANT
            expr = self.eval_cmdline(val=self.default)

        else:
            raise NotImplementedError(f'cannot autofill a {self.itype}')
        
        return expr


    ### FORMATTING METHODS BY IType
    def flag_true(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = self.generic_varname
        prescript = PS_FLAG_TRUE.format(
            name=new_varname, 
            src=self.current_varname, 
            prefix=self.prefix
        )
        self.update_varname(new_varname)
        script = SC_FLAG_TRUE.format(var=self.generic_varname)
        return prescript, script

    def flag_false(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = self.generic_varname
        prescript = PS_FLAG_FALSE.format(
            name=new_varname, 
            src=self.current_varname, 
            prefix=self.prefix
        )
        self.update_varname(new_varname)
        script = SC_FLAG_FALSE.format(var=self.current_varname)
        return prescript, script

    def pos_basic(self) -> Tuple[Optional[str], Optional[str]]:
        """
        A basic positional. Has no prefix, and is mandatory.
        Will have either a process input or param input, or will be fed a value via InputSelector.
        """
        prescript = None
        script = SC_POS_BASIC.format(var=self.current_varname)
        return prescript, script

    def pos_default(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = self.generic_varname
        prescript = PS_POS_DEFAULT.format(
            name=new_varname,
            src=self.current_varname,
            default=self.default
        )
        self.update_varname(new_varname)
        script = SC_POS_DEFAULT.format(var=self.current_varname)
        return prescript, script

    def pos_optional(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = self.generic_varname
        prescript = PS_POS_OPTIONAL.format(
            name=new_varname, 
            src=self.current_varname
        )
        self.update_varname(new_varname)
        script = SC_POS_OPTIONAL.format(var=self.current_varname)
        return prescript, script

    def pos_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = f'{self.generic_varname}_joined'
        prescript = PS_POS_BASIC_ARR.format(
            name=new_varname, 
            arr_join=self.arr_join
        )
        self.update_varname(new_varname)
        script = SC_POS_BASIC_ARR.format(var=self.current_varname)
        return prescript, script

    def pos_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = f'{self.generic_varname}_joined'
        prescript = PS_POS_DEFAULT_ARR.format(
            name=new_varname, 
            src=self.current_varname, 
            arr_join=self.arr_join, 
            default=self.default
        )
        self.update_varname(new_varname)
        script = SC_POS_DEFAULT_ARR.format(var=self.current_varname)
        return prescript, script

    def pos_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = f'{self.generic_varname}_joined'
        prescript = PS_POS_OPTIONAL_ARR.format(
            name=new_varname, 
            src=self.current_varname, 
            arr_join=self.arr_join
        )
        self.update_varname(new_varname)
        script = SC_POS_OPTIONAL_ARR.format(var=self.current_varname)
        return prescript, script

    def opt_basic(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = None
        script = SC_OPT_BASIC.format(prefix=self.prefix, var=self.current_varname)
        return prescript, script

    def opt_default(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = self.generic_varname
        prescript = PS_OPT_DEFAULT.format(
            name=new_varname,
            src=self.current_varname,
            default=self.default
        )
        self.update_varname(new_varname)
        script = SC_OPT_DEFAULT.format(prefix=self.prefix, var=self.current_varname)
        return prescript, script

    def opt_optional(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = self.generic_varname
        prescript = PS_OPT_OPTIONAL.format(
            name=new_varname,
            src=self.current_varname,
            prefix=self.prefix
        )
        self.update_varname(new_varname)
        script = SC_OPT_OPTIONAL.format(var=self.current_varname)
        return prescript, script

    def opt_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = f'{self.generic_varname}_joined'
        prescript = PS_OPT_BASIC_ARR.format(
            name=new_varname, 
            arr_join=self.arr_join
        )
        self.update_varname(new_varname)
        script = SC_OPT_BASIC_ARR.format(prefix=self.prefix, var=self.current_varname)
        return prescript, script

    def opt_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = f'{self.generic_varname}_joined'
        prescript = PS_OPT_DEFAULT_ARR.format(
            name=new_varname, 
            src=self.current_varname, 
            arr_join=self.arr_join, 
            default=self.default
        )
        self.update_varname(new_varname)
        script = SC_OPT_DEFAULT_ARR.format(prefix=self.prefix, var=self.current_varname)
        return prescript, script
    
    def opt_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = f'{self.generic_varname}_joined'
        prescript = PS_OPT_OPTIONAL_ARR.format(
            name=new_varname, 
            src=self.current_varname, 
            prefix=self.prefix,
            arr_join=self.arr_join
        )
        self.update_varname(new_varname)
        script = SC_OPT_OPTIONAL_ARR.format(var=self.current_varname)
        return prescript, script
    
    def opt_basic_arr_prefixeach(self) -> Tuple[Optional[str], Optional[str]]:
        # no prefix on front
        new_varname = f'{self.generic_varname}_items'
        prescript = PS_OPT_BASIC_ARR_PREFIXEACH.format(
            name=new_varname, 
            arr_join=self.arr_join
        )
        self.update_varname(new_varname)
        script = SC_OPT_BASIC_ARR_PREFIXEACH.format(prefix=self.prefix, var=self.current_varname)
        return prescript, script

    def opt_default_arr_prefixeach(self) -> Tuple[Optional[str], Optional[str]]:
        # no prefix on front
        new_varname = f'{self.generic_varname}_items'
        prescript = PS_OPT_DEFAULT_ARR_PREFIXEACH.format(
            name=new_varname, 
            src=self.current_varname, 
            arr_join=self.arr_join, 
            default=self.default
        )
        self.update_varname(new_varname)
        script = SC_OPT_DEFAULT_ARR_PREFIXEACH.format(prefix=self.prefix, var=self.current_varname)
        return prescript, script

    def opt_optional_arr_prefixeach(self) -> Tuple[Optional[str], Optional[str]]:
        new_varname = f'{self.generic_varname}_items'
        prescript = PS_OPT_OPTIONAL_ARR_PREFIXEACH.format(
            name=new_varname, 
            src=self.current_varname, 
            prefix=self.prefix,
            arr_join=self.arr_join
        )
        self.update_varname(new_varname)
        script = SC_OPT_OPTIONAL_ARR_PREFIXEACH.format(var=self.current_varname)
        return prescript, script