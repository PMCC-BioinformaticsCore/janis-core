


from typing import Optional, Any, Tuple

from janis_core import ToolInput, CommandTool
from janis_core.types import Boolean, Array, Filename, DataType

from ... import naming
from ... import nfgen_utils
from ...unwrap import unwrap_expression
from .itype import IType, get_itype


PS_FLAG_TRUE          = 'def {name} = {src} == false ? "" : "{prefix}"'
PS_FLAG_FALSE         = 'def {name} = {src} ? "{prefix}" : ""'
PS_POS_BASIC          = ''
PS_POS_DEFAULT        = 'def {name} = {src} ? {src} : {default}'
PS_POS_OPTIONAL       = 'def {name} = {src} ? {src} : ""'
PS_POS_BASIC_ARR      = 'def {name} = {arr_join}'
PS_POS_DEFAULT_ARR    = 'def {name} = {src} ? {arr_join} : "{default}"'
PS_POS_OPTIONAL_ARR   = 'def {name} = {src} ? {arr_join} : ""'
PS_OPT_BASIC          = ''
PS_OPT_DEFAULT        = 'def {name} = {src} ? {src} : {default}'
PS_OPT_OPTIONAL       = 'def {name} = {src} ? "{prefix}${{{src}}}" : ""'
PS_OPT_BASIC_ARR      = 'def {name} = {arr_join}'
PS_OPT_DEFAULT_ARR    = 'def {name} = {src} && {src}[0] != \'\' ? {arr_join} : "{default}"'
PS_OPT_OPTIONAL_ARR   = 'def {name} = {src} && {src}[0] != \'\' ? {arr_join} : ""'


SC_FLAG_TRUE          = '${{{var}}}'
SC_FLAG_FALSE         = '${{{var}}}'
SC_POS_BASIC          = '${{{var}}}'
SC_POS_DEFAULT        = '${{{var}}}'
SC_POS_OPTIONAL       = '${{{var}}}'
SC_POS_BASIC_ARR      = '${{{var}}}'
SC_POS_DEFAULT_ARR    = '${{{var}}}'
SC_POS_OPTIONAL_ARR   = '${{{var}}}'
SC_OPT_BASIC          = '{prefix}${{{var}}}'
SC_OPT_DEFAULT        = '{prefix}${{{var}}}'
SC_OPT_OPTIONAL       = '${{{var}}}'
SC_OPT_BASIC_ARR      = '{prefix}${{{var}}}'
SC_OPT_DEFAULT_ARR    = '{prefix}${{{var}}}'
SC_OPT_OPTIONAL_ARR   = '${{{var}}}'



class ScriptFormatter:
    def __init__(self, 
        tinput: ToolInput, 
        tool: CommandTool,
        process_inputs: set[str], 
        param_inputs: set[str], 
        internal_inputs: set[str], 
        sources: dict[str, Any]
    ) -> None:

        self.tinput = tinput
        self.tool = tool
        self.process_inputs = process_inputs    # ToolInput which has corresponding process input
        self.param_inputs = param_inputs        # ToolInput which has corresponding global param
        self.internal_inputs = internal_inputs  # ToolInput which has no corresponding process input or param
        self.sources = sources
        self.itype = get_itype(tinput)

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

            IType.OPT_BASIC_ARR:    self.opt_basic_arr,
            IType.OPT_DEFAULT_ARR:  self.opt_default_arr,
            IType.OPT_OPTIONAL_ARR: self.opt_optional_arr,
        }

    
    # PUBLIC METHOD
    def format(self) -> Tuple[Optional[str], Optional[str]]:
        prescript: Optional[str] = None
        script: Optional[str] = None
        
        if self.should_ignore:
            pass

        elif self.should_autofill:
            script = self.autofill_script_expr()
        
        else:
            formatting_func = self.func_map[self.itype]
            prescript, script = formatting_func()
        
        if nfgen_utils.is_array_secondary_type(self.dtype):
            assert(prescript)
            prescript = self.prepend_function_call(prescript, 'get_primary_files')
            print()

        return prescript, script

    def prepend_function_call(self, prescript: str, func_name: str) -> str:
        func_call = f'def {self.src} = {func_name}({self.src})'
        return f'{func_call}\n{prescript}'

    # HELPER PROPERTIES / METHODS
    @property
    def is_process_input(self) -> bool:
        if self.tinput.id() in self.process_inputs:
            return True
        return False
    
    @property
    def is_param_input(self) -> bool:
        if self.tinput.id() in self.param_inputs:
            return True
        return False
    
    @property
    def is_internal_input(self) -> bool:
        # ToolInput which references another ToolInput using InputSelector
        if self.tinput.id() in self.internal_inputs:
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
        
        Other cases for autofill exist including where values are driven using InputSelectors. 
        """
        if self.is_internal_input:
            if self.tinput.default is not None:
                return True
            elif isinstance(self.tinput.input_type, Filename):
                return True
        return False

    @property
    def src(self) -> Optional[str]:
        return naming.get_varname_toolinput(
            self.tinput,
            self.process_inputs,
            self.param_inputs,
            self.sources
        )
    
    @property
    def internal_name(self) -> Optional[str]:
        return naming.process_input_name(self.tinput)
    
    @property
    def dtype(self) -> DataType:
        return self.tinput.input_type
    
    @property
    def basetype(self) -> DataType:
        return nfgen_utils.get_base_type(self.dtype)

    @property
    def prefix(self) -> Optional[str]:
        if self.tinput.prefix is None:
            return None
        if isinstance(self.tinput.input_type, Boolean):
            return self.tinput.prefix
        elif self.tinput.separate_value_from_prefix == False:
            return self.tinput.prefix
        else:
            return f'{self.tinput.prefix}{self.delim}'
        
    @property
    def delim(self) -> str:
        return self.tinput.separator if self.tinput.separator else ' '
    
    @property
    def default(self) -> Any:
        default = self.tinput.default
        if default is None and isinstance(self.basetype, Filename):
            default = self.unwrap(self.basetype)
        elif self.itype in [IType.POS_DEFAULT_ARR, IType.OPT_DEFAULT_ARR]:
            default = self.eval_cmdline(default)
        else:
            default = self.unwrap(default)
        return default

    def unwrap(self, val: Any) -> Any:
        return unwrap_expression(
            val=val,
            tool=self.tool,
            sources=self.sources,
            process_inputs=self.process_inputs,
            param_inputs=self.param_inputs,
            internal_inputs=self.internal_inputs,
            in_shell_script=True,
        )
    
    @property
    def arr_join(self) -> Optional[str]:
        ARR_JOIN_BASIC      = "{src}.join('{delim}')"
        ARR_JOIN_PREFIX     = "\"{prefix}\" + {src}.join('{delim}')"
        ARR_JOIN_PREFIXEACH = "{src}.collect{{ \"{prefix}\" + it }}" + ".join('{delim}')"
        if isinstance(self.tinput.input_type, Array):
            if self.tinput.prefix and self.tinput.prefix_applies_to_all_elements:
                arr_join = ARR_JOIN_PREFIXEACH.format(src=self.src, prefix=self.prefix, delim=self.delim)
            elif self.tinput.prefix:
                arr_join = ARR_JOIN_PREFIX.format(prefix=self.prefix, src=self.src, delim=self.delim)
            else:
                arr_join = ARR_JOIN_BASIC.format(src=self.src, delim=self.delim)
        else:
            arr_join = None
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
            raise NotImplementedError


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
        prescript = PS_FLAG_TRUE.format(
            name=self.internal_name, 
            src=self.src, 
            prefix=self.prefix
        )
        script = SC_FLAG_TRUE.format(var=self.internal_name)
        return prescript, script

    def flag_false(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_FLAG_FALSE.format(
            name=self.internal_name, 
            src=self.src, 
            prefix=self.prefix
        )
        script = SC_FLAG_FALSE.format(var=self.internal_name)
        return prescript, script

    def pos_basic(self) -> Tuple[Optional[str], Optional[str]]:
        """
        A basic positional. Has no prefix, and is mandatory.
        Will have either a process input or param input, or will be fed a value via InputSelector.
        """
        prescript = None
        script = SC_POS_BASIC.format(var=self.src)
        return prescript, script

    def pos_default(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_POS_DEFAULT.format(
            name=self.internal_name,
            src=self.src,
            default=self.default
        )
        script = SC_POS_DEFAULT.format(var=self.internal_name)
        return prescript, script

    def pos_optional(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_POS_OPTIONAL.format(
            name=self.internal_name, 
            src=self.src
        )
        script = SC_POS_OPTIONAL.format(var=self.internal_name)
        return prescript, script

    def pos_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_POS_BASIC_ARR.format(
            name=self.internal_name, 
            arr_join=self.arr_join
        )
        script = SC_POS_BASIC_ARR.format(var=self.internal_name)
        return prescript, script

    def pos_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
        assert(self.src)
        assert(self.arr_join)
        prescript = PS_POS_DEFAULT_ARR.format(
            name=self.internal_name, 
            src=self.src, 
            arr_join=self.arr_join, 
            default=self.default
        )
        script = SC_POS_DEFAULT_ARR.format(var=self.internal_name)
        return prescript, script

    def pos_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_POS_OPTIONAL_ARR.format(
            name=self.internal_name, 
            src=self.src, 
            arr_join=self.arr_join
        )
        script = SC_POS_OPTIONAL_ARR.format(var=self.internal_name)
        return prescript, script

    def opt_basic(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = None
        script = SC_OPT_BASIC.format(prefix=self.prefix, var=self.src)
        return prescript, script

    def opt_default(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_OPT_DEFAULT.format(
            name=self.internal_name,
            src=self.src,
            default=self.default
        )
        script = SC_OPT_DEFAULT.format(prefix=self.prefix, var=self.internal_name)
        return prescript, script

    def opt_optional(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_OPT_OPTIONAL.format(
            name=self.internal_name,
            src=self.src,
            prefix=self.prefix
        )
        script = SC_OPT_OPTIONAL.format(var=self.internal_name)
        return prescript, script

    def opt_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_OPT_BASIC_ARR.format(
            name=self.internal_name, 
            arr_join=self.arr_join
        )
        script = SC_OPT_BASIC_ARR.format(prefix=self.prefix, var=self.internal_name)
        return prescript, script

    def opt_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_OPT_DEFAULT_ARR.format(
            name=self.internal_name, 
            src=self.src, 
            arr_join=self.arr_join, 
            default=self.default
        )
        script = SC_OPT_DEFAULT_ARR.format(prefix=self.prefix, var=self.internal_name)
        return prescript, script

    def opt_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = PS_OPT_OPTIONAL_ARR.format(
            name=self.internal_name, 
            src=self.src, 
            prefix=self.prefix,
            arr_join=self.arr_join
        )
        script = SC_OPT_OPTIONAL_ARR.format(var=self.internal_name)
        return prescript, script