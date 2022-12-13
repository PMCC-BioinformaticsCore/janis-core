

from typing import Optional, Any, Tuple
from janis_core import ToolInput, TInput, CommandTool
from janis_core.types import Boolean, Array, File, Filename
from .. import nfgen_utils
from .. import secondaries
from .. import params
from .. unwrap import unwrap_expression

from enum import Enum, auto


class IType(Enum):
    FLAG_TRUE           = auto()
    FLAG_FALSE          = auto()

    POS_BASIC           = auto()
    POS_DEFAULT         = auto()
    POS_OPTIONAL        = auto()

    POS_BASIC_ARR       = auto()
    POS_DEFAULT_ARR     = auto()
    POS_OPTIONAL_ARR    = auto()

    OPT_BASIC           = auto()
    OPT_DEFAULT         = auto()
    OPT_OPTIONAL        = auto()

    OPT_BASIC_ARR       = auto()
    OPT_DEFAULT_ARR     = auto()
    OPT_OPTIONAL_ARR    = auto()
    

### prescript formats
prescript_template_map = {
    IType.FLAG_TRUE:    'def {name} = {src} == false ? "" : "{prefix}"',
    IType.FLAG_FALSE:   'def {name} = {src} ? "{prefix}" : ""',

    IType.POS_BASIC:    '',
    IType.POS_DEFAULT:  'def {name} = {src} ? {src} : {default}',
    IType.POS_OPTIONAL: 'def {name} = {src} ? {src} : ""',

    IType.POS_BASIC_ARR:    'def {name} = {arr_join}',
    IType.POS_DEFAULT_ARR:  'def {name} = {src} ? {arr_join} : "{default}"',
    IType.POS_OPTIONAL_ARR: 'def {name} = {src} ? {arr_join} : ""',

    IType.OPT_BASIC:    '',
    IType.OPT_DEFAULT:  'def {name} = {src} ? {src} : {default}',
    IType.OPT_OPTIONAL: 'def {name} = {src} ? "{prefix}${{{src}}}" : ""',

    IType.OPT_BASIC_ARR:    'def {name} = {arr_join}',
    IType.OPT_DEFAULT_ARR:  'def {name} = {src} ? {arr_join} : "{default}"',
    IType.OPT_OPTIONAL_ARR: 'def {name} = {src} ? "{prefix}" + {arr_join} : ""',
}

ARR_JOIN_BASIC      = "{src}.join('{delim}')"
ARR_JOIN_PREFIXEACH = "{src}.collect{{ \"{prefix}\" + it }}.join('{delim}')"


def get_itype(tinput: ToolInput) -> IType:
    # flags
    if is_flag_true(tinput):
        return IType.FLAG_TRUE
    elif is_flag_false(tinput):
        return IType.FLAG_FALSE

    # positional arrays
    elif is_array(tinput) and is_positional(tinput) and is_basic(tinput):
        return IType.POS_BASIC_ARR
    elif is_array(tinput) and is_positional(tinput) and has_default(tinput):
        return IType.POS_DEFAULT_ARR
    elif is_array(tinput) and is_positional(tinput) and is_optional(tinput):
        return IType.POS_OPTIONAL_ARR
    
    # option arrays
    elif is_array(tinput) and is_option(tinput) and is_basic(tinput):
        return IType.OPT_BASIC_ARR
    elif is_array(tinput) and is_option(tinput) and has_default(tinput):
        return IType.OPT_DEFAULT_ARR
    elif is_array(tinput) and is_option(tinput) and is_optional(tinput):
        return IType.OPT_OPTIONAL_ARR
    
    # positionals
    elif is_positional(tinput) and is_basic(tinput):
        return IType.POS_BASIC
    elif is_positional(tinput) and has_default(tinput):
        return IType.POS_DEFAULT
    elif is_positional(tinput) and is_optional(tinput):
        return IType.POS_OPTIONAL
    
    # options
    elif is_option(tinput) and is_basic(tinput):
        return IType.OPT_BASIC
    elif is_option(tinput) and has_default(tinput):
        return IType.OPT_DEFAULT
    elif is_option(tinput) and is_optional(tinput):
        return IType.OPT_OPTIONAL

    else:
        raise NotImplementedError



### helper methods

def get_nf_variable_name(
    inp: ToolInput | TInput, 
    process_inputs: set[str], 
    param_inputs: set[str],
    sources: dict[str, Any]
    ) -> str:
    # get nextflow variable name which feeds data for this tool input
    if inp.id() in process_inputs:
        return get_src_process_input(inp)
    elif inp.id() in param_inputs:
        return get_src_param_input(inp, sources)
    else:
        raise NotImplementedError

def get_src_process_input(inp: ToolInput | TInput) -> str:
    # data fed via process input
    dtype = inp.input_type if isinstance(inp, ToolInput) else inp.intype # type: ignore
    basetype = nfgen_utils.get_base_type(dtype)
    # secondary files (name mapped to ext of primary file)
    if isinstance(basetype, File) and basetype.has_secondary_files():
        exts = secondaries.get_names(basetype)
        name = exts[0]
    # everything else
    else:
        name = inp.id()
    return name

def get_src_param_input(inp: ToolInput | TInput, sources: dict[str, Any]) -> str: 
    # data fed via global param
    src = sources[inp.id()]
    sel = src.source_map[0].source
    param = params.get(sel.input_node.uuid)
    return f'params.{param.name}'




class InputFormatter:
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
        self.prescript_template = prescript_template_map[self.itype]

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


    ### PUBLIC METHODS
    def format(self) -> Tuple[Optional[str], Optional[str]]:
        if self.should_ignore:
            prescript = None
            script = None
        
        elif self.should_autofill:
            prescript = None
            script = self.autofill_script_expr()
        
        else:
            formatting_func = self.func_map[self.itype]
            prescript, script = formatting_func()
        
        return prescript, script
    
    ### HELPER PROPERTIES
    @property
    def name(self) -> str:
        return self.tinput.id()

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
        dtype = self.tinput.input_type
        basetype = nfgen_utils.get_base_type(dtype)
        
        # get unwrapped default
        if default is None and isinstance(basetype, Filename):
            default = self.unwrap(basetype)
        else:
            default = self.unwrap(default)
        
        # if default is not None:
        #     default = nfgen_utils.to_groovy(default, self.tinput.input_type)
        
        return default
    
    @property
    def is_process_input(self) -> bool:
        if self.name in self.process_inputs:
            return True
        return False
    
    @property
    def is_param_input(self) -> bool:
        if self.name in self.param_inputs:
            return True
        return False
    
    @property
    def is_internal_input(self) -> bool:
        # ToolInput which references another ToolInput using InputSelector
        """
        ToolArgument(
            value=StringFormatter(
                "my string {VALUE}",
                VALUE=InputSelector(
                    input_to_select="ref_fasta"
                )
            )
        )
        ToolOutput(
            tag="outp",
            output_type=VcfTabix(),
            selector=WildcardSelector(
                wildcard=BasenameOperator(
                    InputSelector(input_to_select="input_file", type_hint=File())
                )
            ),
            doc=OutputDocumentation(doc=None),
        )

        ToolArgument: value
        ToolOutput: selector

        """
        if self.name in self.internal_inputs:
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
    def src(self) -> str:
        return get_nf_variable_name(
            inp=self.tinput,
            process_inputs=self.process_inputs,
            param_inputs=self.param_inputs,
            sources=self.sources
        )
    
    @property
    def arr_join(self) -> Optional[str]:
        if isinstance(self.tinput.input_type, Array):
            if self.tinput.prefix_applies_to_all_elements:
                return ARR_JOIN_PREFIXEACH.format(src=self.src, prefix=self.prefix, delim=self.delim)
            else:
                return ARR_JOIN_BASIC.format(src=self.src, delim=self.delim)
        return None

    
    ### HELPER METHODS
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

    def eval_cmdline(self, inp: ToolInput, val: Any) -> str:
        if isinstance(val, list):
            basetype = nfgen_utils.get_base_type(inp.input_type)
            
            if inp.prefix_applies_to_all_elements:
                vals_groovy = [nfgen_utils.to_groovy(elem, basetype) for elem in val]
                elems = [f'{self.prefix}{elem}' for elem in vals_groovy]
                cmdline = ' '.join(elems)
                return cmdline
            
            else:
                prefix = self.prefix if self.prefix else ''
                vals_groovy = [nfgen_utils.to_groovy(elem, basetype) for elem in val]
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
            expr = self.eval_cmdline(
                inp=self.tinput, 
                val=self.default
            )

        elif self.itype == IType.OPT_DEFAULT:
            expr = '{prefix}{default}'.format(
                prefix=self.prefix, 
                default=self.default
            )
        
        elif self.itype == IType.OPT_DEFAULT_ARR:
            # TODO TEST IMPORTANT
            expr = self.eval_cmdline(
                inp=self.tinput, 
                val=self.default
            )

        else:
            raise NotImplementedError(f'cannot autofill a {self.itype}')
        
        return expr
        



    ### FORMATTING METHODS BY IType
    def flag_true(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            src=self.src, 
            prefix=self.prefix
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def flag_false(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            src=self.src, 
            prefix=self.prefix
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_basic(self) -> Tuple[Optional[str], Optional[str]]:
        """
        A basic positional. Has no prefix, and is mandatory.
        Will have either a process input or param input, or will be fed a value via InputSelector.
        """
        prescript = None
        script = '${{{var}}}'.format(var=self.src)
        return prescript, script

    def pos_default(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name,
            src=self.src,
            default=self.default
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_optional(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            src=self.src
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            arr_join=self.arr_join
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            src=self.src, 
            arr_join=self.arr_join, 
            default=self.default
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            src=self.src, 
            arr_join=self.arr_join
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def opt_basic(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = None  
        script = '{prefix}${{{var}}}'.format(prefix=self.prefix, var=self.src)
        return prescript, script

    def opt_default(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name,
            src=self.src,
            default=self.default
        )
        script = '{prefix}${{{var}}}'.format(prefix=self.prefix, var=self.name)
        return prescript, script

    def opt_optional(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name,
            src=self.src,
            prefix=self.prefix
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def opt_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            arr_join=self.arr_join
        )
        script = '{prefix}${{{var}}}'.format(prefix=self.prefix, var=self.name)
        return prescript, script

    def opt_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            src=self.src, 
            arr_join=self.arr_join, 
            default=self.default
        )
        script = '{prefix}${{{var}}}'.format(prefix=self.prefix, var=self.name)
        return prescript, script

    def opt_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
        prescript = self.prescript_template.format(
            name=self.name, 
            src=self.src, 
            prefix=self.prefix,
            arr_join=self.arr_join
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script





# bool ToolInput type identities
def is_positional(tinput: ToolInput) -> bool:
    if tinput.prefix == None:
        return True
    return False

def is_flag(tinput: ToolInput) -> bool:
    if isinstance(tinput.input_type, Boolean):
        return True
    return False

def is_option(tinput: ToolInput) -> bool:
    if not is_positional(tinput) and not is_flag(tinput):
        return True
    return False


# bool features
def is_basic(tinput: ToolInput) -> bool:
    if not is_optional(tinput):
        if tinput.default == None:
            return True
    return False

def is_optional(tinput: ToolInput) -> bool:
    if isinstance(tinput.input_type, Filename):
        return False
    elif tinput.input_type.optional == True:
        return True
    return False

def has_default(tinput: ToolInput) -> bool:
    if tinput.default is not None:
        return True
    return False

def is_flag_true(tinput: ToolInput) -> bool:
    if is_flag(tinput) and str(tinput.default) != 'False':
        return True
    return False

def is_flag_false(tinput: ToolInput) -> bool:
    if is_flag(tinput) and not is_flag_true(tinput):
        return True
    return False

def is_array(tinput: ToolInput) -> bool:
    if tinput.input_type.is_array():
        return True
    return False


# module entry points

def format_input(
    tinput: ToolInput, 
    tool: CommandTool,
    process_inputs: set[str], 
    param_inputs: set[str],
    internal_inputs: set[str],
    sources: dict[str, Any]
    ) -> Tuple[Optional[str], Optional[str]]:
    return InputFormatter(
        tinput=tinput, 
        tool=tool, 
        process_inputs=process_inputs, 
        param_inputs=param_inputs, 
        internal_inputs=internal_inputs, 
        sources=sources
    ).format()

