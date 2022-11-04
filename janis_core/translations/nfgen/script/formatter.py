

from typing import Optional, Any, Tuple
from janis_core import ToolInput
from janis_core.types import Boolean, Array
from janis_core.translations.nfgen import utils
from janis_core.translations.nfgen import params

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
    IType.POS_DEFAULT_ARR:  'def {name} = {src} ? {arr_join} : {default}',
    IType.POS_OPTIONAL_ARR: 'def {name} = {src} ? {arr_join} : ""',

    IType.OPT_BASIC:    '',
    IType.OPT_DEFAULT:  'def {name} = {src} ? {src} : {default}',
    IType.OPT_OPTIONAL: 'def {name} = {src} ? "{prefix}${{{src}}}" : ""',

    IType.OPT_BASIC_ARR:    'def {name} = {arr_join}',
    IType.OPT_DEFAULT_ARR:  'def {name} = {src} ? {arr_join} : {default}',
    IType.OPT_OPTIONAL_ARR: 'def {name} = {src} ? "{prefix}" + {arr_join} : ""',
}

ARR_JOIN_BASIC      = '{src}.join({delim})'
ARR_JOIN_PREFIXEACH = '{src}.collect{{ {prefix} + it }}.join({delim})'


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


class InputFormatter:
    def __init__(self, 
        tinput: ToolInput, 
        process_inputs: set[str], 
        param_inputs: set[str], 
        scope: list[str]
    ) -> None:

        self.tinput = tinput
        self.process_inputs = process_inputs
        self.param_inputs = param_inputs
        self.scope = scope
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
        func = self.func_map[self.itype]
        prescript, script = func()
        return prescript, script

    
    ### HELPER PROPERTIES
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
    def name(self) -> str:
        return self.tinput.id()
        
    @property
    def delim(self) -> str:
        return self.tinput.separator if self.tinput.separator else ' '
    
    @property
    def default(self) -> Any:
        if self.tinput.default is not None:
            return utils.wrap_value(self.tinput.default, self.tinput)
        return None
    
    @property
    def src(self) -> str:
        if self.name in self.process_inputs:
            return self.name
        elif self.name in self.param_inputs:
            param = params.get(self.name, self.scope)
            return f'params.{param.name}'
        else:
            raise NotImplementedError
    
    @property
    def arr_join(self) -> Optional[str]:
        if isinstance(self.tinput, Array):
            if self.tinput.prefix_applies_to_all_elements:
                return ARR_JOIN_PREFIXEACH.format(self.src, self.prefix, self.delim)
            else:
                return ARR_JOIN_BASIC.format(self.src, self.delim)
        return None


    ### PRESCRIPT TEXT
    def flag_true(self) -> Tuple[Optional[str], Optional[str]]:
        # if fed by param or process input, there will be prescript var definition
        if not self.is_param_input and not self.is_process_input:
            prescript = None
            script = self.prefix
        else:
            prescript = self.prescript_template.format(
                name=self.name, 
                src=self.src, 
                prefix=self.prefix
            )
            script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def flag_false(self) -> Tuple[Optional[str], Optional[str]]:
        # if fed by param or process input, there will be prescript var definition
        if not self.is_param_input and not self.is_process_input:
            prescript = None
            script = self.prefix
        else:
            prescript = self.prescript_template.format(
                name=self.name, 
                src=self.src, 
                prefix=self.prefix
            )
            script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_basic(self) -> Tuple[Optional[str], Optional[str]]:
        # Mandatory: always fed value via process input or param input
        # TODO dubious about this
        prescript = None
        script = '${{{var}}}'.format(var=self.src)
        return prescript, script

    def pos_default(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.is_param_input and not self.is_process_input:
            # if internal input, value is directly injected
            prescript = None  
            script = self.default
        else:
            # if has param or is process input, there will be prescript var definition
            prescript = self.prescript_template.format(
                name=self.name,
                src=self.src,
                default=self.default
            )
            script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_optional(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.is_param_input and not self.is_process_input:
            # if internal input
            prescript = None  # no prescript needed - input is ignored
            script = None  # no script needed - input is ignored
        else:
            # there will be prescript var definition
            prescript = self.prescript_template.format(
                name=self.name, 
                src=self.src
            )
            script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
        # always has a prescript var to format array
        prescript = self.prescript_template.format(
            name=self.name, 
            arr_join=self.arr_join
        )
        script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.is_param_input and not self.is_process_input:
            # internal input
            prescript = None  # no prescript needed - value is directly injected
            script = self.default
        else:
            # fed from param or process input
            prescript = self.prescript_template.format(
                name=self.name, 
                src=self.src, 
                arr_join=self.arr_join, 
                default=self.default
            )
            script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def pos_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.is_param_input and not self.is_process_input:
            # internal input
            prescript = None  # no prescript needed - input is ignored
            script = None  # no script needed - input is ignored
        else:
            # fed by param or process input
            prescript = self.prescript_template.format(
                name=self.name, 
                src=self.src, 
                arr_join=self.arr_join
            )
            script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def opt_basic(self) -> Tuple[Optional[str], Optional[str]]:
        # Mandatory: always fed value via process input or param input
        # TODO DUBIOUSSSSSS
        prescript = None   # no prescript needed
        script = '{prefix}${{{var}}}'.format(prefix=self.prefix, var=self.src)
        return prescript, script

    def opt_default(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.is_param_input and not self.is_process_input:
            # internal
            prescript = None  # no prescript needed - value is directly injected
            script = '{prefix}{default}'.format(prefix=self.prefix, default=self.default)
        else:
            # fed by param or process input
            prescript = self.prescript_template.format(
                name=self.name,
                src=self.src,
                default=self.default
            )
            script = '{prefix}${{{var}}}'.format(prefix=self.prefix, var=self.name)
        return prescript, script

    def opt_optional(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.is_param_input and not self.is_process_input:
            # internal
            prescript = None  # no prescript needed - input is ignored
            script = None
        else:
            # external
            prescript = self.prescript_template.format(
                name=self.name,
                src=self.src,
                prefix=self.prefix
            )
            script = '${{{var}}}'.format(var=self.name)
        return prescript, script

    def opt_basic_arr(self) -> Tuple[Optional[str], Optional[str]]:
        # Mandatory: always fed value via process input or param input
        # always has a prescript var to format array
        # TODO DUBIOUSSSSSS
        prescript = self.prescript_template.format(
            name=self.name, 
            arr_join=self.arr_join
        )
        script = '{prefix}${{{var}}}'.format(prefix=self.prefix, var=self.name)
        return prescript, script

    def opt_default_arr(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.is_param_input and not self.is_process_input:
            # internal
            prescript = None  # no prescript needed - value is directly injected
            script = '{prefix}{default}'.format(prefix=self.prefix, default=self.default)
        else:
            # external
            prescript = self.prescript_template.format(
                name=self.name, 
                src=self.src, 
                arr_join=self.arr_join, 
                default=self.default
            )
            script = '{prefix}${{{var}}}'.format(prefix=self.prefix, var=self.name)
        return prescript, script

    def opt_optional_arr(self) -> Tuple[Optional[str], Optional[str]]:
        if not self.is_param_input and not self.is_process_input:
            # internal
            prescript = None  # no prescript needed
            script = None     # input is ignored
        else:
            # external
            prescript = self.prescript_template.format(
                name=self.name, 
                src=self.src, 
                prefix=self.prefix,
                arr_join=self.arr_join
            )
            script = '${{{var}}}'.format(var=self.name)
        return prescript, script



# bool input type identity
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
    if not tinput.input_type.optional:
        if tinput.default == None:
            return True
    return False

def is_optional(tinput: ToolInput) -> bool:
    if tinput.input_type.optional == True:
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
    process_inputs: set[str], 
    param_inputs: set[str],
    scope: list[str]
    ) -> Tuple[Optional[str], Optional[str]]:
    return InputFormatter(tinput, process_inputs, param_inputs, scope).format()

