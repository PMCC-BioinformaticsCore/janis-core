
from typing import Any, Tuple, Optional

from janis_core.workflow.workflow import InputNode, StepNode
from janis_core import (
    File,
    CommandTool,
    ToolArgument,
    ToolInput,
    TInput,
)

from .. import ordering
from .. import settings
from .. import channels
from .. import params
from .. import nfgen_utils as nfgen_utils
from janis_core.translations.nfgen.unwrap import unwrap_expression

from .script_formatting import format_input



def gen_script_for_cmdtool(
    tool: CommandTool,
    input_in_selectors: dict[str, Any],
    stdout_filename: str,
    scope: list[str],
    values: dict[str, Any],
) -> Tuple[Optional[str], str]:
    return ProcessScriptGenerator(
        tool=tool,
        input_in_selectors=input_in_selectors,
        stdout_filename=stdout_filename,
        scope=scope,
        values=values,
    ).generate()



class ProcessScriptGenerator:
    def __init__(
        self,
        tool: CommandTool, 
        input_in_selectors: dict[str, Any],
        stdout_filename: str,
        scope: list[str],
        values: Optional[dict[str, Any]]=None,
    ):
        assert(tool)
        self.tool = tool
        self.scope = scope
        self.process_name = scope[-1] if scope else tool.id()
        self.input_in_selectors = input_in_selectors
        self.stdout_filename = stdout_filename

        # think this is ok?
        self.values = values if values is not None else {}
        self.process_inputs = get_process_inputs(self.values)
        self.param_inputs = get_param_inputs(self.values)
        self.internal_inputs = get_internal_input_ids(tool, self.values)

        self.prescript: list[str] = []
        self.script: list[str] = []

    def generate(self) -> Tuple[Optional[str], str]:
        """Generate the script content of a Nextflow process for Janis command line tool"""
        self.handle_cmdtool_preprocessing()
        self.handle_cmdtool_base_command()
        self.handle_cmdtool_inputs()
        prescript = self.finalise_prescript()
        script = self.finalise_script()
        return prescript, script

    def handle_cmdtool_preprocessing(self) -> None:
        for dirpath in self.tool.directories_to_create() or []:
            unwrapped_dir = unwrap_expression(
                value=dirpath, 
                input_in_selectors=self.input_in_selectors,
                inputs_dict=self.tool.inputs_map(),
                tool=self.tool, 
                in_shell_script=True
            ) 
            line = f"mkdir -p '{unwrapped_dir}'"
            self.script.append(line)

    def handle_cmdtool_base_command(self) -> None:
        bc = self.tool.base_command()
        if bc is not None:
            if isinstance(bc, list):
                self.script += [' '.join([str(cmd) for cmd in bc])]
            else:
                self.script += [str(bc)]

    def handle_cmdtool_inputs(self) -> None:
        for inp in ordering.cmdtool_inputs_arguments(self.tool):
            match inp:
                case ToolInput():
                    prescript_ln, script_ln = format_input(
                        inp, 
                        self.process_inputs, 
                        self.param_inputs, 
                        scope=self.scope,
                        sources=self.values
                    )
                    if prescript_ln is not None:
                        self.prescript.append(prescript_ln)
                    if script_ln is not None:
                        self.script.append(script_ln)
                    
                # arguments
                case _:
                    self.handle_tool_argument(inp)
        
    def handle_tool_argument(self, arg: ToolArgument) -> None:
        expression = unwrap_expression(
            value=arg.value,
            input_in_selectors=self.input_in_selectors,
            tool=self.tool,
            inputs_dict=self.tool.inputs_map(),
            skip_inputs_lookup=True,
            quote_string=False,
            in_shell_script=True,
        )
        if arg.prefix is not None:
            space = " " if arg.separate_value_from_prefix else ""
            line = f'{arg.prefix}{space}"{expression}"'
        else:
            line = expression
        self.script.append(line)

    def finalise_prescript(self) -> Optional[str]:
        if self.prescript:
            return '\n'.join(self.prescript)
        return None

    def finalise_script(self) -> str:
        script = self.script
        if settings.JANIS_ASSISTANT:
            script = script + [f'| tee {self.stdout_filename}_{self.process_name}']
        script = [f'{ln} \\' for ln in script]
        return '\n'.join(script)







### process inputs

def get_process_inputs(sources: dict[str, Any]) -> set[str]:
    """
    determine the tool inputs which should remnain as process inputs
    """
    if settings.MODE == 'workflow':
        return get_process_inputs_workflowmode(sources)
    elif settings.MODE == 'tool':  # type: ignore
        return get_process_inputs_toolmode(sources)
    else:
        raise RuntimeError

def get_process_inputs_workflowmode(sources: dict[str, Any]) -> set[str]:
    """
    inputs which are fed (via step inputs) using a file type workflow input
    inputs which are fed (via step inputs) using a connection
    inputs which are involved in scatter
    """
    channel_wfinp_ids = get_channel_process_inputs(sources)
    step_conn_ids = get_connection_process_inputs(sources)
    scatter_wfinp_ids = get_scatter_process_inputs(sources)
    surviving_ids = channel_wfinp_ids | step_conn_ids | scatter_wfinp_ids
    return surviving_ids

def get_process_inputs_toolmode(sources: dict[str, Any]) -> set[str]:
    """
    inputs which are file types
    non-file types usually fed values from params instead.
    
    if MINIMAL_PROCESS:
        - remove inputs which are optional
        - remove inputs with defaults
    """
    all_inputs: list[TInput] = list(tool.inputs_map().values())
    
    surviving_ids = get_all_input_ids(all_inputs)
    file_ids = get_file_input_ids(all_inputs)
    optional_ids = get_optional_input_ids(all_inputs)
    default_ids = get_default_input_ids(all_inputs)
    
    if settings.MINIMAL_PROCESS:
        surviving_ids = surviving_ids & file_ids
        surviving_ids = surviving_ids - optional_ids
        surviving_ids = surviving_ids - default_ids
    else:
        surviving_ids = surviving_ids & file_ids

    return surviving_ids




### param inputs

def get_param_inputs(sources: dict[str, Any]) -> set[str]:
    """
    determine the tool inputs which should be fed a value via params
    """
    if settings.MODE == 'workflow':
        return get_param_inputs_workflowmode(sources)
    elif settings.MODE == 'tool':  # type: ignore
        return get_param_inputs_toolmode(sources)
    else:
        raise RuntimeError

def get_param_inputs_workflowmode(sources: dict[str, Any]) -> set[str]:
    """
    inputs which are fed (via step inputs) using a non-File type workflow input
    """
    if settings.MINIMAL_PROCESS:
        surviving_ids = get_param_process_inputs(sources)
        surviving_ids = surviving_ids - get_process_inputs(sources)
    else:
        all_inputs: list[TInput] = list(tool.inputs_map().values())
        all_ids = get_all_input_ids(all_inputs)
        process_ids = get_process_inputs(sources)
        surviving_ids = all_ids - process_ids
    return surviving_ids

def get_param_inputs_toolmode(sources: dict[str, Any]) -> set[str]:
    """
    nonfile types 
    
    if MINIMAL_PROCESS:
        - remove inputs which are optional
        - remove inputs with defaults
    """
    
    all_inputs: list[TInput] = list(tool.inputs_map().values())
    
    surviving_ids = get_all_input_ids(all_inputs)
    file_ids = get_file_input_ids(all_inputs)
    optional_ids = get_optional_input_ids(all_inputs)
    default_ids = get_default_input_ids(all_inputs)
    
    if settings.MINIMAL_PROCESS:
        surviving_ids = surviving_ids - file_ids
        surviving_ids = surviving_ids - optional_ids
        surviving_ids = surviving_ids - default_ids
    else:
        surviving_ids = surviving_ids - file_ids

    return surviving_ids


### internal inputs

def get_internal_input_ids(tool: CommandTool, sources: dict[str, Any]) -> set[str]:
    """
    internal inputs = all inputs - process inputs - param inputs
    """
    all_inputs: list[TInput] = list(tool.inputs_map().values())

    surviving_ids = get_all_input_ids(all_inputs)
    process_inputs = get_process_inputs(sources)
    param_inputs = get_param_inputs(sources)
    
    surviving_ids = surviving_ids - process_inputs
    surviving_ids = surviving_ids - param_inputs
    return surviving_ids



###  helper methods

# based on TInput info
def get_all_input_ids(tinputs: list[TInput]) -> set[str]:
    return {x.id() for x in tinputs}

def get_file_input_ids(tinputs: list[TInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs which are File types"""
    return {x.id() for x in tinputs if isinstance(nfgen_utils.get_base_type(x.intype), File)}

def get_optional_input_ids(tinputs: list[TInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs which are optional"""
    out: set[str] = set()
    for tinp in tinputs:
        basetype = nfgen_utils.get_base_type(tinp.intype)
        if basetype and basetype.optional:
            out.add(tinp.id())
    return out

def get_default_input_ids(tinputs: list[TInput]) -> set[str]:
    """get tool inputs (ids) for tool inputs with a default value"""
    return {x.id() for x in tinputs if x.default is not None}


# based on step info
def get_channel_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a channel"""
    out: set[str] = set()
    for tag, src in sources.items():
        node = nfgen_utils.resolve_node(src)
        if isinstance(node, InputNode):
            if channels.exists(node.uuid):
                out.add(tag)
    return out

def get_param_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a channel"""
    out: set[str] = set()
    for tag, src in sources.items():
        node = nfgen_utils.resolve_node(src)
        if isinstance(node, InputNode):
            if params.exists(node.uuid):
                out.add(tag)
    return out

def get_connection_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being fed a value from a step connection"""
    out: set[str] = set()
    for tag, src in sources.items():
        node = nfgen_utils.resolve_node(src)
        if isinstance(node, StepNode):
            out.add(tag)
    return out

def get_scatter_process_inputs(sources: dict[str, Any]) -> set[str]:
    """get tool inputs (ids) which are being scattered on"""
    out: set[str] = set()
    for inname, src in sources.items():
        scatter = src.source_map[0].scatter
        node = nfgen_utils.resolve_node(src)
        if scatter and isinstance(node, InputNode):
            out.add(inname)
    return out




### DEPRECATED ###


# def get_connection_process_inputs(sources: dict[str, Any]) -> set[str]:
#     """get tool inputs (ids) which are being fed a value from a step connection"""
#     out: set[str] = set()
#     for inname, src in sources.items():
#         node = nfgen_utils.resolve_node(src)
#         if isinstance(node, StepNode):
#             out.add(inname)
#     return out

# def get_nonfile_wfinp_connected_input_ids(sources: dict[str, Any]) -> set[str]:
#     """get tool inputs (ids) which are being fed a value from a non-File type workflow input"""
#     out: set[str] = set()
#     for inname, src in sources.items():
#         node = nfgen_utils.resolve_node(src) 
#         if isinstance(node, InputNode):
#             if not isinstance(nfgen_utils.get_base_type(node.datatype), File):
#                 out.add(inname)
#     return out











    # def get_src_varname(self, inp: ToolInput) -> str:
    #     # get variable name of feeder (process input or param)
    #     if inp.id() in self.process_inputs:
    #         return inp.id()
    #     elif inp.id() in self.param_inputs:
    #         return f'params.{inp.id()}'
    #     else:
    #         raise RuntimeError


#                 # positionals
#                 case ToolInput(prefix=None): 
#                     if inp.default is not None:
#                         self.handle_positional_default(inp)
#                     elif inp.input_type.optional == True:
#                         self.handle_positional_optional(inp)
#                     elif inp.input_type.optional == False and inp.default is None:
#                         self.handle_positional(inp)
#                     else:
#                         raise RuntimeError
                
#                 # flags
#                 case ToolInput(input_type=Boolean()):
#                     if str(inp.default) == 'False':
#                         self.handle_flag_false_default(inp)
#                     else:
#                         self.handle_flag_true_default(inp)

#                 # options
#                 case ToolInput():
#                     if inp.default is not None:
#                         self.handle_option_default(inp)
#                     elif inp.input_type.optional == True:
#                         self.handle_option_optional(inp)
#                     elif inp.input_type.optional == False and inp.default is None:
#                         self.handle_option(inp)
#                     else:
#                         raise RuntimeError


#     def should_autofill(self, inp: ToolInput) -> bool:
#         if settings.MINIMAL_PROCESS:
#             if not inp.id() in self.process_inputs and not inp.id() in self.param_inputs:
#                 return True
#         return False

#     def handle_positional_default(self, inp: ToolInput) -> None:
#         if self.should_autofill(inp):
#             self.script.append(f'{inp.default}')
#         else:
#             src_name = self.get_src_varname(inp)
#             self.prescript.append(f'def {inp.id()} = {src_name} ? {src_name} : "{inp.default}"')
#             self.script.append(f'${{{inp.id()}}}')
    
#     def handle_positional_optional(self, inp: ToolInput) -> None:
#         if self.should_autofill(inp):
#             # if its not a process input, 
#             # and doesn't have a default,
#             # add nothing to script.
#             pass 
#         else:
#             src_name = self.get_src_varname(inp)
#             self.prescript.append(f'def {inp.id()} = {src_name} ? "${{{src_name}}}" : ""')
#             self.script.append(f'${{{inp.id()}}}')
    
#     def handle_positional(self, inp: ToolInput) -> None:
#         src_name = self.get_src_varname(inp)
#         self.script.append(f'${{{src_name}}}')

#     def handle_flag_false_default(self, inp: ToolInput) -> None:
#         if self.should_autofill(inp):
#             # if its not a process input, 
#             # and its false by default, 
#             # add nothing to script.
#             pass  
#         else:
#             prefix = inp.prefix
#             src_name = self.get_src_varname(inp)
#             self.prescript.append(f'def {inp.id()} = {src_name} ? "{prefix}" : ""')
#             self.script.append(f'${{{src_name}}}')

#     def handle_flag_true_default(self, inp: ToolInput) -> None:
#         prefix = inp.prefix
#         assert(prefix)
#         if self.should_autofill(inp):
#             self.script.append(prefix)
#         else:
#             src_name = self.get_src_varname(inp)
#             self.prescript.append(f'def {inp.id()} = {src_name} == false ? "" : "{prefix}"')
#             self.script.append(f'${{{src_name}}}')

#     def handle_option_default(self, inp: ToolInput) -> None:
#         if self.should_autofill(inp):
#             prefix = self.get_option_prefix(inp)
#             self.script.append(f'{prefix}{inp.default}')
#         else:
#             prefix = self.get_option_prefix(inp)
#             src_name = self.get_src_varname(inp)
#             self.prescript.append(f'def {inp.id()} = {src_name} ? {src_name} : "{inp.default}"')
#             self.script.append(f'{prefix}${{{inp.id()}}}')

#     def handle_option_optional(self, inp: ToolInput) -> None:
#         if self.should_autofill(inp):
#             # if its not a process input, 
#             # and doesn't have a default,
#             # add nothing to script.
#             pass 
#         else:
#             prefix = self.get_option_prefix(inp)
#             src_name = self.get_src_varname(inp)
#             self.prescript.append(f'def {inp.id()} = {src_name} ? "{prefix}${{{src_name}}}" : ""')
#             self.script.append(f'${{{inp.id()}}}')

#     def handle_option(self, inp: ToolInput) -> None:
#         prefix = self.get_option_prefix(inp)
#         src_name = self.get_src_varname(inp)
#         self.script.append(f'{prefix}${{{src_name}}}')

#     def get_option_prefix(self, inp: ToolInput) -> str:
#         assert(inp.prefix)
#         if inp.separate_value_from_prefix == False:
#             return inp.prefix 
#         else:
#             delim = inp.separator if inp.separator else ' '
#             return inp.prefix + delim






# """
# single example


# // these will be overridden by workflow params if exist
# params.opt                         = 500
# params.opt_optional                = 500
# params.opt_default                 = 500
# params.flag_true_default           = true
# params.flag_false_default          = true

# process SINGLE_ITEMS_PARAMS {
#     input:
#       path pos

#     output:
#       stdout
    
#     script:
#       def opt_optional = params.opt_optional ? "--min-length ${params.opt_optional}" : ""
#       def opt_default = params.opt_default ?: '100'
#       def flag_true_default = params.flag_true_default == false ? "" : "--trim-n"
#       def flag_false_default = params.flag_false_default ? "--trim-n" : ""
#       '''
#       echo --- SINGLES ---
#       echo POSITIONAL: ${pos}
#       echo OPTION: --min-length ${params.opt}
#       echo OPTION_OPTIONAL: ${opt_optional}
#       echo OPTION_DEFAULT: --min-length ${opt_default}
#       echo FLAG_TRUE_DEFAULT: ${flag_true_default}
#       echo FLAG_FALSE_DEFAULT: ${flag_false_default}
#       '''
# }

# """



    ### ARRAYS

# can ignore ToolArguments
# can ignore Boolean ToolInputs
"""
### pre-script fmt
PS_NONE     = None
PS_SIMPLE   = 'def {name} = {derived}'
PS_TERNARY  = 'def {name} = {name} ? {derived} : {fallback}'
PS_FLAG_TRUE   = 'def {name} = {name} == false ? {fallback} : "{derived}"'
PS_FLAG_FALSE  = 'def {name} = {name} ? "{derived}" : {fallback}'

### derived expr fmt
DE_DIRECT           = '{src}'  # injected into script
DE_ARRAY            = '{src}.join{separator}'
DE_ARRAY_PREFIXEACH = '{src}.collect{{ {prefix} + it }}.join({separator})'

### prefix fmt
PF_SEPARATE = '{prefix}{separator}'  # separator = inp.separator if inp.separator else ' '
PF_JOINED   = '{prefix}'

### src variable fmt
SC_PROCESS_INPUT   = '{name}'
SC_PARAM           = 'params.{name}'

### script fmt
ST_DIRECT   = '{varname}'
ST_PREFIXED = '{prefix}{varname}'
"""

"""
The following is true for MINIMAL_PROCESS=True and MINIMAL_PROCESS=False.
Each tool input will either have a process input or param, or can be 
autofilled (is optional=ignored, has default=templated)

note: handle booleans seperately to positionals / options

PRE-SCRIPT FORMATS ---
(positionals / options)
basic:          None -> {expr} = {name} or {params.name} so gets injected into script
basic (arr):    def {name} = {expr} 
optional:       def {name} = {name} ? {expr} : ''
default:        def {name} = {name} ? {expr} : {default}

(flags)
true_default    def {name} = {name} == false ? "" : "{prefix}"
false_default   def {name} = {name} ? "{prefix}" : ""

EXPR ---
(doesn't matter if optional or default etc, the {expr} will always be the same)
basic                   {src} -> injected into script
arr               {src}.join{delim}
arr prefixeach    {src}.collect{ {prefix} + it }.join({delim})}

DE_DIRECT           = '{src}'  # injected into script
DE_ARRAY            = '{src}.join{delim}'
DE_ARRAY_PREFIXEACH = '{src}.collect{{ {prefix} + it }}.join({delim})'

SRC ---
process input
param

### prefix fmt
PF_SEPARATE = '{prefix}{separator}'  # separator = inp.separator if inp.separator else ' '
PF_JOINED   = '{prefix}'

### src variable fmt
SC_PROCESS_INPUT   = '{name}'
SC_PARAM           = 'params.{name}'

### script fmt
ST_DIRECT   = '{varname}'
ST_PREFIXED = '{prefix}{varname}'

FINAL VALUE ---
(all = positional / option)
optional (all)          {prefix}{src}
optional (all)(arr)     
default (all)           {src}   # the prefix will be in the script body
default (all)(arr)      

SCRIPT FORMATS ---
positional (all)            {final_value}
positional (all)(arr)       {final_value}

option (basic)              {prefix}{final_value}  # prefix due to mandatory
option (optional)           {final_value}
option (default)            {prefix}{final_value}  # prefix due to always having value

option (basic)(arr)                 {prefix}{final_value}  # prefix due to mandatory
option (optional)(arr)              {final_value}
option (default)(arr)               {prefix}{final_value}  # prefix due to always having value
option (default)(prefixeach)(arr)   {final_value}          

flag (true default)         {final_value}
flag (false default)        {final_value}

"""
"""
The following is true for MINIMAL_PROCESS=True and MINIMAL_PROCESS=False.
Each tool input will either have a process input or param, or can be 
autofilled (is optional=ignored, has default=templated)

note: handle booleans seperately to positionals / options


"""