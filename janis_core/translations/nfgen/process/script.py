
from typing import Any, Tuple, Optional
from janis_core import CommandTool, ToolArgument, ToolInput

from janis_core.translations.nfgen.unwrap import unwrap_expression
from janis_core.translations.nfgen import ordering
from janis_core.translations.nfgen import utils
from janis_core.translations.nfgen import settings

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
        self.process_inputs = utils.get_process_input_ids(tool, self.values)
        self.param_inputs = utils.get_param_input_ids(tool, self.values)
        self.internal_inputs = utils.get_internal_input_ids(tool, self.values)

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

    def get_src_varname(self, inp: ToolInput) -> str:
        # get variable name of feeder (process input or param)
        if inp.id() in self.process_inputs:
            return inp.id()
        elif inp.id() in self.param_inputs:
            return f'params.{inp.id()}'
        else:
            raise RuntimeError

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