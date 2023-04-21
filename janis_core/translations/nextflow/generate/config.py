

from __future__ import annotations
from collections import defaultdict
from typing import Optional

from janis_core.types import File
from janis_core import settings
from janis_core import translation_utils as utils

from ..params import Param, getall
from ..casefmt import to_case

USE_CLOSED_FORM = True
BOILERPLATE_LINES = ['docker.enabled = true']
INDENT = settings.translate.nextflow.NF_INDENT
MAX_LINE_WIDTH = 80
TEMPLATE = """\
{boilerplate}

params {{

{param_block}
}}
"""

def generate_config() -> str:
    # boilerplate 
    boilerplate_str = '\n'.join(BOILERPLATE_LINES)

    # params    
    param_block_str = ''
    params = getall()
    groups = ParamGrouper(params).group()
    for g in groups:
        param_block_str += f'{g.to_string()}\n'

    # final config text formatting
    config = TEMPLATE.format(boilerplate=boilerplate_str, param_block=param_block_str)
    return config



class ParamGroup:
    def __init__(self, heading: str, params: list[Param]) -> None:
        self.heading = heading
        self.params = params

    @property
    def name_width(self) -> int:
        return max([len(name(p)) for p in self.params]) + 1
    
    @property
    def value_width(self) -> int:
        return max([len(value(p)) for p in self.params]) + 1
    
    @property
    def dtype_width(self) -> int:
        return max([len(dtype_label(p)) for p in self.params]) + 1

    def to_string(self) -> str:
        heading = f'{INDENT}// {self.heading}\n'
        out: str = ''
        out += heading
        for param in self.params:
            formatter = ParamFormatterClosedForm(param, self.name_width, self.value_width, self.dtype_width)
            out += f'{formatter.to_string()}\n'
        return out



# helper functions to format params

def name(param: Param) -> str:
    return param.name

def value(param: Param) -> str:
    dtype = param.dtype
    basetype = utils.get_base_type(dtype)
    
    # secondary array optional
    if utils.is_array_secondary_type(dtype) and dtype.optional:
        exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
        val = ["'NO_FILE'"] * len(exts)
        val = ', '.join(val)
        val = f'[[{val}]]'
    
    # secondary array
    elif utils.is_array_secondary_type(dtype):
        val = '[[]]'
    
    # secondary optional
    elif utils.is_secondary_type(dtype) and dtype.optional:
        exts = utils.get_extensions(dtype, remove_prefix_symbols=True)
        val = ["'NO_FILE'"] * len(exts)
        val = ', '.join(val)
        val = f'[{val}]'
    
    # secondary
    elif utils.is_secondary_type(basetype):
        val = '[]'
    
    # file pair array optional
    elif utils.is_array_file_pair_type(dtype) and dtype.optional:
        val = "[['NO_FILE', 'NO_FILE']]"
    
    # file pair array
    elif utils.is_array_file_pair_type(dtype):
        val = '[[]]'

    # file pair optional
    elif utils.is_file_pair_type(basetype) and dtype.optional:
        val = "['NO_FILE', 'NO_FILE']"

    # file pair
    elif utils.is_file_pair_type(basetype):
        val = '[]'

    # file array optional
    elif isinstance(basetype, File) and dtype.is_array() and dtype.optional:
        val = "['NO_FILE]"
    
    # file array
    elif isinstance(basetype, File) and dtype.is_array():
        val = '[]'
    
    # file optional
    elif isinstance(basetype, File) and dtype.optional:
        val = "'NO_FILE'"
    
    # file 
    elif isinstance(basetype, File):
        val = 'null'
    
    # array default
    elif dtype.is_array() and param.default is not None:
        val = param.groovy_value
    
    # array
    elif dtype.is_array():
        val = '[]'
    
    # default
    elif param.default is not None:
        val = param.groovy_value

    # generic
    else:
        val = 'null'
        
    return val

def format_label(param: Param) -> str:
    dtype = param.dtype
    basetype = utils.get_base_type(dtype)
    
    # secondary array
    if utils.is_array_secondary_type(dtype):
        exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
        exts_str = ', '.join(exts)
        label = f'[[{exts_str}]]'
    
    # secondary
    elif utils.is_secondary_type(dtype):
        exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
        exts_str = ', '.join(exts)
        label = f'[{exts_str}]'
    
    # file pair array
    elif utils.is_array_file_pair_type(dtype):
        label = "[[pair1, pair2]]"
    
    # file pair
    elif utils.is_file_pair_type(basetype):
        label = "[pair1, pair2]"

    # file array
    elif isinstance(basetype, File) and dtype.is_array():
        name = basetype.name().lower()
        label = f'[{name}1, ...]'
    
    # file 
    elif isinstance(basetype, File):
        label = None
    
    # array
    elif dtype.is_array():
        name = basetype.name().lower()
        label = f'[{name}1, ...]'

    # generic
    else:
        label = None
    
    if label is None:
        return ''
    
    return f'eg. {label}'

def dtype_label(param: Param) -> str:
    name = param.dtype.name().lower()
    if name == 'file':
        name = 'generic file'
    return f'({name})'



class ParamFormatterClosedForm:
    def __init__(self, param: Param, name_width: int, value_width: int, dtype_width: int) -> None:
        self.param = param
        self.name_width = name_width
        self.value_width = value_width
        self.dtype_width = dtype_width

    @property
    def name(self) -> str:
        return name(self.param)
    
    @property
    def value(self) -> str:
        return value(self.param)
    
    @property
    def format_label(self) -> Optional[str]:
        return format_label(self.param)

    @property
    def dtype_label(self) -> Optional[str]:
        return dtype_label(self.param)

    def to_string(self) -> str:
        return f"""\
{INDENT}\
{self.name:<{self.name_width}} \
= \
{self.value:<{self.value_width}}  \
// \
{self.dtype_label:<{self.dtype_width}} \
{self.format_label}\
"""



class ParamGrouper:
    def __init__(self, params: list[Param]):
        self.params = params

    def group(self) -> list[ParamGroup]:
        # awful code but no time. 
        outdir: list[Param] = []
        mandatory_inputs: list[Param] = []
        optional_inputs: list[Param] = []
        subtools: dict[str, list[Param]] = defaultdict(list)
        subworkflows: dict[str, list[Param]] = defaultdict(list)
        for param in self.params:
            heading = self.get_group_heading(param)
            if heading == 'OUTPUT DIRECTORY':
                outdir.append(param)
            elif heading == 'INPUTS (MANDATORY)':
                mandatory_inputs.append(param)
            elif heading == 'INPUTS (OPTIONAL)':
                optional_inputs.append(param)
            elif heading.startswith('PROCESS:'):
                subtools[heading].append(param)
            else:
                subworkflows[heading].append(param)

        # sorting
        out: list[ParamGroup] = []
        if outdir:
            out.append(ParamGroup('OUTPUT DIRECTORY', outdir))
        if mandatory_inputs:
            out.append(ParamGroup('INPUTS (MANDATORY)', mandatory_inputs))
        if optional_inputs:
            out.append(ParamGroup('INPUTS (OPTIONAL)', optional_inputs))
        if subtools:
            # alphabetical sorted tools
            for heading, param_list in sorted(subtools.items()):
                out.append(ParamGroup(heading, param_list))
        if subworkflows:
            # alphabetical sorted subworkflows
            for heading, param_list in sorted(subworkflows.items()):
                out.append(ParamGroup(heading, param_list))
        return out

    def get_group_heading(self, param: Param) -> str:
        if param.name == 'outdir':
            name = 'OUTPUT DIRECTORY'
        elif param.subtype == 'main_workflow':
            if param.dtype.optional:
                name = 'INPUTS (OPTIONAL)'
            else:
                name = 'INPUTS (MANDATORY)'
        elif param.subtype == 'sub_workflow':
            name = f'SUBWORKFLOW: {to_case(param.task_id, settings.translate.nextflow.NF_PROCESS_CASE)}'
        elif param.subtype == 'sub_tool':
            name = f'PROCESS: {to_case(param.task_id, settings.translate.nextflow.NF_PROCESS_CASE)}'
        else:
            print()
        return name

        



    # def param_to_string(self, param: Param) -> str:
    #     dtype = param.dtype
    #     basetype = utils.get_base_type(dtype)
    #     basetype = utils.ensure_single_type(basetype)

    #     if utils.is_array_secondary_type(dtype):
    #         return self.format_param_array_secondary(param)
        
    #     elif utils.is_secondary_type(basetype):
    #         return self.format_param_secondary(param)
        
    #     elif utils.is_array_file_pair_type(dtype):
    #         return self.format_param_array_file_pair(param)

    #     elif utils.is_file_pair_type(basetype):
    #         return self.format_param_file_pair(param)

    #     elif isinstance(basetype, File) and dtype.is_array():
    #         return self.format_param_file_array(param)
        
    #     elif isinstance(basetype, File):
    #         return self.format_param_file(param)
        
    #     elif dtype.is_array():
    #         return self.format_param_val_array(param)
        
    #     else:
    #         return self.format_param_val(param)
        
    # def get_dtype_label(self, param: Param) -> str:
    #     if param.dtype.name() == 'File':
    #         return 'generic file'
    #     else:
    #         return param.dtype.name().lower()

    # def format_param_array_secondary(self, param: Param) -> str:
    #     dtype = param.dtype
    #     basetype = utils.get_base_type(dtype)
    #     exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
    #     exts = [x.replace("\'", "") for x in exts]
    #     dtype_label = self.get_dtype_label(param)
    #     format_label = self.get_format_label(param)

    #     text: str = ''
    #     # type heading
    #     # text += f'{INDENT}// array of {basetype.name()}\n'

    #     if USE_CLOSED_FORM:
    #         if dtype.optional:
    #             val = ["'NO_FILE'"] * len(exts)
    #             val = ', '.join(val)
    #             val = f'[[{val}]]'
    #         else:
    #             exts_str = ', '.join(exts)
    #             val = f'[[]]  // format: [[{exts_str}]]'
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = {val} ({dtype_label})'

    #     else:
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
    #         text += f'{INDENT}{INDENT}[\n'
    #         for ext in exts:
    #             val = "'NO_FILE'," if dtype.optional else ext
    #             text += f'{INDENT}{INDENT}{INDENT}// {val}\n'
    #         text += f'{INDENT}{INDENT}],\n'
    #         text += f'{INDENT}]'

    #     return text

    # def format_param_secondary(self, param: Param) -> str:
    #     dtype = param.dtype
    #     basetype = utils.get_base_type(dtype)
    #     exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
    #     exts = [x.replace("\'", "") for x in exts]
        
    #     text: str = ''
    #     # type heading
    #     # text += f'{INDENT}// {basetype.name()}\n'
    #     if USE_CLOSED_FORM:
    #         if dtype.optional:
    #             val = ["'NO_FILE'"] * len(exts)
    #             val = ', '.join(val)
    #             val = f'[{val}]'
    #             print()
    #         else:
    #             exts_str = ', '.join(exts)
    #             val = f'[]  // format: [{exts_str}]'
    #             print()
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = {val}'
        
    #     else:
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
    #         for ext in exts:
    #             val = "'NO_FILE'," if dtype.optional else ext
    #             text += f'{INDENT}{INDENT}{INDENT}// {val}\n'
    #         text += f'{INDENT}]'
        
    #     return text

    # def format_param_array_file_pair(self, param: Param) -> str:
    #     dtype = param.dtype
    #     text: str = ''

    #     if USE_CLOSED_FORM:
    #         val = "[['NO_FILE', 'NO_FILE']]" if dtype.optional else '[[]]  // format: [pair1, pair2]'
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = {val}'
        
    #     else:
    #         val1 = "'NO_FILE'," if dtype.optional else '// read1'
    #         val2 = "'NO_FILE'" if dtype.optional else '// read2'
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
    #         text += f'{INDENT}{INDENT}[\n'
    #         text += f'{INDENT}{INDENT}{INDENT}{val1},\n'
    #         text += f'{INDENT}{INDENT}{INDENT}{val2}\n'
    #         text += f'{INDENT}{INDENT}],\n'
    #         text += f'{INDENT}]'
        
    #     return text

    # def format_param_file_pair(self, param: Param) -> str:
    #     dtype = param.dtype
    #     text: str = ''

        
    #     if USE_CLOSED_FORM:
    #         val = "['NO_FILE', 'NO_FILE']" if dtype.optional else '[]  // format: [pair1, pair2]'
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = {val}'
         
    #     else:
    #         val1 = "'NO_FILE'," if dtype.optional else '// read1'
    #         val2 = "'NO_FILE'" if dtype.optional else '// read2'
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
    #         text += f'{INDENT}{INDENT}{val1},\n'
    #         text += f'{INDENT}{INDENT}{val2}\n'
    #         text += f'{INDENT}]'
        
    #     return text
    
    # def format_param_file_array(self, param: Param) -> str:
    #     dtype = param.dtype
    #     text: str = ''
    #     basetype = utils.get_base_type(dtype)
    #     assert(isinstance(basetype, File))
    #     if basetype.name() == 'File':
    #         example = 'example.ext'
    #     else:
    #         example = f'example.{basetype.extension}'

    #     if USE_CLOSED_FORM:
    #         val = "['NO_FILE']" if dtype.optional else f'[]  // format: [{example}]'
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = {val}'
         
    #     else:
    #         val = "'NO_FILE'" if dtype.optional else f'// format: [{example}]'
    #         text += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
    #         text += f'{INDENT}{INDENT}{val},\n'
    #         text += f'{INDENT}]'

    #     return text
    
    # def format_param_file(self, param: Param) -> str:
    #     dtype = param.dtype
    #     val = "'NO_FILE'" if dtype.optional else 'null'
    #     return f'{INDENT}{param.name:<{self.linewidth}} = {val}'
    
    # def format_param_val_array(self, param: Param) -> str:
    #     if param.default is not None:
    #         single_line = f'{INDENT}{param.name:<{self.linewidth}} = {param.groovy_value}'
            
    #         multi_line = ''
    #         multi_line += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
    #         for val in param.groovy_value.strip('[]').split(', '):
    #             multi_line += f'{INDENT}{INDENT}{val},\n'
    #         multi_line += f'{INDENT}]\n'
    #         print(multi_line)
            
    #         return single_line if len(single_line) <= MAX_LINE_WIDTH else multi_line

    #     else:
    #         return f'{INDENT}{param.name:<{self.linewidth}} = []  // list values here'

    # def format_param_val(self, param: Param) -> str:
    #     return f'{INDENT}{param.name:<{self.linewidth}} = {param.groovy_value}'


