

from __future__ import annotations
from collections import defaultdict
from typing import Optional

from janis_core.types import File
from janis_core import settings
from janis_core import translation_utils as utils

from ..params import Param, getall
from ..casefmt import to_case
# from .. import nulls

USE_CLOSED_FORM = True
INDENT = settings.translate.nextflow.NF_INDENT
MAX_LINE_WIDTH = 80

TEMPLATE = """\

nextflow.enable.dsl = 2
docker.enabled = true

params {{
    
    // Placeholder for null values.
    // Do not alter unless you know what you are doing.
    NULL_VALUE = 'NULL'

    // WORKFLOW OUTPUT DIRECTORY
    outdir  = './outputs'

{param_block}
}}
"""

def generate_config() -> str:

    # params    
    param_block_str = ''
    params = getall()
    groups = ParamGrouper(params).group()
    for g in groups:
        param_block_str += f'{g.to_string()}\n'

    # final config text formatting
    config = TEMPLATE.format(param_block=param_block_str)
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
        non_default_widths = [len(value(p)) for p in self.params if p.default is None]
        default_widths = [len(value(p)) for p in self.params if p.default is not None]
        if non_default_widths:
            return max(non_default_widths) + 1
        else:
            return max(default_widths) + 1
    
    @property
    def dtype_width(self) -> int:
        non_default_widths = [len(dtype_label(p)) for p in self.params if p.default is None]
        default_widths = [len(dtype_label(p)) for p in self.params if p.default is not None]
        if non_default_widths:
            return max(non_default_widths) + 1
        else:
            return max(default_widths) + 1

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
    return param.groovy_value
    # dtype = param.dtype
    # # has default
    # if param.default is not None:
    #     val = param.groovy_value
    # else:
    #     val = nulls.get_null_value(dtype)
    # return val

def format_label(param: Param) -> str:
    dtype = param.dtype
    basetype = utils.get_base_type(dtype)
    
    # secondary array
    if utils.is_secondary_array_type(dtype):
        exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
        exts_str = ', '.join(exts)
        label = f'[[{exts_str}]]'
    
    # secondary
    elif utils.is_secondary_type(dtype):
        exts = utils.get_extensions(basetype, remove_prefix_symbols=True)
        exts_str = ', '.join(exts)
        label = f'[{exts_str}]'
    
    # file pair array
    elif utils.is_file_pair_array_type(dtype):
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
    # get datatype name
    name = param.dtype.name().lower()
    if name == 'file':
        dtype_str = 'generic file'
    else:
        dtype_str = param.dtype.name().lower()

    # return label formatted with optionality
    if param.dtype.optional:
        return f'(optional {dtype_str})'
    else:
        return f'(MANDATORY {dtype_str})'



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
        if self.param.default is not None:
            return self.to_string_without_labels()
        else:
            return self.to_string_with_labels()

    def to_string_without_labels(self) -> str:
        return f'{INDENT}{self.name:<{self.name_width}} = {self.value}'

    def to_string_with_labels(self) -> str:
        return f"""\
{INDENT}\
{self.name:<{self.name_width}} \
= \
{self.value:<{self.value_width}} \
// \
{self.dtype_label:<{self.dtype_width}} \
{self.format_label}\
"""



class ParamGrouper:
    def __init__(self, params: list[Param]):
        self.params = params

    def group(self) -> list[ParamGroup]:
        # awful code but no time. 
        mandatory_inputs: list[Param] = []
        optional_inputs: list[Param] = []
        subtools: dict[str, list[Param]] = defaultdict(list)
        subworkflows: dict[str, list[Param]] = defaultdict(list)
        for param in self.params:
            heading = self.get_group_heading(param)
            if heading == 'INPUTS (MANDATORY)':
                mandatory_inputs.append(param)
            elif heading == 'INPUTS (OPTIONAL)':
                optional_inputs.append(param)
            elif heading.startswith('PROCESS:'):
                subtools[heading].append(param)
            else:
                subworkflows[heading].append(param)

        # ordering
        out: list[ParamGroup] = []
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
        # else:
        #     print()
        return name

 