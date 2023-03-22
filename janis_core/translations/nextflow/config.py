

from __future__ import annotations
from collections import defaultdict
from janis_core.types import File

from janis_core import settings
from janis_core import translation_utils as utils
from .params import Param, getall
from .casefmt import to_case


BOILERPLATE_LINES = ['docker.enabled = true']
INDENT = settings.translate.nextflow.NF_INDENT
MAX_LINE_WIDTH = 80
COMMENTER = '// '
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
    def linewidth(self) -> int:
        return max([p.width for p in self.params]) + 1

    def to_string(self) -> str:
        heading = f'{INDENT}{COMMENTER}{self.heading}\n'
        out: str = ''
        out += heading
        for param in self.params:
            out += f'{self.param_to_string(param)}\n'
        return out
    
    def param_to_string(self, param: Param) -> str:
        dtype = param.janis_type
        basetype = utils.get_base_type(dtype)
        basetype = utils.ensure_single_type(basetype)

        if utils.is_array_secondary_type(dtype):
            return self.format_param_array_secondary(param)
        
        elif utils.is_secondary_type(basetype):
            return self.format_param_secondary(param)
        
        elif utils.is_array_file_pair_type(dtype):
            return self.format_param_array_file_pair(param)

        elif utils.is_file_pair_type(basetype):
            return self.format_param_file_pair(param)

        elif isinstance(basetype, File) and dtype.is_array():
            return self.format_param_file_array(param)
        
        elif isinstance(basetype, File):
            return self.format_param_file(param)
        
        elif dtype.is_array():
            return self.format_param_val_array(param)
        
        else:
            return self.format_param_val(param)
        

    def format_param_array_secondary(self, param: Param) -> str:
        dtype = param.janis_type
        basetype = utils.get_base_type(dtype)
        exts = utils.get_extensions(basetype, remove_prefix_symbols=True)

        text: str = ''
        text += f'{INDENT}// array of {basetype.name()}\n'
        text += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
        text += f'{INDENT}{INDENT}[\n'
        for ext in exts:
            text += f'{INDENT}{INDENT}{INDENT}// {ext}\n'
        text += f'{INDENT}{INDENT}],\n'
        text += f'{INDENT}]\n'
        print(text)
        return text

    def format_param_secondary(self, param: Param) -> str:
        dtype = param.janis_type
        basetype = utils.get_base_type(dtype)
        exts = utils.get_extensions(basetype, remove_prefix_symbols=True)

        text: str = ''
        text += f'{INDENT}// {basetype.name()}\n'
        text += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
        for ext in exts:
            text += f'{INDENT}{INDENT}// {ext}\n'
        text += f'{INDENT}]\n'
        print(text)
        return text

    def format_param_array_file_pair(self, param: Param) -> str:
        text = f"""\
{INDENT}{param.name:<{self.linewidth}} = [
{INDENT}{INDENT}[
{INDENT}{INDENT}{INDENT}// read 1
{INDENT}{INDENT}{INDENT}// read 2
{INDENT}{INDENT}],
{INDENT}]"""
        return text

    def format_param_file_pair(self, param: Param) -> str:
        text = f"""\
{INDENT}{param.name:<{self.linewidth}} = [
{INDENT}{INDENT}// read 1
{INDENT}{INDENT}// read 2
{INDENT}]"""
        return text
    
    def format_param_file_array(self, param: Param) -> str:
        comment = 'list files here'
        text = f"""\
{INDENT}{param.name:<{self.linewidth}} = [
{INDENT}{INDENT}// {comment}
{INDENT}]"""
        return text
    
    def format_param_file(self, param: Param) -> str:
        return f'{INDENT}{param.name:<{self.linewidth}} = {param.groovy_value}'
    
    def format_param_val_array(self, param: Param) -> str:
        if param.default is not None:
            single_line = f'{INDENT}{param.name:<{self.linewidth}} = {param.groovy_value}'
            
            multi_line = ''
            multi_line += f'{INDENT}{param.name:<{self.linewidth}} = [\n'
            for val in param.groovy_value.strip('[]').split(', '):
                multi_line += f'{INDENT}{INDENT}{val},\n'
            multi_line += f'{INDENT}]\n'
            print(multi_line)
            
            return single_line if len(single_line) <= MAX_LINE_WIDTH else multi_line

        else:
            return f'{INDENT}{param.name:<{self.linewidth}} = []  // list values here'

    def format_param_val(self, param: Param) -> str:
        return f'{INDENT}{param.name:<{self.linewidth}} = {param.groovy_value}'



class ParamGrouper:
    def __init__(self, params: list[Param]):
        self.params = params

    def group(self) -> list[ParamGroup]:
        # awful code but no time. 
        outdir: list[Param] = []
        inputs: list[Param] = []
        processes: list[Param] = []
        subwf: dict[str, list[Param]] = defaultdict(list)
        for param in self.params:
            heading = self.get_group_heading(param)
            if heading == 'OUTPUT DIRECTORY':
                outdir.append(param)
            elif heading == 'INPUTS':
                inputs.append(param)
            elif heading == 'PROCESSES':
                processes.append(param)
            else:
                subwf[heading].append(param)

        # sorting
        out: list[ParamGroup] = []
        if outdir:
            out.append(ParamGroup('OUTPUT DIRECTORY', outdir))
        if inputs:
            out.append(ParamGroup('INPUTS', inputs))
        if processes:
            out.append(ParamGroup('PROCESSES', processes))
        if subwf:
            # alphabetical sorted subworkflow groups
            for heading, param_list in sorted(subwf.items()):
                out.append(ParamGroup(heading, param_list))
        return out

    def get_group_heading(self, param: Param) -> str:
        if param.name == 'outdir':
            name = 'OUTPUT DIRECTORY'
        elif len(param.scope.labels) == 1:
            name = 'INPUTS'
        elif len(param.scope.labels) == 2:
            if param.scope.subtypes[-1] == 'tool':
                name = 'PROCESSES'
            else:
                name = self.get_group_heading_subwf(param)
        elif len(param.scope.labels) > 2:
            name = self.get_group_heading_subwf(param)
        else:
            raise NotImplementedError
        return name

    def get_group_heading_subwf(self, param: Param) -> str:
        subname = param.scope.current_entity # safety net. kinda weird.
        for label, subtype in zip(reversed(param.scope.labels), reversed(param.scope.subtypes)):
            if subtype == 'workflow':
                subname = label
                break
        subname = to_case(subname, settings.translate.nextflow.NF_PROCESS_CASE)
        return f'SUBWORKFLOW: {subname}'
        
