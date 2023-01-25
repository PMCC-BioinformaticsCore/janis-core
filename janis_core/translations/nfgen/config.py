

from __future__ import annotations
from collections import defaultdict
from janis_core.types import File

from . import settings
from . import nfgen_utils
from .params import Param, getall
from .casefmt import to_case


BOILERPLATE_LINES = ['docker.enabled = true']
INDENT = settings.NF_INDENT
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

        if nfgen_utils.is_array_secondary_type(dtype):
            return self.format_param_array_secondary(param)
        
        elif nfgen_utils.is_secondary_type(dtype):
            return self.format_param_secondary(param)
        
        elif nfgen_utils.is_array_file_pair_type(dtype):
            return self.format_param_array_file_pair(param)

        elif nfgen_utils.is_file_pair_type(dtype):
            return self.format_param_file_pair(param)

        elif dtype.is_array():
            return self.format_param_array(param)
        
        else:
            return self.format_param_single(param)


    def format_param_array_secondary(self, param: Param) -> str:
        dtype = param.janis_type
        basetype = nfgen_utils.get_base_type(dtype)
        exts = nfgen_utils.get_extensions(basetype, remove_symbols=True)

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
        basetype = nfgen_utils.get_base_type(dtype)
        exts = nfgen_utils.get_extensions(basetype, remove_symbols=True)

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
    
    def format_param_array(self, param: Param) -> str:
        basetype = nfgen_utils.get_base_type(param.janis_type)
        if isinstance(basetype, File):
            comment = 'list files here'
        else:
            comment = 'list values here'
        text = f"""\
{INDENT}{param.name:<{self.linewidth}} = [
{INDENT}{INDENT}// {comment}
{INDENT}]"""
        print(text)
        return text

    def format_param_single(self, param: Param) -> str:
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
        subname = param.scope.labels[-1] # safety net. kinda weird.
        for label, subtype in zip(reversed(param.scope.labels), reversed(param.scope.subtypes)):
            if subtype == 'workflow':
                subname = label
                break
        subname = to_case(subname, settings.NF_PROCESS_CASE)
        return f'SUBWORKFLOW: {subname}'
        
