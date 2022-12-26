

"""

docker.enabled = true

params {

    // INPUTS
    fastqs          = []
    reference_fasta = null
    reference_amb   = null

    // PROCESSES
    sample_name           = null
    allele_freq_threshold = 0.05
    min_mapping_qual      = null

    // SUBWORKFLOW: ALIGN_AND_SORT
    align_and_sort.bwamem.mark_shorter_splits    = true
    align_and_sort.cutadapt.quality_cutoff       = 15
    align_and_sort.cutadapt.minimum_length       = 50
    align_and_sort.sortsam.sort_order            = 'coordinate'
    align_and_sort.sortsam.create_index          = true
    align_and_sort.sortsam.validation_stringency = 'SILENT'
    align_and_sort.sortsam.max_records_in_ram    = 5000000

    // SUBWORKFLOW: MERGE_AND_MARKDUPS
    merge_and_markdups.create_index                                          = true
    merge_and_markdups.max_records_in_ram                                    = 5000000
    merge_and_markdups.merge_sam_files.use_threading                         = true
    merge_and_markdups.merge_sam_files.validation_stringency                 = 'SILENT'

"""

from __future__ import annotations
from collections import defaultdict
from typing import Optional
from janis_core.types import File, String
from .params import Param, getall
from .casefmt import to_case
from . import settings

DEFAULT_LINES = ['docker.enabled = true']
INDENT = settings.NF_INDENT
COMMENTER = '// '
TEMPLATE = """\
{defaults}

params {{

{body}
}}
"""

def generate_config() -> str:
    params = getall()
    groups = ParamGrouper(params).group()
    body_str = groups_to_string(groups)
    defaults_str = defaults_to_string()
    config = TEMPLATE.format(defaults=defaults_str, body=body_str)
    return config

def groups_to_string(groups: list[ParamGroup]) -> str:
    group_strs: list[str] = [g.to_string() for g in groups]
    return '\n'.join(group_strs) 

def defaults_to_string() -> str:
    return '\n'.join(DEFAULT_LINES)


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
        comment: Optional[str] = None
        if dtype:
            if dtype.is_array() and isinstance(dtype.subtype(), File):
                comment = 'list files here'
            if dtype.is_array() and isinstance(dtype.subtype(), String):
                comment = 'list strings here'
        # everything else
        return self.param_to_string_inline(param, comment)

    def param_to_string_inline(self, param: Param, comment: Optional[str]=None) -> str:
        line = f'{INDENT}{param.name:<{self.linewidth}} = {param.groovy_value}'
        if comment:
            line += f'{INDENT}// {comment}'
        return line

    def param_to_string_multiline(self, param: Param, comment: str) -> str:
        lines: list[str] = []
        if isinstance(param.default, list):
            lines.append(f'{param.name:<{self.linewidth}} = [')
            for val in param.default:
                if isinstance(val, str):
                    val = f"'{val}'"
                lines.append(f'{INDENT}{val},')
            lines.append(']')
            lines = [f'{INDENT}{ln}' for ln in lines]
            return '\n'.join(lines)
        else:
            lines.append(f'{param.name:<{self.linewidth}} = [')
            lines.append(f'{INDENT}// {comment}')
            lines.append(']')
            lines = [f'{INDENT}{ln}' for ln in lines]
            return '\n'.join(lines)



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
        
