

from typing import Optional
from janis_core.types import File, String
from .params import Param, getall
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
    body = generate_config_body(params)
    body_str = body.to_string()
    defaults_str = defaults_to_string()
    config = TEMPLATE.format(defaults=defaults_str, body=body_str)
    return config



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


class ConfigBody:
    def __init__(self):
        self.outdir: Optional[Param] = None
        self.wfinput_params: list[Param] = []
        self.process_params: list[Param] = []
        self.subworkflow_params: dict[str, list[Param]] = {}
    
    def to_string(self) -> str:
        groups: list[str] = []
        
        if self.outdir:
            group_str = self.group_to_string('OUTPUT DIRECTORY', [self.outdir])
            groups.append(group_str)
        
        if self.wfinput_params:
            group_str = self.group_to_string('INPUTS', self.wfinput_params)
            groups.append(group_str)
        
        if self.process_params:
            group_str = self.group_to_string('PROCESSES', self.process_params)
            groups.append(group_str)
        
        if self.subworkflow_params:
            for subwf_label, group in self.subworkflow_params.items():
                label = f'SUBWORKFLOW: {subwf_label}'
                group_str = self.group_to_string(label, group)
                groups.append(group_str)
        
        return '\n'.join(groups)

    def group_to_string(self, heading: str, group: list[Param]) -> str:
        heading = f'{INDENT}{COMMENTER}{heading}\n'
        width = max([p.width for p in group]) + 1
        
        out: str = ''
        out += heading
        for param in group:
            out += f'{self.param_to_string(param, width)}\n'
        return out
    
    def param_to_string(self, param: Param, group_width: int) -> str:
        dtype = param.dtype
        # # secondary files
        # if isinstance(dtype, File) and dtype.has_secondary_files():
        #     return self.param_to_string_secondary_files(param, group_width)
        # # secondary files array
        # if isinstance(dtype, Array):
        #     if isinstance(dtype.subtype(), File) and dtype.subtype().has_secondary_files():
        #         return self.param_to_string_secondary_files_array(param, group_width)
        # file arrays
        if dtype:
            if dtype.is_array() and isinstance(dtype.subtype(), File):
                return self.param_to_string_multiline(
                    param, group_width, comment='list files here'
                )
            # string arrays
            if dtype.is_array() and isinstance(dtype.subtype(), String):
                return self.param_to_string_multiline(
                    param, group_width, comment='list strings here'
                )
        # everything else
        return self.param_to_string_inline(param, group_width)

    def param_to_string_inline(self, param: Param, group_width: int) -> str:
        return f'{INDENT}{param.name:<{group_width}} = {param.groovy_value}'

    def param_to_string_multiline(self, param: Param, group_width: int, comment: str) -> str:
        lines: list[str] = []
        if isinstance(param.default, list):
            lines.append(f'{param.name:<{group_width}} = [')
            for val in param.default:
                if isinstance(val, str):
                    val = f"'{val}'"
                lines.append(f'{INDENT}{val},')
            lines.append(']')
            lines = [f'{INDENT}{ln}' for ln in lines]
            return '\n'.join(lines)
        else:
            lines.append(f'{param.name:<{group_width}} = [')
            lines.append(f'{INDENT}// {comment}')
            lines.append(']')
            lines = [f'{INDENT}{ln}' for ln in lines]
            return '\n'.join(lines)


    

    
    

def generate_config_body(params: list[Param]) -> ConfigBody:
    cbody = ConfigBody()
    for p in params:
        if p.name == 'outdir':
            cbody.outdir = p
        elif p.is_subworkflow_param:
            # TODO
            pass
        elif p.is_channel_input:
            cbody.wfinput_params.append(p)
        else:
            cbody.process_params.append(p)
    return cbody

def defaults_to_string() -> str:
    return '\n'.join(DEFAULT_LINES)

# @dataclass
# class ConfigGroup:
#     params: list[Param]
#     heading: Optional[str]=None

#     @property
#     def width(self) -> int:
#         return max([p.width for p in self.params]) 

#     def param_to_string(self, param: Param) -> str:
#         return f'{param.name:<{self.width}} = {param.value}\n'
    
#     def to_string(self, indent: str) -> str:
#         out: str = ''
#         if self.heading:
#             out += f'{indent}{COMMENTER}{self.heading}\n'
#         for param in self.params:
#             out += f'{indent}{self.param_to_string(param)}'
#         return out

# def generate_blocks(groups: dict[str, list[Param]]) -> list[ConfigGroup]:
#     blocks: list[ConfigBlock] = []
#     for heading, params in groups.items():
#         new_block = ConfigBlock(params, heading)
#         blocks.append(new_block)
#     return blocks

# def sort_groups(groups: dict[str, list[Param]]) -> list[Tuple[str, list[Param]]]:
#     out: list[Tuple[str, list[Param]]] = []
#     for heading, params in groups.items():
#         if 
#     return out

# def blocks_to_string(blocks: list[ConfigGroup]) -> str:
#     out: str = ''
#     for block in blocks:
#         out += f'{block.to_string(INDENT)}\n'
#     out = out.rstrip('\n')
#     return out

    
    
