

from collections import defaultdict
from typing import Optional
from .params import Param, getall


DEFAULT_LINES = ['docker.enabled = true']
INDENT = ' ' * 4
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


class ConfigBody:
    def __init__(self):
        self.outdir: Optional[Param] = None
        self.wfinput_params: list[Param] = []
        self.process_params: dict[str, list[Param]] = defaultdict(list)
    
    def to_string(self) -> str:
        groups: list[str] = []
        
        if self.outdir:
            group_str = self.group_to_string('OUTPUT DIRECTORY', [self.outdir])
            groups.append(group_str)
        
        if self.wfinput_params:
            group_str = self.group_to_string('WORKFLOW INPUTS', self.wfinput_params)
            groups.append(group_str)
        
        for scope_name, group in self.process_params.items():
            group_str = self.group_to_string(scope_name, group)
            groups.append(group_str)

        return '\n'.join(groups)

    def group_to_string(self, heading: str, group: list[Param]) -> str:
        heading = f'{INDENT}{COMMENTER}{heading}\n'
        width = max([p.width for p in group]) + 1
        
        out: str = ''
        out += heading
        for param in group:
            out += f'{INDENT}{self.param_to_string(param, width)}\n'
        return out
    
    def param_to_string(self, param: Param, group_width: int) -> str:
        return f'{param.name:<{group_width}} = {param.default}'
    

def generate_config_body(params: list[Param]) -> ConfigBody:
    cbody = ConfigBody()
    for p in params:
        if p.name == 'outdir':
            cbody.outdir = p
        elif p.is_wf_input:
            cbody.wfinput_params.append(p)
        elif p.scope:
            scope_name = f'{"_".join(p.scope)}'
            cbody.process_params[scope_name].append(p)
        else:
            raise NotImplementedError
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
#         return f'{param.name:<{self.width}} = {param.default}\n'
    
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

    
    
