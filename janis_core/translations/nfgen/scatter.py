
from . import settings

INDENT = settings.NEXTFLOW_INDENT
CROSS_CHANNEL_NAME = 'ch_cartesian_cross'

"""
ch_bams
.combine(ch_bais)
.multiMap { it ->
    bams: it[0]
    bais: it[1]
}
.set { ch_cartesian_cross }
"""

def cartesian_cross_subname(channel: str) -> str:
    # eg: ch_bams -> ch_cartesian_cross.bams
    ch_subname = channel.replace('ch_', '')
    return f'{CROSS_CHANNEL_NAME}.{ch_subname}'   

def cartesian_cross_operation(channels: list[str]) -> str:
    # .flatten()?????

    lines: list[str] = []
    lines.append(channels[0])                   # ch_bams
    
    # cartesian cross
    for ch_name in channels[1:]:
        lines.append(f'.combine({ch_name})')    # .combine(ch_bais)
    
    # multiMap separation
    lines.append('.multiMap { it ->')           # .multiMap { it ->
    for i, ch_name in enumerate(channels):
        ch_subname = ch_subname.replace('ch_', '')     
        lines.append(f'{INDENT}{ch_subname}: it[{i}]')     # bams: it[0] \n bais: it[1]
    lines.append('}')

    # set new channel names
    lines.append(f'.set {{ {CROSS_CHANNEL_NAME} }}')
    text = '\n'.join(lines)
    return text



