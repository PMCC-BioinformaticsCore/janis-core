


from janis_core.types import DataType

from .. import settings

INDENT = settings.NF_INDENT



def is_scatter_relationship(srcscatter: bool, destscatter: bool) -> bool:
    raise NotImplementedError

def handle_scatter_relationship(srcscatter: bool, destscatter: bool, srctype: DataType, desttype: DataType) -> bool:
    raise NotImplementedError




def is_array_singular_mismatch():
    # TODO implement
    raise NotImplementedError

def handle_array_singular_mismatch():
    raise NotImplementedError



"""

which two entities are involved? 
what is scatter doing? 
are the datatypes singular vs array? 
are the datatypes both same type? does this work with scatter?

CONVERSIONS 

# Array 
Array(Singular) -> Array(Singular)
Array(Singular) -> Singular
Singular -> Array(Singular)

# Secondary Array
Array(Secondary) -> Array(Secondary)
Array(Secondary) -> Secondary
Secondary -> Array(Secondary)

# Scatter Singular
Scatter(Singular) -> Scatter(Singular)
Scatter(Singular) -> Singular
Singular -> Scatter(Singular)

# Scatter Secondary
Scatter(Secondary) -> Scatter(Secondary)
Scatter(Secondary) -> Secondary
Secondary -> Scatter(Secondary)

# Scatter Singular Array
Scatter(Array(Singular)) -> Scatter(Array(Singular)) 
Scatter(Array(Singular)) -> Array(Singular)
Array(Singular) -> Scatter(Array(Singular)) 

# Scatter Secondary Array
Scatter(Array(Secondary)) -> Scatter(Array(Secondary)) 
Scatter(Array(Secondary)) -> Array(Secondary)
Scatter(Array(Secondary)) -> Secondary
Array(Secondary) -> Scatter(Array(Secondary))


"""



# SecondaryType -> Array(SecondaryType)





# SCATTER: CARTESIAN CROSS

CROSS_CHANNEL_NAME = 'ch_cartesian_cross'

def cartesian_cross_subname(channel: str) -> str:
    # eg: ch_bams -> ch_cartesian_cross.bams
    ch_subname = channel.replace('ch_', '')
    return f'{CROSS_CHANNEL_NAME}.{ch_subname}'   

def cartesian_cross_operation(channels: list[str]) -> str:
    """
    * just an example, you wouldnt actually do this cuz it mixes up the bams and bais

    - janis - 
    scatter=ScatterDescription(
        ["bams", "bais"],
        method=ScatterMethod.cross
    ),

    - nextflow - 
    ch_bams
    .combine(ch_bais)
    .multiMap { it ->
        bams: it[0]
        bais: it[1]
    }
    .set { ch_cartesian_cross }
    """
    # .flatten()?????

    lines: list[str] = []
    lines.append(channels[0])                   # ch_bams
    
    # cartesian cross
    for ch_name in channels[1:]:
        lines.append(f'.combine({ch_name})')    # .combine(ch_bais)
    
    # multiMap separation
    lines.append('.multiMap { it ->')           # .multiMap { it ->
    for i, ch_name in enumerate(channels):
        ch_subname = ch_name.replace('ch_', '')     
        lines.append(f'{INDENT}{ch_subname}: it[{i}]')     # bams: it[0] \n bais: it[1]
    lines.append('}')

    # set new channel names
    lines.append(f'.set {{ {CROSS_CHANNEL_NAME} }}')
    text = '\n'.join(lines)
    return text






