

from janis_core import Array
from .types import Bam, BamBai, Cram, CramCrai


# function that takes the inputs of a tool and changes every bam type into the respecitve cram type
# while also keeping anything else the same
def cast_input_bams_to_crams(inputs):
    from copy import deepcopy

    retval = deepcopy(inputs)
    for inp in retval:

        # we need to store it the input was optional originally
        is_optional = inp.input_type.optional

        # we need to check for BamBai first, as due to inheritance, the bambai is also a bam
        if isinstance(inp.input_type, BamBai):
            inp.input_type = CramCrai(optional=is_optional)
        elif isinstance(inp.input_type, Bam):
            inp.input_type = Cram(optional=is_optional)
        elif isinstance(inp.input_type, Array):
            internal = None
            if isinstance(inp.input_type.subtype(), BamBai):
                internal = CramCrai
            elif isinstance(inp.input_type.subtype(), Bam):
                internal = Cram

            if internal is not None:
                inp.input_type = Array(
                    internal(optional=inp.input_type.subtype().optional),
                    optional=is_optional,
                )

        # and now that we changed things, we set the optional state again
        inp.input_type.optional = is_optional

    return retval
