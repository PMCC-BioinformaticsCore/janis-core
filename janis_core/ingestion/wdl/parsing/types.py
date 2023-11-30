

from typing import Optional

import WDL
import janis_core as j
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory

def parse_type(t: WDL.Type.Base, uuid: Optional[str]=None):
    optional = t.optional
    if isinstance(t, WDL.Type.Int):
        return j.Int(optional=optional)
    elif isinstance(t, WDL.Type.String):
        return j.String(optional=optional)
    elif isinstance(t, WDL.Type.Float):
        return j.Float(optional=optional)
    elif isinstance(t, WDL.Type.Boolean):
        return j.Boolean(optional=optional)
    elif isinstance(t, WDL.Type.File):
        return j.File(optional=optional)
    elif isinstance(t, WDL.Type.Directory):
        return j.Directory(optional=optional)
    elif isinstance(t, WDL.Type.Array):
        return j.Array(parse_type(t.item_type, uuid), optional=optional)
    elif isinstance(t, WDL.Type.StructInstance):
        if not uuid:
            raise Exception('add uuid here')
        msg = 'WDL Struct type unsupported. Has been cast to File type.'
        log_message(uuid, msg, ErrorCategory.DATATYPES)
        return j.File(optional=optional)

    raise Exception(f"Didn't handle WDL type conversion for '{t}' ({type(t)})")