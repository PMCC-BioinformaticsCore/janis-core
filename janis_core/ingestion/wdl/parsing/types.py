


import WDL
import janis_core as j
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory


def parse_type(wdl_type: WDL.Type.Base, wdl_tool: WDL.Tree.Workflow | WDL.Tree.Task, uuid: str):
    parser = WDLTypeParser(wdl_tool, uuid)
    return parser.parse(wdl_type)


class WDLTypeParser:
    
    def __init__(self, wdl_tool: WDL.Tree.Workflow | WDL.Tree.Task, uuid: str):
        self.wdl_tool = wdl_tool # wdl tool or workflow the type appears in
        self.uuid = uuid # uuid of respective janis entity (logging)
        
    def parse(self, t: WDL.Type.Base) -> j.Type:
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
            return j.Array(self.parse(t.item_type), optional=optional)
        elif isinstance(t, WDL.Type.StructInstance):
            if not self.uuid:
                raise Exception('add self.uuid here')
            msg = 'WDL Struct type unsupported. Has been cast to File type.'
            log_message(self.uuid, msg, ErrorCategory.DATATYPES)
            return j.File(optional=optional)
        raise Exception(f"Didn't handle WDL type conversion for '{t}' ({type(t)})")