
import WDL 

from janis_core import CommandToolBuilder
from janis_core import ToolInput
from janis_core import ToolOutput

from .requirements import (
    parse_container_requirement, 
    parse_memory_requirement, 
    parse_disk_requirement, 
    parse_cpus_requirement
)

from .expressions import parse_expr
from .types import parse_type
from .command import parse_command


def parse_task(task: WDL.Tree.Task):
    # metadata
    name = task.name
    version = "DEV"
    
    # requirements
    container = parse_container_requirement(task)
    memory = parse_memory_requirement(task)
    disk = parse_disk_requirement(task)
    cpus = parse_cpus_requirement(task)
    
    # io
    inputs = [parse_command_tool_input(x) for x in task.inputs if not x.name.startswith("runtime_")]
    outputs = [parse_command_tool_output(o) for o in task.outputs]

    # generate first-pass CommandToolBuilder
    internal = CommandToolBuilder(
        tool=name,
        version=version,
        base_command=None,
        container=container,
        memory=memory,
        disk=disk,
        cpus=cpus,
        inputs=inputs,
        outputs=outputs,
        files_to_create=None
    )
    
    # augment inputs / outputs using command block
    internal = parse_command(internal, task)
    return internal

def parse_command_tool_input(inp: WDL.Tree.Decl) -> ToolInput:
    default = None
    if inp.expr:
        default = parse_expr(inp.expr)
    tinput = ToolInput(inp.name, parse_type(inp.type, uuid=inp.name), default=default)
    return tinput

def parse_command_tool_output(out: WDL.Tree.Decl) -> ToolOutput:
    if out.expr is None:
        raise Exception(f"Output {out.name} has no expression")
    sel = parse_expr(out.expr)
    tout = ToolOutput(out.name, parse_type(out.type, uuid=out.name), selector=sel)
    return tout