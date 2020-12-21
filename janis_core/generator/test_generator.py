import janis_core as jc
from janis_core import *


CLT = jc.CommandToolBuilder(
    tool="devtool",
    base_command="echo",
    inputs=[jc.ToolInput("inp", str, position=1)],
    outputs=[jc.ToolOutput("out", jc.Stdout)],
    version="v1.0",
    container="ubuntu",
)

if __name__ == "__main__":
    from janis_core.generator.util import convert_generic_class

    print(convert_generic_class(CLT))
