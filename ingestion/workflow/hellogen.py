from datetime import datetime
from typing import List, Optional, Dict, Any

from janis_core import *
from janis_core.types.common_data_types import File, Array, Boolean, Stdout

Cat_V1_0_0 = CommandToolBuilder(
    tool="cat",
    base_command="cat",
    inputs=[
        ToolInput(
            tag="file", input_type=File(optional=True), doc=InputDocumentation(doc=None)
        ),
        ToolInput(
            tag="files",
            input_type=Array(t=File(), optional=True),
            position=1,
            doc=InputDocumentation(doc=None),
        ),
        ToolInput(
            tag="number_output",
            input_type=Boolean(optional=True),
            prefix="-n",
            doc=InputDocumentation(doc="Number the output lines, starting at 1."),
        ),
        ToolInput(
            tag="number_non_blank",
            input_type=Boolean(optional=True),
            prefix="-b",
            doc=InputDocumentation(
                doc="Number the non-blank output lines, starting at 1."
            ),
        ),
        ToolInput(
            tag="disable_output_buffer",
            input_type=Boolean(optional=True),
            prefix="-u",
            doc=InputDocumentation(doc="Disable output buffering."),
        ),
        ToolInput(
            tag="squeeze",
            input_type=Boolean(optional=True),
            prefix="-s",
            doc=InputDocumentation(
                doc="Squeeze multiple adjacent empty lines, causing the output to be single spaced."
            ),
        ),
        ToolInput(
            tag="display_nonprint_and_eol_chars",
            input_type=Boolean(optional=True),
            prefix="-e",
            doc=InputDocumentation(
                doc="Display non-printing characters (see the -v option), and display a dollar sign (`$') at the end of each line."
            ),
        ),
        ToolInput(
            tag="display_nonprint_and_tab_chars",
            input_type=Boolean(optional=True),
            prefix="-t",
            doc=InputDocumentation(
                doc="Display non-printing characters (see the -v option), and display tab characters as `^I'."
            ),
        ),
        ToolInput(
            tag="display_nonprint_chars",
            input_type=Boolean(optional=True),
            prefix="-v",
            doc=InputDocumentation(
                doc="Display non-printing characters so they are visible.  Control characters print as `^X' for control-X; the delete character (octal 0177) prints as `^?'.  Non-ASCII characters (with the high bit set) are printed as `M-' (for meta) followed by the character for the low 7 bits."
            ),
        ),
    ],
    outputs=[
        ToolOutput(
            tag="out",
            output_type=Stdout(subtype=File(), optional=False),
            doc=OutputDocumentation(doc=None),
        )
    ],
    container="ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956",
    version="v1.0.0",
    friendly_name="Concatenate",
    tool_module="unix",
    metadata=ToolMetadata(
        contributors=[],
        dateUpdated=datetime(2019, 7, 26),
        documentation="The cat utility reads files sequentially, writing them to the standard output. The file operands are processed in command-line order. If file is a single dash (`-') or absent,cat reads from the standard input. If file is a UNIX domain socket, cat connects to it and then reads it until EOF. This complements the UNIX domain binding capability available in inetd(8).",
    ),
)

Cat_V1_0_0().translate("wdl")
