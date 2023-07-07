
from janis_core import String, ToolInput, ToolOutput, Stdout, Boolean, ToolMetadata, CommandTool


class Echo(CommandTool):

    def tool_module(self):
        return "unix"
    
    def tool(self):
        return "echo"

    def friendly_name(self):
        return "Echo"

    def container(self):
        return "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"

    def version(self):
        return "v1.0.0"

    def base_command(self):
        return "echo"

    def inputs(self):
        return [
            ToolInput("inp", String(), position=1),
            ToolInput(
                "include_newline",
                Boolean(optional=True),
                prefix="-n",
                doc="Do not print the trailing newline character.  This may also be achieved by appending `\c' to "
                "the end of the string, as is done by iBCS2 compatible systems.  Note that this option as well as the "
                "effect of `\c' are implementation-defined in IEEE Std 1003.1-2001 (``POSIX.1'') as amended by "
                "Cor. 1-2002.  Applications aiming for maximum portability are strongly encouraged to use printf(1) "
                "to suppress the newline character.",
            ),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    def bind_metadata(self):
        self.metadata.documentation = """\
The echo utility writes any specified operands, separated by single blank (` ') characters \
and followed by a newline (`\n') character, to the standard output.

Some shells may provide a builtin echo command which is similar or identical to this utility. \
Most notably, the builtin echo in sh(1) does not accept the -n option. Consult the builtin(1) manual page."""

    def time(self, hints):
        return 60  # seconds
    