from janis_core import ToolOutput, Stdout

from janis_core.tool.commandtool import ToolInput
from janis_core.types.common_data_types import Array, File
from janis_core.unix.tools.unixtool import UnixTool


class Cat(UnixTool):
    @staticmethod
    def tool():
        return "cat"

    def friendly_name(self):
        return "Concatenate"

    @staticmethod
    def base_command():
        return "cat"

    def inputs(self):
        return [ToolInput("files", Array(File()))]

    def outputs(self):
        return [ToolOutput("out", Stdout())]

    @staticmethod
    def docker():
        return "ubuntu:latest"
