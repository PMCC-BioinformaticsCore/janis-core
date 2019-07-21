from core import ToolOutput, Stdout

from core.tool.commandtool import ToolInput
from core.types.common_data_types import Array, File
from core.unix.tools.unixtool import UnixTool


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
