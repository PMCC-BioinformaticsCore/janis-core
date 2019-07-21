#
# Untar a file
from core import Array, ToolInput, ToolOutput, InputSelector, File, Filename
from core.unix.data_types.tarfile import TarFile
from core.unix.tools.unixtool import UnixTool


class Tar(UnixTool):
    @staticmethod
    def tool():
        return "Tar"

    def friendly_name(self):
        return "Tar archive"

    @staticmethod
    def base_command():
        return ["tar", "cvf"]

    def inputs(self):
        return [
            ToolInput("files", Array(File()), position=2),
            ToolInput("files2", Array(File(), optional=True), position=3),
            ToolInput("outputFilename", Filename(extension=".tar"), position=1),
        ]

    def outputs(self):
        return [ToolOutput("out", TarFile(), glob=InputSelector("outputFilename"))]


if __name__ == "__main__":
    print(Tar().help())
