

from janis_core import Workflow

from janis_core.types import (
    String,
    File,
)

from ..testtools import FilenameTestTool
from ..testtools import FilenameInputSelectorTestTool


class FilenameTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)

        self.step(
            "stp1", 
            FilenameTestTool(
                inp1=self.inFile,
                inp2=self.inStr,
            )
        )
        self.step(
            "stp2", 
            FilenameInputSelectorTestTool(
                inp1=self.inFile,
            )
        )

        self.output("outFilename", File, source=self.stp1.out)
        self.output("outFilenameInputSelector", File, source=self.stp2.out)

    def friendly_name(self):
        return "TEST: BasicIOTestWF"

    def id(self) -> str:
        return self.__class__.__name__
