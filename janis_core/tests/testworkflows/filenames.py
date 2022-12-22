

from janis_core import Workflow

from janis_core.types import (
    String,
    File,
)

from ..testtools import FilenameTestTool
from ..testtools import FilenameGeneratedTool
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



class FilenameGeneratedTestWF(Workflow):

    def constructor(self):
        self.input('inStr', String)
        self.input('inStrOpt', String(optional=True))
        self.input('inFile', File)
        self.input('inFileOpt', File(optional=True))

        self.step(
            "stp1", 
            FilenameGeneratedTool(
                inp=self.inStr,
                inpOptional=self.inStrOpt,
                fileInp=self.inFile,
                fileInpOptional=self.inFileOpt,
            )
        )

        self.output("out", File, source=self.stp1.out)

    def friendly_name(self):
        return "TEST: FilenameGeneratedTestWF"

    def id(self) -> str:
        return self.__class__.__name__

    
