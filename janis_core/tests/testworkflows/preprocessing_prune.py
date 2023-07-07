

from typing import Optional

from janis_core import (
    Workflow, 
    CommandTool, 
    ToolInput, 
    ToolArgument, 
    ToolOutput, 
    InputSelector,
)

from janis_core.types import (
    String,
    Filename,
    Stdout,
    Array,
    File,
    Int,
)


### WORKFLOWS ###

class PruneFlatTW(Workflow):
    
    def tool(self):
        return "PruneTW"
    
    def friendly_name(self):
        return "TEST: PruneTW"

    def id(self) -> str:
        return self.__class__.__name__

    def constructor(self):
        self.input('inFile1', File())
        self.input('inFile2', File())
        self.input('inFile3', File())

        self.input('inStr1', String)
        self.input('inStr2', String)
        self.input('inStr3', String)

        self.input('inFileOpt1', File(optional=True))
        self.input('inFileOpt2', File(optional=True))
        self.input('inFileOpt3', File(optional=True))
        
        self.input('inStrOpt1', String(optional=True))
        self.input('inStrOpt2', String(optional=True))
        self.input('inStrOpt3', String(optional=True))

        self.step(
            "stp1", 
            PruneMandatoryTT(
                inFile1=self.inFileOpt1,
                inFile2=self.inFileOpt2,
                inStr1=self.inStrOpt1,
                inStr2=self.inStrOpt2,
            )
        )
        self.step(
            "stp2", 
            PruneInputRefTT(
                inFile1=self.inFileOpt1,
            )
        )
        self.step(
            "stp3", 
            PruneOutputRefTT()
        )
        self.step(
            "stp4", 
            PruneOptionalTT(
                inFileOpt1=self.stp3.outFile1,
                inStrOpt1=self.stp3.outStr1,
            )
        )
        self.step(
            "stp5", 
            PruneOptionalTT(
                inFileOpt2=self.inFileOpt2,
                inStrOpt2=self.inStrOpt2,
            )
        )
        self.step(
            "stp6", 
            PruneOptional2TT(
                inStrOpt1="hello",
                inStrOpt2="there",
                inStrOpt3="friend",
            )
        )
        
        self.output("outFile1", File, source=self.stp1.outFile1)
        self.output("outFile2", File, source=self.stp1.outFile2)


class PruneNestedTW(Workflow):
    
    def tool(self):
        return "PruneSubTW"
    
    def friendly_name(self):
        return "TEST: PruneSubTW"

    def id(self) -> str:
        return self.__class__.__name__

    def constructor(self):
        self.input('inFile1', File())
        self.input('inFile2', File())
        self.input('inFile3', File())

        self.input('inStr1', String)
        self.input('inStr2', String)
        self.input('inStr3', String)

        self.input('inFileOpt1', File(optional=True))
        self.input('inFileOpt2', File(optional=True))
        self.input('inFileOpt3', File(optional=True))
        
        self.input('inStrOpt1', String(optional=True))
        self.input('inStrOpt2', String(optional=True))
        self.input('inStrOpt3', String(optional=True))

        self.step(
            "stp1", 
            PruneOptionalTT(
                inFileOpt1=self.inFileOpt1,
                inStrOpt1=self.inStrOpt1,
            )
        )
        
        self.step(
            "stp2", 
            Sub1TW(
                inFileOpt2=self.inFileOpt2,
                inFileOpt3=self.inFileOpt3,
                inStrOpt2=self.inStrOpt2,
                inStrOpt3=self.inStrOpt3,
            )
        )
        
        self.output("outFile1", File, source=self.stp1.outFile1)
        self.output("outFile2", File, source=self.stp1.outFile2)


class Sub1TW(Workflow):
    
    def tool(self):
        return "Sub1TW"
    
    def friendly_name(self):
        return "TEST: Sub1TW"

    def id(self) -> str:
        return self.__class__.__name__

    def constructor(self):
        self.input('inFileOpt2', File(optional=True))
        self.input('inFileOpt3', File(optional=True))
        
        self.input('inStrOpt2', String(optional=True))
        self.input('inStrOpt3', String(optional=True))

        self.step(
            "stp1", 
            PruneOptionalTT(
                inFileOpt2=self.inFileOpt2,
                inStrOpt2=self.inStrOpt2,
            )
        )
        
        self.step(
            "stp2", 
            Sub2TW(
                inFileOpt3=self.inFileOpt3,
                inStrOpt3=self.inStrOpt3,
            )
        )
        
        self.output("outFile1", File, source=self.stp1.outFile1)
        self.output("outFile2", File, source=self.stp1.outFile2)
        self.output("outFile3", File, source=self.stp1.outFile3)


class Sub2TW(Workflow):
    
    def tool(self):
        return "Sub2TW"
    
    def friendly_name(self):
        return "TEST: Sub2TW"

    def id(self) -> str:
        return self.__class__.__name__

    def constructor(self):
        self.input('inFileOpt3', File(optional=True))
        self.input('inStrOpt3', String(optional=True))

        self.step(
            "stp1", 
            PruneOptionalTT(
                inFileOpt3=self.inFileOpt3,
                inStrOpt3=self.inStrOpt3,
            )
        )
        
        self.output("outFile1", File, source=self.stp1.outFile1)
        self.output("outFile2", File, source=self.stp1.outFile2)
        self.output("outFile3", File, source=self.stp1.outFile3)



### TOOLS ###

class PruneTTBase(CommandTool):
    
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"
    
    def id(self) -> str:
        return self.__class__.__name__
    
    def tool(self) -> str:
        return self.__class__.__name__
    
    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
    
class PruneMandatoryTT(PruneTTBase):
    
    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inFile1", File(), position=1),
            ToolInput("inFile2", File(), position=2),
            ToolInput("inStr1", String(), position=3),
            ToolInput("inStr2", String(), position=4),
        ]

    def outputs(self):
        return [
            ToolOutput('outFile1', File(), glob='out1*'),
            ToolOutput('outFile2', File(), glob='out2*'),
        ]

class PruneOptionalTT(PruneTTBase):
    
    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inFileOpt1", File(optional=True), position=1),
            ToolInput("inFileOpt2", File(optional=True), position=2),
            ToolInput("inFileOpt3", File(optional=True), position=3),
            ToolInput("inStrOpt1", String(optional=True), position=4),
            ToolInput("inStrOpt2", String(optional=True), position=5),
            ToolInput("inStrOpt3", String(optional=True), position=6),
        ]

    def outputs(self):
        return [
            ToolOutput('outFile1', File(), glob='out1*'),
            ToolOutput('outFile2', File(), glob='out2*'),
            ToolOutput('outFile3', File(), glob='out3*'),
        ]

class PruneOptional2TT(PruneTTBase):
    
    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inFileOpt1", File(optional=True), position=1),
            ToolInput("inFileOpt2", File(optional=True), position=2),
            ToolInput("inFileOpt3", File(optional=True), position=3),
            ToolInput("inStrOpt1", String(optional=True), position=4),
            ToolInput("inStrOpt2", String(optional=True), position=5),
            ToolInput("inStrOpt3", String(optional=True), position=6),
        ]

    def outputs(self):
        return [
            ToolOutput('outFile1', File(), glob='out1*'),
            ToolOutput('outFile2', File(), glob='out2*'),
            ToolOutput('outFile3', File(), glob='out3*'),
        ]

class PruneInputRefTT(PruneTTBase):
    
    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inFile1", File(), position=1),
            ToolInput(
                "OutputName", 
                Filename(
                    prefix=InputSelector('inFile1', remove_file_extension=True),
                    extension='.csv'
                ), 
                position=2
            ),
        ]

    def outputs(self):
        return [ToolOutput('outFile1', File(), glob='out1*')]

class PruneOutputRefTT(PruneTTBase):

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inFileOpt1", File(optional=True), position=1),
            ToolInput("inStrOpt1", String(optional=True), position=2),
        ]

    def outputs(self):
        return [
            ToolOutput("outFile1", File(), selector=InputSelector('inFileOpt1')),
            ToolOutput("outStr1", String(), selector=InputSelector('inStrOpt1')),
        ]


    

