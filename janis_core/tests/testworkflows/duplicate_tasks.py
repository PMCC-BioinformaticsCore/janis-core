



from janis_core import Workflow

from janis_core.types import (
    File,
    String,
    Int,
    Stdout, 
)

from typing import Optional
from janis_core import CommandTool, ToolInput, ToolOutput


# TODO add filename test

class DuplicateTasksTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr', String())
        self.input('inInt', Int())

        self.step(
            "stp1", 
            EchoTestTool(
                inFile=self.inFile,
                inStr='hello',
                inInt=5,
            )
        )
        self.step(
            "stp2", 
            EchoTestTool(
                inFile=self.inFile,
                inStr='there',
                inInt=10,
            )
        )
        self.step(
            "stp3", 
            EchoTestWorkflow1(
                inFile=self.inFile,
                inStr='my',
                inInt=15,
            )
        )
        self.step(
            "stp4", 
            EchoTestWorkflow1(
                inFile=self.inFile,
                inStr='friend',
                inInt=20,
            )
        )
        self.step(
            "stp5", 
            EchoTestWorkflow2(
                inFile=self.inFile,
                inStr='ok',
                inInt=99,
            )
        )

    def friendly_name(self):
        return "TEST: DuplicateTasksTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class EchoTestTool(CommandTool):
    
    def tool(self) -> str:
        return "EchoTestTool"
        
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('inFile', File(), position=1),
            ToolInput('inStr', String(), position=2),
            ToolInput('inInt', Int(), position=3),
        ]

    def outputs(self):
        return [ToolOutput('out', Stdout())]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
    


class EchoTestWorkflow1(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr', String())
        self.input('inInt', Int())

        self.step(
            "stp1", 
            EchoTestTool(
                inFile=self.inFile,
                inStr='good',
                inInt=25,
            )
        )
        self.step(
            "stp2", 
            EchoTestTool(
                inFile=self.inFile,
                inStr='night',
                inInt=30,
            )
        )

        self.output('out', File(), source=self.stp2.out)

    def friendly_name(self):
        return "TEST: EchoTestWorkflow"

    def id(self) -> str:
        return self.__class__.__name__


class EchoTestWorkflow2(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr', String())
        self.input('inInt', Int())

        self.step(
            "lol", 
            EchoTestWorkflow1(
                inFile=self.inFile,
                inStr='good',
                inInt=25,
            )
        )

    def friendly_name(self):
        return "TEST: EchoTestWorkflow"

    def id(self) -> str:
        return self.__class__.__name__


# class EchoTestWorkflow(Workflow):

#     def constructor(self):
#         self.input('inFile', File())
#         self.input('inStr', String())
#         self.input('inInt', Int())

#         self.step(
#             "stp1", 
#             EchoTestTool(
#                 inFile=self.inFile,
#                 inStr=self.inStr,
#                 inInt=self.inInt,
#             )
#         )
#         self.step(
#             "stp2", 
#             EchoTestTool(
#                 inFile=self.inFile,
#                 inStr=self.inStr,
#                 inInt=self.inInt,
#             )
#         )

#     def friendly_name(self):
#         return "TEST: EchoTestWorkflow"

#     def id(self) -> str:
#         return self.__class__.__name__
    

    
