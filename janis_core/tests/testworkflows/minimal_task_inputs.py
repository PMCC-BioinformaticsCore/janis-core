

from typing import Optional

from janis_core import Workflow
from janis_core import CommandTool, ToolInput, ToolOutput
from janis_core.types import (
    String,
    File,
    Int,
    Stdout
)


class MinimalTaskInputsTestWF1(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr1', String())
        self.input('inStr2', String())
        self.input('inStr3', String())
        self.input('inInt1', Int())
        self.input('inInt2', Int())
        self.input('inInt3', Int())

        self.step(
            "stp1", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr1=self.inStr1, 
                inInt1=5, 
                inInt2=self.inInt1, 
            )
        )
        
    def friendly_name(self):
        return "TEST: MinimalTaskInputsTestWF1"

    def id(self) -> str:
        return self.__class__.__name__


class MinimalTaskInputsTestWF2(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr1', String())
        self.input('inStr2', String())
        self.input('inStr3', String())
        self.input('inInt1', Int())
        self.input('inInt2', Int())
        self.input('inInt3', Int())

        self.step(
            "stp1", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr1=self.inStr1,
                inStr2=self.inStr2,
                inStr3='hello',
                inStr4=self.inStr3,
                inInt1=5, 
                inInt2=self.inInt1,
                inInt3=100,
            )
        )
        self.step(
            "stp2", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr1=self.inStr1,
                inStr4='there',
                inInt1=5, 
                inInt2=self.inInt2, 
                inInt3=200, 
            )
        )
        
    def friendly_name(self):
        return "TEST: MinimalTaskInputsTestWF2"

    def id(self) -> str:
        return self.__class__.__name__



class MinimalTaskInputsTestWF3(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr1', String())
        self.input('inStr2', String())
        self.input('inStr3', String())
        self.input('inInt1', Int())
        self.input('inInt2', Int())
        self.input('inInt3', Int())

        self.step(
            "stp1", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr1=self.inStr1,
                inStr2=self.inStr2,
                inStr3='hello',
                inStr4=self.inStr3,
                inInt1=5, 
                inInt2=self.inInt1,
                inInt3=100,
            )
        )
        self.step(
            "stp2", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr1=self.inStr1,
                inStr4='there',
                inInt1=5, 
                inInt2=self.inInt2, 
                inInt3=200, 
            )
        )
        self.step(
            "stp3", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr4='friend',
                inInt4=300, 
            )
        )
        
    def friendly_name(self):
        return "TEST: MinimalTaskInputsTestWF3"

    def id(self) -> str:
        return self.__class__.__name__
    


class MinimalTaskInputsTestWF4(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr1', String())
        self.input('inStr2', String())
        self.input('inStr3', String())
        self.input('inInt1', Int())
        self.input('inInt2', Int())
        self.input('inInt3', Int())

        self.step(
            "stp1", 
            SubMinimalTaskInputsTestWF(
                inFile=self.inFile, 
                inStr1=self.inStr1, 
                inInt1=5, 
                inInt2=self.inInt1, 
            )
        )
        
    def friendly_name(self):
        return "TEST: MinimalTaskInputsTestWF4"

    def id(self) -> str:
        return self.__class__.__name__


class MinimalTaskInputsTestWF5(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr1', String())
        self.input('inStr2', String())
        self.input('inStr3', String())
        self.input('inInt1', Int())
        self.input('inInt2', Int())
        self.input('inInt3', Int())

        self.step(
            "stp1", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr1=self.inStr1, 
                inInt1=5, 
                inInt2=self.inInt1, 
            )
        )
        self.step(
            "stp2", 
            SubMinimalTaskInputsTestWF(
                inFile=self.inFile, 
                inStr1=self.inStr1, 
                inInt2=self.inInt1, 
            )
        )
        
    def friendly_name(self):
        return "TEST: MinimalTaskInputsTestWF5"

    def id(self) -> str:
        return self.__class__.__name__


class MinimalTaskInputsTestWF6(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr1', String())
        self.input('inStr2', String())
        self.input('inStr3', String())
        self.input('inInt1', Int())
        self.input('inInt2', Int())
        self.input('inInt3', Int())

        self.step(
            "stp1", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr2=self.inStr2,
                inStr3='hello',
                inStr4=self.inStr3,
                inInt1=5, 
                inInt2=self.inInt1,
                inInt3=100,
            )
        )
        self.step(
            "stp2", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr4='there',
                inInt1=5, 
                inInt2=self.inInt2, 
                inInt3=200, 
            )
        )
        self.step(
            "stp3", 
            SubMinimalTaskInputsTestWF(
                inFile=self.inFile, 
                inInt1=5, 
                inInt2=self.inInt1, 
            )
        )
        
    def friendly_name(self):
        return "TEST: MinimalTaskInputsTestWF5"

    def id(self) -> str:
        return self.__class__.__name__


class SubMinimalTaskInputsTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inStr1', String(optional=True))
        self.input('inStr2', String(optional=True))
        self.input('inStr3', String(optional=True))
        self.input('inInt1', Int(optional=True))
        self.input('inInt2', Int(optional=True))
        self.input('inInt3', Int(optional=True))

        self.step(
            "stp1", 
            TaskInputsTestTool1(
                inFile=self.inFile, 
                inStr2=self.inStr2, 
                inStr3='hi', 
                inInt1=5, 
                inInt2=self.inInt2, 
                inInt4=999, 
            )
        )
        
    def friendly_name(self):
        return "TEST: SubMinimalTaskInputsTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class TaskInputsTestTool1(CommandTool):
    def tool(self) -> str:
        return "TaskInputsTestTool1"
    
    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput('inFile', File(), position=1),
            ToolInput('inStr1', String(optional=True), position=2),
            ToolInput('inStr2', String(optional=True), position=2),
            ToolInput('inStr3', String(optional=True), position=2),
            ToolInput('inStr4', String(optional=True), position=2),
            ToolInput('inInt1', Int(optional=True), position=3),
            ToolInput('inInt2', Int(optional=True), position=3),
            ToolInput('inInt3', Int(optional=True), position=3),
            ToolInput('inInt4', Int(optional=True), position=3),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
    

