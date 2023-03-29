
"""
tricky workflow to test data source categories & data source variables.
this is only interested in checking whether ToolInputs are declared as process_inputs, 
param_inputs, or internal_inputs. 

data source categories:
    for each CommandTool / PythonTool in entity, decides which ToolInputs are:
        - process inputs
        - param inputs
        - internal inputs

data source variables:
    for each CommandTool / PythonTool, decides the variable_name (a process_input, a param_input, or None ) 
    which feeds data for the tinput.
    the variable_name is how the particular tinput_id will be referenced inside its process.    

includes nested subworkflows to check these situations. 

"""

from typing import Optional

from janis_core import (
    Workflow,
    CommandTool,
    ToolInput,
    ToolOutput,
    InputSelector
)
from janis_core.types import (
    String,
    File,
)


class DataSourceTestWF(Workflow):
    def constructor(self):
        self.input('inFile1', File)
        self.input('inFileOpt1', File(optional=True))
        self.input('inStr1', String)
        self.input('inStrOpt1', String(optional=True))

        self.step(
            "stp1", 
            DataSourceTestTool(
                inFile=self.inFile1,
                inFileOpt=self.inFileOpt1,
                inStr=self.inStr1,
                inStrOpt=self.inStrOpt1,
            )
        )
        self.step(
            "stp2", 
            DataSourceTestTool(
                inFile=self.inFile1,
                inFileOpt=self.inFileOpt1,
                inStr='hello',
                inStrOpt='there',
            )
        )
        self.step(
            "stp3", 
            DataSourceTestTool(
                inFile=self.inFile1,
                inStr='hello',
            )
        )
        # fully connected
        self.step(
            "stp4", 
            SubDataSourceTestWF(
                inFile2=self.inFile1,
                inFileOpt2=self.inFileOpt1,
                inStr2=self.inStr1,
                inStrOpt2=self.inStrOpt1,
            )
        )
        # statics
        self.step(
            "stp5", 
            SubDataSourceTestWF(
                inFile2=self.inFile1,
                inFileOpt2=self.inFileOpt1,
                inStr2='hello',
                inStrOpt2='there',
            )
        )
        # omitted + statics
        self.step(
            "stp6", 
            SubDataSourceTestWF(
                inFile2=self.inFile1,
                inStr2='hello',
            )
        )

    def friendly_name(self):
        return "TEST: DataSourceTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class SubDataSourceTestWF(Workflow):
    def constructor(self):
        self.input('inFile2', File)
        self.input('inFileOpt2', File(optional=True))
        self.input('inStr2', String)
        self.input('inStrOpt2', String(optional=True))

        self.step(
            "stp1", 
            DataSourceTestTool(
                inFile=self.inFile2,
                inFileOpt=self.inFileOpt2,
                inStr=self.inStr2,
                inStrOpt=self.inStrOpt2,
            )
        )
        self.step(
            "stp2", 
            DataSourceTestTool(
                inFile=self.inFile2,
                inFileOpt=self.inFileOpt2,
                inStr='hello',
                inStrOpt='there',
            )
        )
        self.step(
            "stp3", 
            DataSourceTestTool(
                inFile=self.inFile2,
                inStr='hello',
            )
        )
        self.step(
            "stp4", 
            SubSubDataSourceTestWF(
                inFile3=self.inFile2,
                inFileOpt3=self.inFileOpt2,
                inStr3=self.inStr2,
                inStrOpt3=self.inStrOpt2,
            )
        )
        self.step(
            "stp5", 
            SubSubDataSourceTestWF(
                inFile3=self.inFile2,
                inFileOpt3=self.inFileOpt2,
                inStr3='hello',
                inStrOpt3='there',
            )
        )
        self.step(
            "stp6", 
            SubSubDataSourceTestWF(
                inFile3=self.inFile2,
                inStr3='hello',
            )
        )

    def friendly_name(self):
        return "TEST: SubDataSourceTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class SubSubDataSourceTestWF(Workflow):
    def constructor(self):
        self.input('inFile3', File)
        self.input('inFileOpt3', File(optional=True))
        self.input('inStr3', String)
        self.input('inStrOpt3', String(optional=True))

        self.step(
            "stp1", 
            DataSourceTestTool(
                inFile=self.inFile3,
                inFileOpt=self.inFileOpt3,
                inStr=self.inStr3,
                inStrOpt=self.inStrOpt3,
            )
        )
        self.step(
            "stp2", 
            DataSourceTestTool(
                inFile=self.inFile3,
                inFileOpt=self.inFileOpt3,
                inStr='hello',
                inStrOpt='there',
            )
        )
        self.step(
            "stp3", 
            DataSourceTestTool(
                inFile=self.inFile3,
                inStr='hello',
            )
        )

    def friendly_name(self):
        return "TEST: SubSubDataSourceTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class DataSourceTestTool(CommandTool):
    
    def tool(self) -> str:
        return "DataSourceTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("inFile", File(), position=1),
            ToolInput("inFileOpt", File(optional=True), position=2),
            ToolInput("inStr", String(), position=3),
            ToolInput("inStrOpt", String(optional=True), position=4),
        ]

    def outputs(self):
        return [
            ToolOutput("outFile", File(), selector=InputSelector('inFile')),
            ToolOutput("outFileOpt", File(optional=True), selector=InputSelector('inFileOpt')),
            ToolOutput("outStr", String(), selector=InputSelector('inStr')),
            ToolOutput("outStrOpt", String(optional=True), selector=InputSelector('inStrOpt')),
        ]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"