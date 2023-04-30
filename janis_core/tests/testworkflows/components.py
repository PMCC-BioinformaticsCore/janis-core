



from janis_core import (
    Workflow
)

from janis_core.types import (
    String,
    File,
    Array,
    Int,
    Boolean,
)

from janis_core.tests.testtools import (
    ComponentsMandatoryTestTool,
    ComponentsOptionalTestTool,
    ComponentsMandatoryArrayTestTool,
    ComponentsOptionalArrayTestTool
)


# StepInputsTestWF
# calling same tool with different step inputs


### WORKFLOWS ### -------------------------------------------------------

class ComponentsMandatoryTestWF(Workflow):

    def friendly_name(self):
        return "TEST: ComponentsMandatoryTestWF"

    def id(self) -> str:
        return self.__class__.__name__

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)
 
        # full inputs
        self.step(
            "stp1", 
            ComponentsMandatoryTestTool(
                pos_basic=self.inFile,
                pos_default=self.inInt,
                flag_true=self.inBool,
                flag_false=self.inBool,
                opt_basic=self.inStr,
                opt_default=self.inInt,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ComponentsMandatoryTestTool(
                pos_basic=self.inFile,
                pos_default=100,
                flag_true=False,
                flag_false=True,
                opt_basic="static",
                opt_default=100,
            )
        )


class ComponentsOptionalTestWF(Workflow):

    def friendly_name(self):
        return "TEST: ComponentsOptionalTestWF"

    def id(self) -> str:
        return self.__class__.__name__

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)
 
        # full inputs
        self.step(
            "stp1", 
            ComponentsOptionalTestTool(
                flag_true=self.inBool,
                flag_false=self.inBool,
                pos_optional=self.inStr,
                opt_optional=self.inStr,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ComponentsOptionalTestTool(
                flag_true=False,
                flag_false=True,
                pos_optional="static",
                opt_optional='static',
            )
        )



class ComponentsMandatoryArrayTestWF(Workflow):

    def friendly_name(self):
        return "TEST: ComponentsMandatoryArrayTestWF"

    def id(self) -> str:
        return self.__class__.__name__

    def constructor(self):
        self.input('inFileArr', Array(File))
        self.input('inStrArr', Array(String))
        self.input('inIntArr', Array(Int))
        self.input('inBoolArr', Array(Boolean))
 
        # full inputs
        self.step(
            "stp1", 
            ComponentsMandatoryArrayTestTool(
                pos_basic_arr=self.inFileArr,
                pos_default_arr=self.inIntArr,
                # flag_true=self.inBool,
                # flag_false=self.inBool,
                opt_basic_arr=self.inStrArr,
                opt_default_arr=self.inIntArr,
                opt_basic_arr_prefixeach=self.inStrArr,
                opt_default_arr_prefixeach=self.inStrArr,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ComponentsMandatoryArrayTestTool(
                pos_basic_arr=self.inFileArr,
                pos_default_arr=[100, 200, 300],
                # flag_true=False,
                # flag_false=True,
                opt_basic_arr=["hi", "there"],
                opt_default_arr=[100, 200, 300],
                opt_basic_arr_prefixeach=["hi", "there"],
                opt_default_arr_prefixeach=["hi", "there"],
            )
        )


class ComponentsOptionalArrayTestWF(Workflow):

    def friendly_name(self):
        return "TEST: ComponentsOptionalArrayTestWF"

    def id(self) -> str:
        return self.__class__.__name__

    def constructor(self):
        self.input('inFileArr', Array(File))
        self.input('inStrArr', Array(String))
        self.input('inIntArr', Array(Int))
        self.input('inBoolArr', Array(Boolean))
 
        # full inputs
        self.step(
            "stp1", 
            ComponentsOptionalArrayTestTool(
                pos_optional_arr=self.inFileArr,
                opt_optional_arr=self.inStrArr,
                opt_optional_arr_prefixeach=self.inStrArr,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ComponentsOptionalArrayTestTool(
                pos_optional_arr=self.inFileArr,
                opt_optional_arr=["hi", "there"],
                opt_optional_arr_prefixeach=["hi", "there"],
            )
        )


