



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
    FileTestTool,
    StringTestTool,
    ArrayFileTestTool,
    ArrayComponentsTestTool,
    ComponentsMandatoryTestTool,
    ComponentsOptionalTestTool,
    ComponentsMandatoryArrayTestTool,
    ComponentsOptionalArrayTestTool
)


class CallWFInputTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inFileOpt', File(optional=True))
        self.input('inFileArr', Array(File()))
        self.input('inFileArrOpt', Array(File(), optional=True))
        self.input('inStr', String())
        self.input('inStrOpt', String())
        self.input('inStrArr', Array(String()))
        self.input('inStrArrOpt', Array(String(), optional=True))
        self.input('inInt', Int())
        self.input('inIntOpt', Int(optional=True))
        self.input('inIntArr', Array(Int()))
        self.input('inIntArrOpt', Array(Int(), optional=True))
        self.input('inBool', Boolean())

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
        self.step(
            "stp2", 
            ComponentsOptionalTestTool(
                flag_true=self.inBool,
                flag_false=self.inBool,
                pos_optional=self.inStr,
                opt_optional=self.inStr,
            )
        )
        self.step(
            "stp3", 
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
        self.step(
            "stp4", 
            ComponentsOptionalArrayTestTool(
                pos_optional_arr=self.inFileArrOpt,
                opt_optional_arr=self.inStrArrOpt,
                opt_optional_arr_prefixeach=self.inStrArrOpt,
            )
        )
        self.output("outFile1", File, source=self.stp1.out)

    def friendly_name(self):
        return "TEST: CallWFInputTestWF"

    def id(self) -> str:
        return self.__class__.__name__



class CallDuplicateUseageWFInputTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inFileOpt', File(optional=True))
        self.input('inFileArr', Array(File()))
        self.input('inFileArrOpt', Array(File(), optional=True))
        self.input('inStr', String())
        self.input('inStrOpt', String())
        self.input('inStrArr', Array(String()))
        self.input('inStrArrOpt', Array(String(), optional=True))
        self.input('inInt', Int())
        self.input('inIntOpt', Int(optional=True))
        self.input('inIntArr', Array(Int()))
        self.input('inIntArrOpt', Array(Int(), optional=True))
        self.input('inBool', Boolean())

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
        self.step(
            "stp1_supp", 
            ComponentsMandatoryTestTool(
                pos_basic=self.inFile,
                pos_default=10,
                flag_true=True,
                flag_false=True,
                opt_basic="hi",
                opt_default=10,
            )
        )
        self.step(
            "stp2", 
            ComponentsOptionalTestTool(
                flag_true=self.inBool,
                flag_false=self.inBool,
                pos_optional=self.inStr,
                opt_optional=self.inStr,
            )
        )
        self.step(
            "stp2_supp", 
            ComponentsOptionalTestTool(
                flag_true=None,
                flag_false=None,
                pos_optional=None,
                opt_optional=None,
            )
        )
        self.step(
            "stp3", 
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
        self.step(
            "stp3_supp", 
            ComponentsMandatoryArrayTestTool(
                pos_basic_arr=self.inFileArr,
                pos_default_arr=[1, 2, 3],
                # flag_true=self.inBool,
                # flag_false=self.inBool,
                opt_basic_arr=['hi', 'there'],
                opt_default_arr=[1, 2, 3],
                opt_basic_arr_prefixeach=['hi', 'there'],
                opt_default_arr_prefixeach=['hi', 'there'],
            )
        )
        self.step(
            "stp4", 
            ComponentsOptionalArrayTestTool(
                pos_optional_arr=self.inFileArrOpt,
                opt_optional_arr=self.inStrArrOpt,
                opt_optional_arr_prefixeach=self.inStrArrOpt,
            )
        )
        self.step(
            "stp4_supp", 
            ComponentsOptionalArrayTestTool(
                pos_optional_arr=None,
                opt_optional_arr=None,
                opt_optional_arr_prefixeach=None,
            )
        )
        self.output("outFile1", File, source=self.stp1.out)

    def friendly_name(self):
        return "TEST: CallDuplicateUseageWFInputTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class CallStaticTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File())
        self.input('inFileOpt', File(optional=True))
        self.input('inFileArr', Array(File()))
        self.input('inFileArrOpt', Array(File(), optional=True))
        self.input('inStr', String())
        self.input('inStrOpt', String())
        self.input('inStrArr', Array(String()))
        self.input('inStrArrOpt', Array(String(), optional=True))
        self.input('inInt', Int())
        self.input('inIntOpt', Int(optional=True))
        self.input('inIntArr', Array(Int()))
        self.input('inIntArrOpt', Array(Int(), optional=True))
        self.input('inBool', Boolean())

        # full inputs static
        self.step(
            "stp1", 
            ComponentsMandatoryTestTool(
                pos_basic=self.inFile,
                pos_default=100,
                flag_true=False,
                flag_false=True,
                opt_basic="static",
                opt_default=100,
            )
        )
        self.step(
            "stp2", 
            ComponentsOptionalTestTool(
                flag_true=None,
                flag_false=None,
                pos_optional=None,
                opt_optional=None,
            )
        )
        self.step(
            "stp3", 
            ComponentsMandatoryArrayTestTool(
                pos_basic_arr=self.inFileArr,
                pos_default_arr=[1, 2, 3],
                # flag_true=self.inBool,
                # flag_false=self.inBool,
                opt_basic_arr=['hi', 'there'],
                opt_default_arr=[1, 2, 3],
                opt_basic_arr_prefixeach=['hi', 'there'],
                opt_default_arr_prefixeach=['hi', 'there'],
            )
        )
        self.step(
            "stp4", 
            ComponentsOptionalArrayTestTool(
                pos_optional_arr=None,
                opt_optional_arr=None,
                opt_optional_arr_prefixeach=None,
            )
        )

        self.output("outFile2", File, source=self.stp2.out)

    def friendly_name(self):
        return "TEST: CallStaticTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class CallPartialStaticTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)

        # partial inputs static
        self.step(
            "stp3", 
            ComponentsMandatoryTestTool(
                pos_basic=self.inFile,
                opt_basic="static",
                opt_default=100,
                opt_optional='',
            )
        )

        self.output("outFile3", File, source=self.stp3.out)

    def friendly_name(self):
        return "TEST: CallPartialStaticTestWF"

    def id(self) -> str:
        return self.__class__.__name__


class CallMinimalTestWF(Workflow):

    def constructor(self):
        self.input('inFile', File)
        self.input('inStr', String)
        self.input('inInt', Int)
        self.input('inBool', Boolean)

        # minimal inputs
        self.step(
            "stp4", 
            ComponentsMandatoryTestTool(
                pos_basic=self.inFile,
                opt_basic=self.inStr,
            )
        )
        self.output("outFile4", File, source=self.stp4.out)

    def friendly_name(self):
        return "TEST: CallMinimalTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# ConnectionsTestWF 
# File, String, Int workflow inputs & tool inputs 

class CallConnectionsTestWF(Workflow):

    def constructor(self):
        self.input('inStr', String)
        self.input('inFile', File)
        self.input('inInt', Int)

        self.step(
            "stp1", 
            StringTestTool(inp=self.inStr)
        )
        self.step(
            "stp2", 
            FileTestTool(inp=self.stp1.out)
        )

        self.output("outString", File, source=self.stp1.out)
        self.output("outFile", File, source=self.stp2.out)

    def friendly_name(self):
        return "TEST: CallConnectionsTestWF"

    def id(self) -> str:
        return self.__class__.__name__



# ArrayStepInputsTestWF
# calling same tool with different step input arrays
class CallArrayInputsTestWF(Workflow):

    def constructor(self):
        self.input('inFileArray', Array(File))
        self.input('inFileArrayOpt', Array(File, optional=True))
        self.input('inStrArray', Array(String))
        self.input('inIntArray', Array(Int))
        # self.input('inBoolArray', Array(Boolean))

        # full inputs
        self.step(
            "stp1",
            ArrayComponentsTestTool(
                pos_basic=self.inFileArray,
                pos_basic2=self.inFileArrayOpt,
                pos_default=self.inIntArray,
                pos_optional=self.inStrArray,
                # flag_true=self.inBoolArray,
                # flag_false=self.inBoolArray,
                opt_basic=self.inStrArray,
                opt_default=self.inIntArray,
                opt_optional=self.inStrArray,
            )
        )
        # full inputs static
        self.step(
            "stp2", 
            ArrayComponentsTestTool(
                pos_basic=self.inFileArray,
                pos_default=[4,5,6],
                pos_optional=["hi", "there", "friend"],
                # flag_true=[True],
                # flag_false=[True],
                opt_basic=["hi", "there", "friend"],
                opt_default=[4,5,6],
                opt_optional=["hi", "there", "friend"],
            )
        )
        # minimal inputs
        self.step(
            "stp3", 
            ArrayComponentsTestTool(
                pos_basic=self.inFileArray,
                opt_basic=self.inStrArray,
            )
        )

        self.output("outFile1", File, source=self.stp1.out)
        self.output("outFile2", File, source=self.stp2.out)
        self.output("outFile3", File, source=self.stp3.out)

    def friendly_name(self):
        return "TEST: CallArrayInputsTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# ArrayStepConnectionsTestWF 
# Array(File) step io connections
# janis does not allow non-file outputs, so can't test
# Array(String) -> Array(String) step connections. 

class CallArrayConnectionsTestWF(Workflow):

    def constructor(self):
        self.input('inStrArray', Array(String))
        self.input('inFileArray', Array(File))
        self.input('inIntArray', Array(Int))

        self.step(
            "stp1",
            ArrayFileTestTool(
                inp=self.inFileArray,
            ),
        )
        self.step(
            "stp2",
            ArrayFileTestTool(
                inp=self.stp1.out,
            ),
        )
        self.output("outFiles1", source=self.stp1.out)
        self.output("outFiles2", source=self.stp2.out)

    def friendly_name(self):
        return "TEST: CallArrayConnectionsTestWF"

    def id(self) -> str:
        return self.__class__.__name__