



from janis_core import (
    Workflow,
)

from janis_core.types import (
    Array,
)
from janis_bioinformatics.data_types.bam import BamBai
from janis_core.tests.testworkflows import (
    SecondariesTestTool,
)


# ------------------- #
#  DISGUSTING COMBOS  #
# ------------------- #


class ScatterSecondariesTestWF(Workflow):
    def id(self) -> str:
        return "ScatterSecondaries"

    def friendly_name(self):
        return "WF which uses Scatter and Secondaries"

    def constructor(self):
        self.input('inAlignments', Array(BamBai))
        
        self.step(
            "stp1", 
            SecondariesTestTool(
                inp=self.inAlignments
            ),
            scatter="inp"
        )

        self.output("outBamBaiArray", Array(BamBai), source=self.stp1.out)
        # self.output("outStdout", source=self.stp1.outStdout)


# class ArrayScatterTestWF(Workflow):
#     def id(self) -> str:
#         return "ArrayScatterTestWF"

#     def friendly_name(self):
#         return "WF which uses Array(File) and Scatter"

#     def constructor(self):
#         self.input('inStrArray', Array(String))
#         self.input('inFileArray', Array(File))
#         self.input('inIntArray', Array(Int))

#         self.step(
#             "stp1", 
#             (testtool=self.inStrArray), 
#             scatter="testtool"
#         )

#         self.output("outStrArray", source=self.stp1.out)

