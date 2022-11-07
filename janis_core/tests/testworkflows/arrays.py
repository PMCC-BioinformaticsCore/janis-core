



class ArraysTestWF(Workflow):
    def constructor(self):
        self.input('inFileArray', Array(File))

        self.step(
            "stp1",
            CoreTypesTestTool(
                inFile=self.inFileArray,
            ),
        )
        self.step(
            "stp2",
            CoreTypesTestTool(
                inFile=self.stp1.outFile,
            ),
        )

        self.output("outFile", source=self.stp2.outFile)

    def friendly_name(self):
        return "TEST: CoreTypesTestWF"

    def id(self) -> str:
        return self.__class__.__name__
