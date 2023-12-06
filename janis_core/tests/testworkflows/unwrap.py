

from typing import Optional

from janis_core import (
    CommandTool,
    ToolInput,
    ToolArgument,
    ToolOutput,
    Workflow
)
from janis_core.types import (
    Filename,
    File,
    String,
    Int,
    Float,
    Array,
)
from janis_core.redefinitions.types import Bam, BamBai
import janis_core as j


# WORKFLOW
class UnwrapTestWF(Workflow):

    def constructor(self):
        self.input("inFile", File())
        self.input("inFileOpt", File(optional=True))
        self.input("inFileArr", Array(File()))
        self.input("inBamBai", BamBai())
        self.input("inBamBaiArr", Array(BamBai()))
        self.input("inStr", String())
        self.input("inStrOpt", String(optional=True))
        self.input("inInt", Int())
        self.input("inIntOpt", Int(optional=True))
        self.input("inFloat", Float())
        self.input("inFloatOpt", Float(optional=True))

        self.step(
            "basics_step", 
            BasicsTestTool(
                reads_1=self.inFile,
                bam_sorted=self.inFile,
            )
        )
        
        self.step(
            "selectors_step", 
            SelectorTestTool(
                inFile=self.inFile,
            )
        )
        
        self.step(
            "standard_step", 
            StandardTestTool(
                inFile=self.selectors_step.out,
                inFileArr=self.inFileArr,
                inBamBai=self.inBamBai,
                inBamBaiArr=self.inBamBaiArr,
                inStr=self.inStr,
                inStrOpt=self.inStrOpt,
                inInt=self.inInt,
                inIntOpt=self.inIntOpt,
            )
        )
        
        self.step(
            "logical_step", 
            LogicalTestTool(
                inFile=self.inFile,
                inFileOpt=self.inFileOpt,
                inFileArr=self.inFileArr,
                inBamBai=self.inBamBai,
                inBamBaiArr=self.inBamBaiArr,
                inStr=self.inStr,
                inStrOpt=self.inStrOpt,
                inInt=self.inInt,
                inIntOpt=self.inIntOpt,
                inFloat=self.inFloat,
                inFloatOpt=self.inFloatOpt,
            )
        )
        
        self.step(
            "complex_step", 
            ComplexOperatorTestTool(
                inFile=self.inFile,
                inFileOpt=self.inFileOpt,
                inFileArr=self.inFileArr,
                inBamBai=self.inBamBai,
                inBamBaiArr=self.inBamBaiArr,
                inStr=self.inStr,
                inStrOpt=self.inStrOpt,
                inInt=self.inInt,
                inIntOpt=self.inIntOpt,
                inFloat=self.inFloat,
                inFloatOpt=self.inFloatOpt,
            )
        )

    def friendly_name(self):
        return "TEST: UnwrapTestWF"

    def id(self) -> str:
        return self.__class__.__name__


# TOOLS
class BasicsTestTool(CommandTool):
    
    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
        
    def base_command(self) -> Optional[str | list[str]]:
        return ['samtools', 'index']

    def friendly_name(self):
        return "TEST: BasicsTestTool"

    def tool(self):
        return "BasicsTestTool"
    
    def inputs(self):
        return [
            ToolInput("reads_1", File()),
            ToolInput("in_filename", Filename(), position=3),
            ToolInput("bam_sorted", Bam(), position=2),
        ]
    
    def arguments(self):
        return [
            ToolArgument('-b', position=1)
        ]

    def outputs(self):
        return [
            ToolOutput(
                "bam_sorted_indexed",
                BamBai(),
                selector=j.BasenameOperator(j.InputSelector('bam_sorted')),
            ),
            ToolOutput(
                "trimmed_reads_1",
                File(),
                selector=j.StringFormatter(
                    format='{token1}.trimmed{token2}', 
                    token1=j.BasenameOperator(j.InputSelector("reads_1")), 
                    token2=j.NameextOperator(j.InputSelector("reads_1"))
                )
            )
        ]


class SelectorTestTool(CommandTool):
    """
    TODO
    AliasSelector?
    ForEachSelector?
    ResourceSelector?
    MemorySelector?
    CpuSelector?
    DiskSelector?
    TimeSelector ?   
    """
    
    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
        
    def base_command(self) -> Optional[str | list[str]]:
        return 'cat'

    def friendly_name(self):
        return "TEST: SelectorTestTool"

    def tool(self):
        return "SelectorTestTool"
    
    def inputs(self):
        return [
            ToolInput(
                "inFile",
                File(),
            ),
            ToolInput(
                "targetFilename",
                Filename(
                    prefix=j.InputSelector("inFile", remove_file_extension=True),  # wtf is this
                    suffix=".fastq",
                    extension=".gz",
                ),
                prefix="--targetFilename",
                position=1,
            ),
        ]
    
    def arguments(self):
        return [
            ToolArgument(
                j.InputSelector("inFile"),
                prefix="--inFile", 
                position=2,
            )
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out",
                File(),
                selector=j.WildcardSelector("*fastq.gz"),
            )
        ]


class StandardTestTool(CommandTool):
    
    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
        
    def base_command(self) -> Optional[str | list[str]]:
        return 'cat'

    def friendly_name(self):
        return "TEST: StandardTestTool"

    def tool(self):
        return "StandardTestTool"
    
    def inputs(self):
        return [
            ToolInput("inFile", File()),
            ToolInput("inFileArr", Array(File())),
            ToolInput("inBamBai", BamBai()),
            ToolInput("inBamBaiArr", Array(BamBai())),
            ToolInput("inStr", String()),
            ToolInput("inStrOpt", String(optional=True)),
            ToolInput("inInt", Int()),
            ToolInput("inIntOpt", Int(optional=True)),
        ]
    
    def arguments(self):
        return [
            # file operators
            ToolArgument(
                value=j.BasenameOperator(j.InputSelector("inFile")),
                prefix="--BasenameOperator",
                position=1,
            ),
            ToolArgument(
                value=j.DirnameOperator(j.InputSelector("inFile")),
                prefix="--DirnameOperator",
                position=2,
            ),
            ToolArgument(
                value=j.NamerootOperator(j.InputSelector("inFile")),
                prefix="--NamerootOperator",
                position=3,
            ),
            ToolArgument(
                value=j.NameextOperator(j.InputSelector("inFile")),
                prefix="--NameextOperator",
                position=4,
            ),
            ToolArgument(
                value=j.FileSizeOperator(j.InputSelector("inFile")),
                prefix="--FileSizeOperator",
                position=5,
            ),
            ToolArgument(
                value=j.ReadContents(j.InputSelector("inFile")),
                prefix="--ReadContents",
                position=6,
            ),
            ToolArgument(
                value=j.ReadJsonOperator(j.InputSelector("inFile")),
                prefix="--ReadJsonOperator",
                position=7,
            ),
            
            # int operators
            ToolArgument(
                value=j.RangeOperator(j.InputSelector("inInt")),
                prefix="--RangeOperator",
                position=8,
            ),

            # string operators
            ToolArgument(
                value=j.SplitOperator(j.InputSelector("inStr"), ' '),
                prefix="--SplitOperator",
                position=9,
            ),
            ToolArgument(
                value=j.ReplaceOperator(j.InputSelector("inStr"), '/[a-z]/', '/[A-Z]/'),
                prefix="--ReplaceOperator",
                position=10,
            ),

            # array operators
            ToolArgument(
                value=j.TransposeOperator(j.InputSelector("inFileArr")),
                prefix="--TransposeOperator",
                position=11,
            ),
            ToolArgument(
                value=j.LengthOperator(j.InputSelector("inFileArr")),
                prefix="--LengthOperator",
                position=12,
            ),
            ToolArgument(
                value=j.SliceOperator(j.InputSelector("inFileArr"), 0, -2),
                prefix="--SliceOperator",
                position=13,
            ),
            ToolArgument(
                value=j.FlattenOperator(j.InputSelector("inFileArr")),
                prefix="--FlattenOperator",
                position=13,
            ),
            ToolArgument(
                value=j.IndexOperator(j.InputSelector("inFileArr"), 1),
                prefix="--IndexOperator",
                position=14,
            ),
            ToolArgument(
                value=j.IndexOperator(j.InputSelector("inBamBaiArr"), 1),
                prefix="--IndexOperator",
                position=14,
            ),
            ToolArgument(
                value=j.ApplyPrefixOperator("--prefix", j.InputSelector("inFileArr")),
                prefix="--ApplyPrefixOperator",
                position=15,
            ),
            ToolArgument(
                value=j.FirstOperator(j.InputSelector("inFileArr")),
                prefix="--FirstOperator",
                position=16,
            ),
            ToolArgument(
                value=j.FilterNullOperator(j.InputSelector("inFileArr")),
                prefix="--FilterNullOperator",
                position=17,
            ),
            ToolArgument(
                value=j.FilterNullOperator([1,None,3,4,None]),
                prefix="--FilterNullOperatorLiteral",
                position=17,
            ),
            ToolArgument(
                value=j.JoinOperator(j.InputSelector("inFileArr"), ','),
                prefix="--JoinOperator",
                position=18,
            ),
        ]

    def outputs(self):
        return []


class LogicalTestTool(CommandTool):

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
        
    def base_command(self) -> Optional[str | list[str]]:
        return 'cat'

    def friendly_name(self):
        return "TEST: LogicalTestTool"

    def tool(self):
        return "LogicalTestTool"

    def inputs(self):
        return [
            ToolInput("inFile", File()),
            ToolInput("inFileOpt", File(optional=True)),
            ToolInput("inFileArr", Array(File())),
            ToolInput("inBamBai", BamBai()),
            ToolInput("inBamBaiArr", Array(BamBai())),
            ToolInput("inStr", String()),
            ToolInput("inStrOpt", String(optional=True)),
            ToolInput("inInt", Int()),
            ToolInput("inIntOpt", Int(optional=True)),
            ToolInput("inFloat", Float()),
            ToolInput("inFloatOpt", Float(optional=True)),
        ]
    
    def arguments(self):
        return [
            # single value 
            ToolArgument(
                j.AssertNotNull(j.InputSelector("inFileOpt")),
                prefix="--AssertNotNull",
                position=1,
            ),
            ToolArgument(
                j.IsDefined(j.InputSelector("inFileOpt")),
                prefix="--IsDefined",
                position=2,
            ),
            ToolArgument(
                j.NotOperator(j.InputSelector("inFileOpt")),
                prefix="--NotOperator",
                position=3,
            ),
            ToolArgument(
                j.FloorOperator(j.InputSelector("inFloat")),
                prefix="--FloorOperator",
                position=4,
            ),
            ToolArgument(
                j.CeilOperator(j.InputSelector("inFloat")),
                prefix="--CeilOperator",
                position=5,
            ),
            ToolArgument(
                j.RoundOperator(j.InputSelector("inFloat")),
                prefix="--RoundOperator",
                position=6,
            ),
            ToolArgument(
                j.GroupOperator(j.AddOperator(5, 10)),
                prefix="--GroupOperator",
                position=7,
            ),

            # two+ value 
            ToolArgument(
                j.If(j.IsDefined(j.InputSelector("inFileOpt")), j.InputSelector("inFileOpt"), ""),
                prefix="--If",
                position=8,
            ),
            ToolArgument(
                j.AndOperator(j.IsDefined(j.InputSelector("inFileOpt")), j.IsDefined(j.InputSelector("inStrOpt"))),
                prefix="--AndOperator",
                position=9,
            ),
            ToolArgument(
                j.OrOperator(j.IsDefined(j.InputSelector("inFileOpt")), j.IsDefined(j.InputSelector("inStrOpt"))),
                prefix="--OrOperator",
                position=10,
            ),
            ToolArgument(
                j.EqualityOperator(j.InputSelector("inStrOpt"), "hello!"),
                prefix="--EqualityOperator",
                position=11,
            ),
            ToolArgument(
                j.InequalityOperator(j.InputSelector("inStrOpt"), "hello!"),
                prefix="--InequalityOperator",
                position=12,
            ),
            ToolArgument(
                j.GtOperator(j.InputSelector("inInt"), 0),
                prefix="--GtOperator",
                position=13,
            ),
            ToolArgument(
                j.GteOperator(j.InputSelector("inInt"), 0),
                prefix="--GteOperator",
                position=14,
            ),
            ToolArgument(
                j.LtOperator(j.InputSelector("inInt"), 0),
                prefix="--LtOperator",
                position=15,
            ),
            ToolArgument(
                j.LteOperator(j.InputSelector("inInt"), 0),
                prefix="--LteOperator",
                position=16,
            ),
            ToolArgument(
                j.AddOperator(9, 10),
                prefix="--AddOperator",
                position=17,
            ),
            ToolArgument(
                j.SubtractOperator(9, 10),
                prefix="--SubtractOperator",
                position=18,
            ),
            ToolArgument(
                j.MultiplyOperator(9, 10),
                prefix="--MultiplyOperator",
                position=19,
            ),
            ToolArgument(
                j.DivideOperator(9, 10),
                prefix="--DivideOperator",
                position=20,
            ),
        ]

    def outputs(self):
        return []


class ComplexOperatorTestTool(CommandTool):
    
    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"
        
    def base_command(self) -> Optional[str | list[str]]:
        return 'cat'

    def friendly_name(self):
        return "TEST: ComplexOperatorTestTool"

    def tool(self):
        return "ComplexOperatorTestTool"
    
    def inputs(self):
        return [
            ToolInput("inFile", File()),
            ToolInput("inFileOpt", File(optional=True)),
            ToolInput("inFileArr", Array(File())),
            ToolInput("inBamBai", BamBai()),
            ToolInput("inBamBaiArr", Array(BamBai())),
            ToolInput("inStr", String()),
            ToolInput("inStrOpt", String(optional=True)),
            ToolInput("inInt", Int()),
            ToolInput("inIntOpt", Int(optional=True)),
            ToolInput("inFloat", Float()),
            ToolInput("inFloatOpt", Float(optional=True)),
        ]
    
    def arguments(self):
        return [
            # $(inputs.reads_1.basename).trimmed$(inputs.reads_1.nameext)
            ToolArgument(
                value=j.StringFormatter(
                    format='{token1}.trimmed{token2}',
                    token1=j.BasenameOperator(j.InputSelector("inFile")),
                    token2=j.NameextOperator(j.InputSelector("inFile")),
                ),
                prefix="--StringFormatter",
                position=1,
            ),

            # $(inputs.molecule_info_h5[0].basename.split('.').slice(0,-1).join('.'))
            ToolArgument(
                value=j.JoinOperator(
                        j.SliceOperator(
                            j.SplitOperator(
                                j.BasenameOperator(j.IndexOperator(j.InputSelector("inFileArr"), 0)), '.'),
                            0, -1), 
                        '.'),
                prefix="--ArrayMethodChain",
                position=2,
            ),

            # $('> ' + inputs.index_base_name + '.log')
            ToolArgument(
                value=j.AddOperator(j.AddOperator('> ', j.InputSelector("inStr")), '.log'),
                prefix="--Concat",
                position=3,
            ),
            
            # $(Math.ceil(inputs.target.size/(1024*1024*1024) + 20))
            ToolArgument(
                value=j.CeilOperator(
                        j.AddOperator(
                            j.DivideOperator(
                                j.FileSizeOperator(j.InputSelector("inFile")),
                                j.GroupOperator(
                                    j.MultiplyOperator(
                                        j.MultiplyOperator(1024, 1024),
                                        1024,
                                    )
                                )
                            ), 
                            20),
                        ),
                prefix="--Math",
                position=4,
            ),
        ]

    def outputs(self):
        return []


