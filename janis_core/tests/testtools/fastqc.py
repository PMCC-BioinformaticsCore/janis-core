
from janis_core.types.common_data_types import Boolean
from janis_core.types.common_data_types import String
from janis_core.types.common_data_types import File
from janis_core.types.common_data_types import Int
from janis_core import CommandToolBuilder
from janis_core import WildcardSelector
from janis_core import ToolMetadata
from janis_core import ToolOutput
from janis_core import ToolInput


metadata = ToolMetadata(
    short_documentation="Read Quality reports",
    keywords=[],
    contributors=['gxtool2janis'],
    dateCreated="2022-08-29 10:49:40",
    dateUpdated="2022-08-29 10:49:40",
    version="0.72",
    doi="None",
    citation="tool xml missing citation",
    documentationUrl=None,
    documentation="""
.. class:: infomark

**Purpose**

FastQC aims to provide a simple way to do some quality control checks on raw
sequence data coming from high throughput sequencing pipelines.
It provides a set of analyses which you can use to get a quick
impression of whether your data has any problems of
which you should be aware before doing any further analysis.

The main functions of FastQC are:

- Import of data from BAM, SAM or FastQ/FastQ.gz files (any variant),
- Providing a quick overview to tell you in which areas there may be problems
- Summary graphs and tables to quickly assess your data
- Export of results to an HTML based permanent report
- Offline operation to allow automated generation of reports without running the interactive application

-----

.. class:: infomark

**FastQC**

This is a Galaxy wrapper. It merely exposes the external package FastQC_ which is documented at FastQC_
Kindly acknowledge it as well as this tool if you use it.
FastQC incorporates the Picard-tools_ libraries for SAM/BAM processing.

The contaminants file parameter was borrowed from the independently developed
fastqcwrapper contributed to the Galaxy Community Tool Shed by J. Johnson.
Adaption to version 0.11.2 by T. McGowan.

-----

.. class:: infomark

**Inputs and outputs**

FastQC_ is the best place to look for documentation - it's very good.
A summary follows below for those in a tearing hurry.

This wrapper will accept a Galaxy fastq, fastq.gz, sam or bam as the input read file to check.
It will also take an optional file containing a list of contaminants information, in the form of
a tab-delimited file with 2 columns, name and sequence. As another option the tool takes a custom
limits.txt file that allows setting the warning thresholds for the different modules and also specifies
which modules to include in the output.

The tool produces a basic text and a HTML output file that contain all of the results, including the following:

- Basic Statistics
- Per base sequence quality
- Per sequence quality scores
- Per base sequence content
- Per base GC content
- Per sequence GC content
- Per base N content
- Sequence Length Distribution
- Sequence Duplication Levels
- Overrepresented sequences
- Kmer Content

All except Basic Statistics and Overrepresented sequences are plots.
 .. _FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
 .. _Picard-tools: https://broadinstitute.github.io/picard/
    """
)

mandatory_inputs = [
	# Positionals
	ToolInput(
		'inputFile',
		File,
		position=2,
		default=None,
		doc="Short read data from your current history",
	),

	# Flags
	ToolInput(
		'extract',
		Boolean(optional=True),
		prefix='--extract',
		position=1,
		default="False",
	),
	ToolInput(
		'nogroup',
		Boolean(optional=True),
		prefix='--nogroup',
		position=1,
		default="False",
	),
	ToolInput(
		'quiet',
		Boolean(optional=True),
		prefix='--quiet',
		position=1,
		default="False",
	),

	# Options
	ToolInput(
		'adapters',
		File(optional=True),
		prefix='--adapters',
		position=1,
		default=None,
		doc="Adapter list. list of adapters adapter sequences which will be explicity searched against the library. tab delimited file with 2 columns: name and sequence.",
	),
	ToolInput(
		'contaminants',
		File(optional=True),
		prefix='--contaminants',
		position=1,
		default=None,
		doc="Contaminant list. tab delimited file with 2 columns: name and sequence.  For example: Illumina Small RNA RT Primer CAAGCAGAAGACGGCATACGA",
	),
	ToolInput(
		'kmers',
		Int,
		prefix='--kmers',
		position=1,
		default=7,
		doc="length of Kmer to look for. note: the Kmer test is disabled and needs to be enabled using a custom Submodule and limits file",
	),
	ToolInput(
		'limits',
		File(optional=True),
		prefix='--limits',
		position=1,
		default=None,
		doc="Submodule and Limit specifing file. a file that specifies which submodules are to be executed (default=all) and also specifies the thresholds for the each submodules warning parameter",
	),
	ToolInput(
		'minLength',
		Int(optional=True),
		prefix='--min_length',
		position=1,
		default=None,
		doc="Lower limit on the length of the sequence to be shown in the report.  As long as you set this to a value greater or equal to your longest read length then this will be the sequence length used to create your read groups.  This can be useful for making directly comaparable statistics from datasets with somewhat variable read lengths.",
	),
	ToolInput(
		'optionF',
		String(optional=True),
		prefix='-f',
		position=1,
		default=None,
	),
	ToolInput(
		'outdir',
		String(optional=True),
		prefix='--outdir',
		position=1,
		default=None,
	),

]

outputs = [
		ToolOutput(
		'outHtmlFile',
		File,
		selector=WildcardSelector(r"output.html"),
		doc="Webpage",
	),
		ToolOutput(
		'outTextFile',
		File,
		selector=WildcardSelector(r"output.txt"),
		doc="RawData",
	),

]

FastqcTestTool = CommandToolBuilder(
    tool="fastqc",
    version="0.72",
    metadata=metadata,
    container="quay.io/biocontainers/fastqc:0.11.8--2",
    base_command=['fastqc'],
    inputs=mandatory_inputs,
    outputs=outputs
)


if __name__ == "__main__":
    FastqcTestTool().translate(
        "wdl", to_console=True
    )


