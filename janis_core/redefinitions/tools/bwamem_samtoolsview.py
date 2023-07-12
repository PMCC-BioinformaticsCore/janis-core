from typing import List, Dict, Any
from datetime import datetime
from janis_core import (
    ToolInput,
    Int,
    Float,
    Boolean,
    String,
    ToolOutput,
    Filename,
    File,
    InputSelector,
    ToolArgument,
    CaptureType,
    CpuSelector,
    Array,
    StringFormatter,
    ToolMetadata,
    get_value_for_hints_and_ordered_resource_tuple
)

from ..types import FastaWithDict, FastqGzPair, Bam, Bed
from ..tools import BioinformaticsTool

BWA_MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 8,
            CaptureType.EXOME: 12,
            CaptureType.CHROMOSOME: 12,
            CaptureType.THIRTYX: 16,
            CaptureType.NINETYX: 20,
            CaptureType.THREEHUNDREDX: 24,
        },
    )
]

BWA_CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.TARGETED: 16,
            CaptureType.EXOME: 20,
            CaptureType.CHROMOSOME: 24,
            CaptureType.THIRTYX: 40,
            CaptureType.NINETYX: 40,
            CaptureType.THREEHUNDREDX: 40,
        },
    )
]


class BwaMem_SamToolsView(BioinformaticsTool):
    
    def tool(self) -> str:
        return "BwaMemSamtoolsView"

    def tool_provider(self):
        return "common"

    def version(self):
        return "0.7.17|1.9"

    def container(self):
        return "michaelfranklin/bwasamtools:0.7.17-1.9"

    def base_command(self):
        return None

    def arguments(self):
        return [
            ToolArgument("bwa", position=0, shell_quote=False),
            ToolArgument("mem", position=1, shell_quote=False),
            ToolArgument("|", position=5, shell_quote=False),
            ToolArgument("samtools", position=6, shell_quote=False),
            ToolArgument("view", position=7, shell_quote=False),
            ToolArgument(
                InputSelector("reference"), prefix="-T", position=8, shell_quote=False
            ),
            ToolArgument(
                CpuSelector(),
                position=8,
                shell_quote=False,
                prefix="--threads",
                doc="(-@)  Number of additional threads to use [0]",
            ),
            ToolArgument(
                "-h",
                position=8,
                shell_quote=False,
                doc="Include the header in the output.",
            ),
            ToolArgument(
                "-b", position=8, shell_quote=False, doc="Output in the BAM format."
            ),
            ToolArgument(
                StringFormatter(
                    "@RG\\tID:{name}\\tSM:{name}\\tLB:{name}\\tPL:{pl}",
                    name=InputSelector("sampleName"),
                    pl=InputSelector("platformTechnology"),
                ),
                prefix="-R",
                position=2,
                doc="Complete read group header line. ’\\t’ can be used in STR and will be converted to a TAB"
                "in the output SAM. The read group ID will be attached to every read in the output. "
                "An example is ’@RG\\tID:foo\\tSM:bar’. (Default=null) "
                "https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups",
            ),
            ToolArgument(
                CpuSelector(),
                prefix="-t",
                position=2,
                shell_quote=False,
                doc="Number of threads. (default = 1)",
            ),
        ]

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("reference", FastaWithDict(), position=2, shell_quote=False),
            ToolInput("reads", FastqGzPair, position=3, shell_quote=False, doc=None),
            ToolInput(
                "mates",
                FastqGzPair(optional=True),
                separator=" ",
                position=4,
                shell_quote=False,
                doc=None,
            ),
            ToolInput(
                "outputFilename",
                Filename(prefix=InputSelector("sampleName"), extension=".bam"),
                position=8,
                shell_quote=False,
                prefix="-o",
                doc="output file name [stdout]",
            ),
            # Eventually it would be cool to have like a cascading:
            #   - If readGroupHeaderLine provided, use that,
            #   - If sampleName provided, construct based on that
            #   - Else don't include
            # but this is probbaly a bit hard to do, and for all our purposes we require a readGroupHeaderLine,
            # so we're always going to construct it:
            ToolInput(
                "sampleName",
                String(),
                doc="Used to construct the readGroupHeaderLine with format: "
                "'@RG\\tID:{name}\\tSM:{name}\\tLB:{name}\\tPL:ILLUMINA'",
            ),
            ToolInput(
                "platformTechnology",
                String(optional=True),
                doc="(ReadGroup: PL) Used to construct the readGroupHeaderLine, defaults: ILLUMINA",
                default="ILLUMINA",
            ),
            *self.bwa_additional_inputs,
            *self.samtools_additional_args,
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", Bam(), glob=InputSelector("outputFilename")),
            # ToolOutput("skippedReads", File(optional=True), glob=InputSelector("skippedReadsOutputFilename"))
        ]

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, BWA_MEM_TUPLE)
        if val:
            return val
        return 16

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, BWA_CORES_TUPLE)
        if val:
            return val
        return 16

    def friendly_name(self) -> str:
        return "Bwa mem + Samtools View"

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=datetime(2019, 5, 10),
            dateUpdated=datetime(2020, 7, 14),
        )

    bwa_additional_inputs = [
        ToolInput(
            "minimumSeedLength",
            Int(optional=True),
            prefix="-k",
            position=2,
            shell_quote=False,
            doc="Matches shorter than INT will be missed. The alignment speed is usually "
            "insensitive to this value unless it significantly deviates 20. (Default: 19)",
        ),
        ToolInput(
            "batchSize",
            Int(optional=True),
            prefix="-K",
            position=2,
            shell_quote=False,
            doc="Process INT input bases in each batch regardless of the number of threads in use [10000000*nThreads]. "
            "By default, the batch size is proportional to the number of threads in use. Because the inferred "
            "insert size distribution slightly depends on the batch size, using different number of threads "
            "may produce different output. Specifying this option helps reproducibility.",
        ),
        ToolInput(
            "useSoftClippingForSupplementaryAlignments",
            Boolean(optional=True),
            prefix="-Y",
            position=2,
            shell_quote=False,
            doc="Use soft clipping CIGAR operation for supplementary alignments. By default, BWA-MEM uses soft "
            "clipping for the primary alignment and hard clipping for supplementary alignments.",
        ),
        ToolInput(
            "bandwidth",
            Int(optional=True),
            prefix="-w",
            position=2,
            shell_quote=False,
            doc="Essentially, gaps longer than ${bandWidth} will not be found. Note that the maximum gap length "
            "is also affected by the scoring matrix and the hit length, not solely determined by this option."
            " (Default: 100)",
        ),
        ToolInput(
            "offDiagonalXDropoff",
            Int(optional=True),
            prefix="-d",
            position=2,
            shell_quote=False,
            doc="(Z-dropoff): Stop extension when the difference between the best and the current extension "
            "score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, "
            "respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except "
            "that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only "
            "avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. "
            "(Default: 100)",
        ),
        ToolInput(
            "reseedTrigger",
            Float(optional=True),
            prefix="-r",
            position=2,
            shell_quote=False,
            doc="Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter "
            "for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment "
            "speed but lower accuracy. (Default: 1.5)",
        ),
        ToolInput(
            "occurenceDiscard",
            Int(optional=True),
            prefix="-c",
            position=2,
            shell_quote=False,
            doc="Discard a MEM if it has more than INT occurence in the genome. "
            "This is an insensitive parameter. (Default: 10000)",
        ),
        ToolInput(
            "performSW",
            Boolean(optional=True),
            prefix="-P",
            position=2,
            shell_quote=False,
            doc="In the paired-end mode, perform SW to rescue missing hits only but "
            "do not try to find hits that fit a proper pair.",
        ),
        ToolInput(
            "matchingScore",
            Int(optional=True),
            prefix="-A",
            position=2,
            shell_quote=False,
            doc="Matching score. (Default: 1)",
        ),
        ToolInput(
            "mismatchPenalty",
            Int(optional=True),
            prefix="-B",
            position=2,
            shell_quote=False,
            doc="Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. "
            "(Default: 4)",
        ),
        ToolInput(
            "openGapPenalty",
            Int(optional=True),
            prefix="-O",
            position=2,
            shell_quote=False,
            doc="Gap open penalty. (Default: 6)",
        ),
        ToolInput(
            "gapExtensionPenalty",
            Int(optional=True),
            prefix="-E",
            position=2,
            shell_quote=False,
            doc="Gap extension penalty. A gap of length k costs O + k*E "
            "(i.e. -O is for opening a zero-length gap). (Default: 1)",
        ),
        ToolInput(
            "clippingPenalty",
            Int(optional=True),
            prefix="-L",
            position=2,
            shell_quote=False,
            doc="Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score "
            "reaching the end of query. If this score is larger than the best SW score minus the "
            "clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag "
            "reports the best SW score; clipping penalty is not deducted. (Default: 5)",
        ),
        ToolInput(
            "unpairedReadPenalty",
            Int(optional=True),
            prefix="-U",
            position=2,
            shell_quote=False,
            doc="Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as "
            "scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. "
            "It compares these two scores to determine whether we should force pairing. (Default: 9)",
        ),
        ToolInput(
            "assumeInterleavedFirstInput",
            Boolean(optional=True),
            prefix="-p",
            position=2,
            shell_quote=False,
            doc="Assume the first input query file is interleaved paired-end FASTA/Q. ",
        ),
        ToolInput(
            "outputAlignmentThreshold",
            Int(optional=True),
            prefix="-T",
            position=2,
            shell_quote=False,
            doc="Don’t output alignment with score lower than INT. Only affects output. (Default: 30)",
        ),
        ToolInput(
            "outputAllElements",
            Boolean(optional=True),
            prefix="-a",
            position=2,
            shell_quote=False,
            doc="Output all found alignments for single-end or unpaired paired-end reads. "
            "These alignments will be flagged as secondary alignments.",
        ),
        ToolInput(
            "appendComments",
            Boolean(optional=True),
            prefix="-C",
            position=2,
            shell_quote=False,
            doc="Append append FASTA/Q comment to SAM output. This option can be used to transfer "
            "read meta information (e.g. barcode) to the SAM output. Note that the FASTA/Q comment "
            "(the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). "
            "Malformated comments lead to incorrect SAM output.",
        ),
        ToolInput(
            "hardClipping",
            Boolean(optional=True),
            prefix="-H",
            position=2,
            shell_quote=False,
            doc="Use hard clipping ’H’ in the SAM output. This option may dramatically reduce "
            "the redundancy of output when mapping long contig or BAC sequences.",
        ),
        ToolInput(
            "markShorterSplits",
            Boolean(optional=True),
            prefix="-M",
            position=2,
            shell_quote=False,
            doc="Mark shorter split hits as secondary (for Picard compatibility).",
        ),
        ToolInput(
            "verboseLevel",
            Int(optional=True),
            prefix="-v",
            position=2,
            shell_quote=False,
            doc="Control the verbose level of the output. "
            "This option has not been fully supported throughout BWA. Ideally, a value: "
            "0 for disabling all the output to stderr; "
            "1 for outputting errors only; "
            "2 for warnings and errors; "
            "3 for all normal messages; "
            "4 or higher for debugging. When this option takes value 4, the output is not SAM. (Default: 3)",
        ),
    ]

    samtools_additional_args = [
        ToolInput(
            "skippedReadsOutputFilename",
            String(optional=True),
            position=8,
            shell_quote=False,
            prefix="-U",
            doc="output reads not selected by filters to FILE [null]",
        ),
        ToolInput(
            "referenceIndex",
            File(optional=True),
            position=8,
            shell_quote=False,
            prefix="-t",
            doc="FILE listing reference names and lengths (see long help) [null]",
        ),
        ToolInput(
            "intervals",
            Bed(optional=True),
            position=8,
            shell_quote=False,
            prefix="-L",
            doc="only include reads overlapping this BED FILE [null]",
        ),
        ToolInput(
            "includeReadsInReadGroup",
            String(optional=True),
            position=8,
            shell_quote=False,
            prefix="-r",
            doc="only include reads in read group STR [null]",
        ),
        ToolInput(
            "includeReadsInFile",
            File(optional=True),
            position=8,
            shell_quote=False,
            prefix="-R",
            doc="only include reads with read group listed in FILE [null]",
        ),
        ToolInput(
            "includeReadsWithQuality",
            Int(optional=True),
            position=8,
            shell_quote=False,
            prefix="-q",
            doc="only include reads with mapping quality >= INT [0]",
        ),
        ToolInput(
            "includeReadsInLibrary",
            String(optional=True),
            position=8,
            shell_quote=False,
            prefix="-l",
            doc="only include reads in library STR [null]",
        ),
        ToolInput(
            "includeReadsWithCIGAROps",
            Int(optional=True),
            position=8,
            shell_quote=False,
            prefix="-m",
            doc="only include reads with number of CIGAR operations consuming query sequence >= INT [0]",
        ),
        ToolInput(
            "includeReadsWithAllFLAGs",
            Array(Int(), optional=True),
            position=8,
            shell_quote=False,
            prefix="-f",
            separator=" ",
            doc="only include reads with all of the FLAGs in INT present [0]",
        ),
        ToolInput(
            "includeReadsWithoutFLAGs",
            Array(Int(), optional=True),
            position=8,
            shell_quote=False,
            prefix="-F",
            separator=" ",
            doc="only include reads with none of the FLAGS in INT present [0]",
        ),
        ToolInput(
            "excludeReadsWithAllFLAGs",
            Array(Int(), optional=True),
            position=8,
            shell_quote=False,
            prefix="-G",
            separator=" ",
            doc="only EXCLUDE reads with all of the FLAGs in INT present [0] "
            "fraction of templates/read pairs to keep; INT part sets seed)",
        ),
        ToolInput(
            "useMultiRegionIterator",
            Boolean(optional=True),
            position=8,
            shell_quote=False,
            prefix="-M",
            doc="use the multi-region iterator (increases the speed, removes "
            "duplicates and outputs the reads as they are ordered in the file)",
        ),
        ToolInput(
            "readTagToStrip",
            String(optional=True),
            position=8,
            shell_quote=False,
            prefix="-x",
            doc="read tag to strip (repeatable) [null]",
        ),
        ToolInput(
            "collapseBackwardCIGAROps",
            Boolean(optional=True),
            position=8,
            shell_quote=False,
            prefix="-B",
            doc="collapse the backward CIGAR operation Specify a single input file format "
            "option in the form of OPTION or OPTION=VALUE",
        ),
        ToolInput(
            "outputFmt",
            String(optional=True),
            position=8,
            shell_quote=False,
            prefix="--output-fmt",
            doc="(OPT[, -O)  Specify output format (SAM, BAM, CRAM) Specify a single "
            "output file format option in the form of OPTION or OPTION=VALUE",
        ),
    ]


if __name__ == "__main__":
    print(BwaMem_SamToolsView().translate("wdl"))
