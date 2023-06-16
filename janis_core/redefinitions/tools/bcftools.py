
from abc import ABC
from datetime import datetime, date
from typing import Dict, Any

from janis_core import (
    ToolInput,
    ToolOutput,
    String,
    Boolean,
    InputSelector,
    File,
    Array,
    Int,
    Filename,
    ToolMetadata,
    CaptureType,
    UnionType,
)
from janis_core import get_value_for_hints_and_ordered_resource_tuple
from janis_core.tool.test_classes import TTestCase

from ..types import Vcf, CompressedVcf
from ..tools import BioinformaticsTool



### BASE CLASSES ###

class BcfToolsToolBase(BioinformaticsTool, ABC):
    def tool_provider(self):
        return "bcftools"
    

class BcfTools_1_5:
    def container(self):
        # Todo: Create a docker container with 1_9
        return "biocontainers/bcftools:v1.5_cv2"

    def version(self):
        return "v1.5"


class BcfTools_1_9:
    def container(self):
        return "biocontainers/bcftools:v1.9-1-deb_cv1"

    def version(self):
        return "v1.9"


class BcfTools_1_12:
    def container(self):
        return "quay.io/biocontainers/bcftools:1.12--h45bccc9_1"

    def version(self):
        return "v1.12"



### BCF_TOOLS_ANNOTATE ###

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 1,
            CaptureType.EXOME: 1,
            CaptureType.THIRTYX: 1,
            CaptureType.NINETYX: 1,
            CaptureType.THREEHUNDREDX: 1,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 16,
            CaptureType.EXOME: 16,
            CaptureType.THIRTYX: 16,
            CaptureType.NINETYX: 64,
            CaptureType.THREEHUNDREDX: 64,
        },
    )
]


class BcfToolsAnnotateBase(BcfToolsToolBase, ABC):
    def bind_metadata(self):
        self.metadata.contributors = ["Michael Franklin"]
        self.metadata.dateCreated = date(2019, 1, 24)
        self.metadata.dateUpdated = date(2019, 1, 24)
        self.metadata.doi = "http://www.ncbi.nlm.nih.gov/pubmed/19505943"
        self.metadata.citation = (
            "Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, "
            "and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) "
            "format and SAMtools, Bioinformatics (2009) 25(16) 2078-9"
        )
        self.metadata.documentationUrl = (
            "https://samtools.github.io/bcftools/bcftools.html#annotate"
        )
        self.metadata.documentation = (
            self.metadata.documentation if self.metadata.documentation else ""
        ) + "------------------------------------\n\nAdd or remove annotations."

    def tool(self):
        return "bcftoolsAnnotate"

    def friendly_name(self):
        return "BCFTools: Annotate"

    def base_command(self):
        return ["bcftools", "annotate"]

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 1

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 8

    def inputs(self):
        return [
            ToolInput("vcf", Vcf(), position=10),
            ToolInput(
                "outputFilename",
                Filename(extension=".vcf"),
                prefix="--output",
                doc="[-o] see Common Options",
            ),
            *self.additional_args,
        ]

    def outputs(self):
        return [ToolOutput("out", Vcf, glob=InputSelector("outputFilename"))]

    def docurl():
        return "https://samtools.github.io/bcftools/bcftools.html#annotate"

    additional_args = [
        ToolInput(
            "annotations",
            File(optional=True),
            prefix="--annotations",
            doc="[-a] Bgzip-compressed and tabix-indexed file with annotations. The file can be VCF, BED, "
            "or a tab-delimited file with mandatory columns CHROM, POS (or, alternatively, FROM and TO), "
            "optional columns REF and ALT, and arbitrary number of annotation columns. BED files are "
            'expected to have the ".bed" or ".bed.gz" suffix (case-insensitive), otherwise a '
            "tab-delimited file is assumed. Note that in case of tab-delimited file, the coordinates POS, "
            "FROM and TO are one-based and inclusive. When REF and ALT are present, only matching VCF "
            "records will be annotated. When multiple ALT alleles are present in the annotation file "
            "(given as comma-separated list of alleles), at least one must match one of the alleles in "
            "the corresponding VCF record. Similarly, at least one alternate allele from a multi-allelic "
            "VCF record must be present in the annotation file. Missing values can be added by "
            'providing "." in place of actual value. Note that flag types, such as "INFO/FLAG", '
            'can be annotated by including a field with the value "1" to set the flag, "0" to remove '
            'it, or "." to keep existing flags. See also -c, --columns and -h, --header-lines.',
        ),
        ToolInput(
            "collapse",
            String(optional=True),
            prefix="--collapse",
            doc="(snps|indels|both|all|some|none) Controls how to match records from the annotation file to "
            "the target VCF. Effective only when -a is a VCF or BCF. See Common Options for more.",
        ),
        ToolInput(
            "columns",
            Array(String(), optional=True),
            prefix="--columns",
            doc="[-c] Comma-separated list of columns or tags to carry over from the annotation file "
            "(see also -a, --annotations). If the annotation file is not a VCF/BCF, list describes "
            "the columns of the annotation file and must include CHROM, POS "
            "(or, alternatively, FROM and TO), and optionally REF and ALT. Unused columns which "
            'should be ignored can be indicated by "-". If the annotation file is a VCF/BCF, '
            "only the edited columns/tags must be present and their order does not matter. "
            "The columns ID, QUAL, FILTER, INFO and FORMAT can be edited, where INFO tags "
            'can be written both as "INFO/TAG" or simply "TAG", and FORMAT tags can be written '
            'as "FORMAT/TAG" or "FMT/TAG". The imported VCF annotations can be renamed as '
            '"DST_TAG:=SRC_TAG" or "FMT/DST_TAG:=FMT/SRC_TAG". To carry over all INFO annotations, '
            'use "INFO". To add all INFO annotations except "TAG", use "^INFO/TAG". '
            "By default, existing values are replaced. To add annotations without overwriting "
            "existing values (that is, to add missing tags or add values to existing tags with "
            'missing values), use "+TAG" instead of "TAG". To append to existing values (rather '
            'than replacing or leaving untouched), use "=TAG" (instead of "TAG" or "+TAG"). '
            'To replace only existing values without modifying missing annotations, use "-TAG". '
            "If the annotation file is not a VCF/BCF, all new annotations must be "
            "defined via -h, --header-lines.",
        ),
        ToolInput(
            "exclude",
            String(optional=True),
            prefix="--exclude",
            doc="[-e] exclude sites for which EXPRESSION is true. For valid expressions see EXPRESSIONS.",
        ),
        ToolInput(
            "headerLines",
            File(optional=True),
            prefix="--header-lines",
            doc="[-h] Lines to append to the VCF header, see also -c, --columns and -a, --annotations.",
        ),
        ToolInput(
            "setId",
            String(optional=True),
            prefix="--set-id",
            doc="[-I] assign ID on the fly. The format is the same as in the query command (see below). "
            'By default all existing IDs are replaced. If the format string is preceded by "+", only'
            " missing IDs will be set. For example, one can use "
            "# bcftools annotate --set-id +' % CHROM\_ % POS\_ % REF\_ % FIRST_ALT' file.vcf",
        ),
        ToolInput(
            "include",
            String(optional=True),
            prefix="--include",
            doc="[-i] include only sites for which EXPRESSION is true. For valid expressions see EXPRESSIONS.",
        ),
        ToolInput(
            "keepSites",
            Boolean(optional=True),
            prefix="--keep-sites",
            doc="keep sites wich do not pass -i and -e expressions instead of discarding them(",
        ),
        ToolInput(
            "markSites",
            String(optional=True),
            prefix="--mark-sites",
            doc='[-m] (+|-)annotate sites which are present ("+") or absent ("-") in the -a file with a '
            "new INFO/TAG flag",
        ),
        ToolInput(
            "outputType",
            String(optional=True),
            prefix="--output-type",
            doc="[-O] (b|u|z|v) see Common Options",
        ),
        ToolInput(
            "regions",
            String(optional=True),
            prefix="--regions",
            doc="([-r] chr|chr:pos|chr:from-to|chr:from-[,…]) see Common Options",
        ),
        ToolInput(
            "regionsFile",
            File(optional=True),
            prefix="--regions-file",
            doc="[-R] see Common Options",
        ),
        ToolInput(
            "renameChrs",
            File(optional=True),
            prefix="--rename-chrs",
            doc='rename chromosomes according to the map in file, with "old_name new_name\\n" pairs '
            "separated by whitespaces, each on a separate line.",
        ),
        ToolInput(
            "samples",
            Array(File(), optional=True),
            prefix="--samples",
            doc="[-s] subset of samples to annotate, see also Common Options",
        ),
        ToolInput(
            "samplesFile",
            File(optional=True),
            prefix="--samples-file",
            doc="[-S] subset of samples to annotate. If the samples are named differently in the target "
            'VCF and the -a, --annotations VCF, the name mapping can be given as "src_name dst_name\\n", '
            "separated by whitespaces, each pair on a separate line.",
        ),
        ToolInput(
            "threads", Int(optional=True), prefix="--threads", doc="see Common Options"
        ),
        ToolInput(
            "remove",
            Array(String(), optional=True),
            prefix="--remove",
            doc='[-x] List of annotations to remove. Use "FILTER" to remove all filters or "FILTER/SomeFilter" '
            'to remove a specific filter. Similarly, "INFO" can be used to remove all INFO tags and '
            '"FORMAT" to remove all FORMAT tags except GT. To remove all INFO tags except "FOO" and '
            '"BAR", use "^INFO/FOO,INFO/BAR" (and similarly for FORMAT and FILTER). "INFO" can be '
            'abbreviated to "INF" and "FORMAT" to "FMT".',
        ),
    ]


class BcfToolsAnnotate_1_9(BcfTools_1_9, BcfToolsAnnotateBase):
    pass


class BcfToolsAnnotate_1_5(BcfTools_1_5, BcfToolsAnnotateBase):
    pass


BcfToolsAnnotateLatest = BcfToolsAnnotate_1_9



### BCF_TOOLS_CONCAT ###


class BcfToolsConcatBase(BcfToolsToolBase, ABC):
    def tool(self):
        return "bcftoolsConcat"

    def friendly_name(self):
        return "BCFTools: Concat"

    def base_command(self):
        return ["bcftools", "concat"]

    def inputs(self):
        return [
            ToolInput("vcf", Array(UnionType(Vcf, CompressedVcf)), position=15),
            ToolInput(
                "outputFilename",
                Filename(extension=".vcf.gz"),
                prefix="-o",
                doc="--output: When output consists of a single stream, "
                "write it to FILE rather than to standard output, where it is written by default.",
            ),
            *self.additional_args,
        ]

    def outputs(self):
        return [
            ToolOutput("out", CompressedVcf(), glob=InputSelector("outputFilename"))
        ]

    def bind_metadata(self):
        from datetime import date

        self.metadata.contributors = ["Michael Franklin"]
        self.metadata.dateCreated = date(2019, 9, 9)
        self.metadata.dateUpdated = date(2019, 9, 9)
        self.metadata.doi = "http://www.ncbi.nlm.nih.gov/pubmed/19505943"
        self.metadata.citation = (
            "Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, "
            "and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) "
            "format and SAMtools, Bioinformatics (2009) 25(16) 2078-9"
        )
        self.metadata.documentationUrl = (
            "https://samtools.github.io/bcftools/bcftools.html#concat"
        )
        self.metadata.documentation = """
Concatenate or combine VCF/BCF files. All source files must have the same sample
columns appearing in the same order. The program can be used, for example, to
concatenate chromosome VCFs into one VCF, or combine a SNP VCF and an indel
VCF into one. The input files must be sorted by chr and position. The files
must be given in the correct order to produce sorted VCF on output unless
the -a, --allow-overlaps option is specified. With the --naive option, the files
are concatenated without being recompressed, which is very fast but dangerous
if the BCF headers differ.
"""

    additional_args = [
        ToolInput(
            "allowOverLaps",
            Boolean(optional=True),
            prefix="-a",
            doc="First coordinate of the next file can precede last record of the current file.",
        ),
        ToolInput(
            "compactPS",
            Boolean(optional=True),
            prefix="-c",
            doc="Do not output PS tag at each site, only at the start of a new phase set block.",
        ),
        ToolInput(
            "rmDups",
            String(optional=True),
            prefix="-d",
            doc="Output duplicate records present in multiple files only once: <snps|indels|both|all|none>",
        ),
        ToolInput(
            "rmDupsNone", Boolean(optional=True), prefix="-d", doc="Alias for -d none"
        ),
        ToolInput(
            "fileList",
            File(optional=True),
            prefix="-f",
            doc="Read the list of files from a file.",
        ),
        ToolInput(
            "ligate",
            Boolean(optional=True),
            prefix="-l",
            doc="Ligate phased VCFs by matching phase at overlapping haplotypes",
        ),
        ToolInput(
            "noVersion",
            Boolean(optional=True),
            prefix="--no-version",
            doc="Do not append version and command line information to the output VCF header.",
        ),
        ToolInput(
            "naive",
            Boolean(optional=True),
            prefix="-n",
            doc="Concatenate files without recompression (dangerous, use with caution)",
        ),
        ToolInput(
            "outputType",
            String(optional=True),
            prefix="-O",
            default="z",
            doc="--output-type b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), "
            "compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping "
            "between bcftools subcommands to speed up performance by removing "
            "unnecessary compression/decompression and VCF←→BCF conversion.",
        ),
        ToolInput(
            "minPG",
            Int(optional=True),
            prefix="-q",
            doc="Break phase set if phasing quality is lower than <int> [30]",
        ),
        ToolInput(
            "regions",
            String(optional=True),
            prefix="-r",
            doc="--regions chr|chr:pos|chr:from-to|chr:from-[,…]: Comma-separated list of regions, "
            "see also -R, --regions-file. Note that -r cannot be used in combination with -R.",
        ),
        ToolInput(
            "regionsFile",
            File(optional=True),
            prefix="-R",
            doc="--regions-file: Regions can be specified either on command line or in a VCF, BED, or "
            "tab-delimited file (the default). The columns of the tab-delimited file are: CHROM, POS, "
            "and, optionally, POS_TO, where positions are 1-based and inclusive. The columns of the "
            "tab-delimited BED file are also CHROM, POS and POS_TO (trailing columns are ignored), "
            "but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather "
            "than the 1-based tab-delimited file, the file must have the '.bed' or '.bed.gz' suffix "
            "(case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and "
            "tabix-indexed region files are streamed. Note that sequence names must match exactly, 'chr20'"
            " is not the same as '20'. Also note that chromosome ordering in FILE will be respected, "
            "the VCF will be processed in the order in which chromosomes first appear in FILE. "
            "However, within chromosomes, the VCF will always be processed in ascending genomic coordinate "
            "order no matter what order they appear in FILE. Note that overlapping regions in FILE can "
            "result in duplicated out of order positions in the output. This option requires indexed "
            "VCF/BCF files. Note that -R cannot be used in combination with -r.",
        ),
        ToolInput(
            "threads",
            Int(optional=True),
            prefix="--threads",
            doc="Number of output compression threads to use in addition to main thread. "
            "Only used when --output-type is b or z. Default: 0.",
        ),
    ]


class BcfToolsConcat_1_9(BcfTools_1_9, BcfToolsConcatBase):
    pass

BcfToolsConcatLatest = BcfToolsConcat_1_9


### BCF_TOOLS_SORT ###

CORES_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 1,
            CaptureType.EXOME: 1,
            CaptureType.THIRTYX: 1,
            CaptureType.NINETYX: 1,
            CaptureType.THREEHUNDREDX: 1,
        },
    )
]

MEM_TUPLE = [
    (
        CaptureType.key(),
        {
            CaptureType.CHROMOSOME: 16,
            CaptureType.EXOME: 16,
            CaptureType.THIRTYX: 16,
            CaptureType.NINETYX: 16,
            CaptureType.THREEHUNDREDX: 16,
        },
    )
]



class BcfToolsSortBase(BcfToolsToolBase, ABC):
    def tool(self) -> str:
        return "bcftoolssort"

    def friendly_name(self) -> str:
        return "BCFTools: Sort"

    def base_command(self):
        return ["bcftools", "sort"]

    def cpus(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, CORES_TUPLE)
        if val:
            return val
        return 1

    def memory(self, hints: Dict[str, Any]):
        val = get_value_for_hints_and_ordered_resource_tuple(hints, MEM_TUPLE)
        if val:
            return val
        return 8

    def inputs(self):
        return [
            ToolInput(
                "vcf",
                UnionType(Vcf, CompressedVcf),
                position=1,
                doc="The VCF file to sort",
            ),
            ToolInput(
                "outputFilename",
                Filename(suffix=".sorted", extension=".vcf.gz"),
                prefix="--output-file",
                doc="(-o) output file name [stdout]",
            ),
            ToolInput(
                "outputType",
                String(optional=True),
                prefix="--output-type",
                default="z",
                doc="(-O) b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]",
            ),
            ToolInput(
                "tempDir",
                String(optional=True),
                prefix="--temp-dir",
                doc="(-T) temporary files [/tmp/bcftools-sort.XXXXXX/]",
            ),
        ]

    def outputs(self):
        return [ToolOutput("out", CompressedVcf, glob=InputSelector("outputFilename"))]

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=datetime(2019, 5, 9),
            dateUpdated=datetime(2019, 7, 11),
            institution=None,
            doi=None,
            citation=None,
            keywords=["BCFTools", "sort"],
            documentationUrl="",
            documentation="""About:   Sort VCF/BCF file.
Usage:   bcftools sort [OPTIONS] <FILE.vcf>""",
        )

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "outputType": "z",
                    "vcf": f"{remote_dir}/NA12878-BRCA1.generated.gathered.vcf.gz",
                },
                output=CompressedVcf.basic_test(
                    "out",
                    11602,
                    221,
                    ["GATKCommandLine"],
                    "fcc35adbb0624abc91f6de2e9042f749",
                ),
            )
        ]

class BcfToolsSort_1_9(BcfTools_1_9, BcfToolsSortBase):
    pass


BcfToolsSortLatest = BcfToolsSort_1_9
