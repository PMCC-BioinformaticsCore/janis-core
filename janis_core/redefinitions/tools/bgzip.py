
from typing import List
from abc import abstractmethod, ABC

from janis_core import ToolMetadata
from janis_core import (
    ToolOutput,
    ToolInput,
    Boolean,
    Int,
    File,
    InputSelector,
    Filename,
    ToolArgument,
)
from janis_core.tool.test_classes import TTestCase
from .bioinformaticstool import BioinformaticsTool
from ..types import Gunzipped, CompressedVcf



class HtsLibBase(BioinformaticsTool, ABC):
    def tool_provider(self):
        return "HTSLib"
    
class HTSLib_1_2_1(HtsLibBase):
    def container(self):
        return "biodckrdev/htslib:1.2.1"

    def version(self):
        return "1.2.1"
    
HTSLibLatest = HTSLib_1_2_1

class BGZipBase(HtsLibBase, ABC):
    def tool_provider(self):
        return "htslib"

    def tool(self):
        return "bgzip"

    def base_command(self):
        return "bgzip"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("file", File(), position=100, doc="File to bgzip compress"),
            ToolInput(
                "outputFilename",
                Filename(
                    prefix=InputSelector("file").basename(),
                    extension=".gz",
                ),
                position=102,
            ),
            *self.additional_args,
        ]

    def arguments(self):
        return [ToolArgument(">", position=101, shell_quote=False)]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput(
                "out",
                Gunzipped(),
                glob=InputSelector("outputFilename"),
            )
        ]

    def friendly_name(self):
        return "BGZip"

    def bind_metadata(self):
        from datetime import date

        return ToolMetadata(
            contributors=["Michael Franklin"],
            dateCreated=date(2018, 12, 24),
            dateUpdated=date(2019, 1, 24),
            institution="HTSLib",
            doi=None,
            citation=None,  # "",
            keywords=["htslib", "bgzip", "compression"],
            documentationUrl="http://www.htslib.org/doc/bgzip.html",
            documentation="""bgzip – Block compression/decompression utility

Bgzip compresses files in a similar manner to, and compatible with, gzip(1). The file is compressed 
into a series of small (less than 64K) 'BGZF' blocks. This allows indexes to be built against the 
compressed file and used to retrieve portions of the data without having to decompress the entire file.

If no files are specified on the command line, bgzip will compress (or decompress if the -d option is used) 
standard input to standard output. If a file is specified, it will be compressed (or decompressed with -d). 
If the -c option is used, the result will be written to standard output, otherwise when compressing bgzip 
will write to a new file with a .gz suffix and remove the original. When decompressing the input file must 
have a .gz suffix, which will be removed to make the output name. 
Again after decompression completes the input file will be removed.""".strip(),
        )

    @abstractmethod
    def container(self):
        raise Exception(
            "An error likely occurred when resolving the method order for docker for the tabix classes "
            "or you're trying to execute the docker method of the base class (ie, don't do that). "
            "The method order resolution must preference BwaBase subclasses, "
            "and the subclass must contain a definition for docker."
        )

    additional_args = [
        ToolInput(
            "offset",
            Int(optional=True),
            prefix="--offset",
            doc="b: Decompress to standard output from virtual file position "
            "(0-based uncompressed offset). Implies -c and -d.",
        ),
        ToolInput(
            "stdout",
            Boolean(optional=True),
            default=True,
            prefix="--stdout",
            doc="c: Write to standard output, keep original files unchanged.",
        ),
        ToolInput(
            "decompress",
            Boolean(optional=True),
            prefix="--decompress",
            doc="d: Decompress.",
        ),
        ToolInput(
            "force",
            Boolean(optional=True),
            prefix="--force",
            doc="f: Overwrite files without asking.",
        ),
        ToolInput(
            "help",
            Boolean(optional=True),
            prefix="--help",
            doc="h: Displays a help message.",
        ),
        ToolInput(
            "index",
            Boolean(optional=True),
            prefix="--index",
            doc="i: Create a BGZF index while compressing. Unless the -I option is used, "
            "this will have the name of the compressed file with .gzi appended to it.",
        ),
        ToolInput(
            "indexName",
            File(optional=True),
            prefix="--index-name",
            doc="-I: Index file name.",
        ),
        ToolInput(
            "compress",
            Int(optional=True),
            prefix="--compress",
            doc="l: Compression level to use when compressing. From 0 to 9, or -1 "
            "for the default level set by the compression library. [-1]",
        ),
        ToolInput(
            "reindex",
            Boolean(optional=True),
            prefix="--reindex",
            doc="r: Rebuild the index on an existing compressed file.",
        ),
        ToolInput(
            "rebgzip",
            Boolean(optional=True),
            prefix="--rebgzip",
            doc="g: Try to use an existing index to create a compressed file with matching block offsets. "
            "Note that this assumes that the same compression library and level are in use "
            "as when making the original file. Don't use it unless you know what you're doing.",
        ),
        ToolInput(
            "size",
            Int(optional=True),
            prefix="--size",
            doc="s: Decompress INT bytes (uncompressed size) to standard output. Implies -c.",
        ),
        ToolInput(
            "threads",
            Int(optional=True),
            prefix="--threads",
            doc="@: Number of threads to use [1].",
        ),
    ]

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "file": f"{remote_dir}/NA12878-BRCA1.generated.gathered.vcf",
                },
                output=CompressedVcf.basic_test(
                    "out",
                    11500,
                    221,
                    ["GATKCommandLine"],
                    "b7acb0a9900713cc7da7aeed5160c971",
                ),
            )
        ]


class BGZipLatest(HTSLibLatest, BGZipBase):
    pass