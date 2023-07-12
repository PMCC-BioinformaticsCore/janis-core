

import operator
import os
from datetime import datetime

from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput,
    TTestPreprocessor,
)

from . import UnixTool
from janis_core import (
    ToolOutput,
    ToolInput,
    Boolean,
    String,
    File,
    ToolMetadata,
    InputSelector,
)
from ..types import Gunzipped

class UncompressArchive(UnixTool):
    def tool(self):
        return "UncompressArchive"

    def friendly_name(self):
        return "UncompressArchive"

    def tool_provider(self):
        return "GNU Project"

    def base_command(self):
        return "gunzip"

    def inputs(self):
        return [
            ToolInput("file", Gunzipped(), position=1, localise_file=True),
            *self.additional_inputs,
        ]

    def outputs(self):
        return [
            ToolOutput(
                "out", File, glob=InputSelector("file", remove_file_extension=True)
            )
        ]

    additional_inputs = [
        ToolInput(
            "stdout",
            Boolean(optional=True),
            prefix="-c",
            doc="write on standard output, keep original files unchanged",
        ),
        ToolInput(
            "decompress",
            Boolean(optional=True),
            prefix="-d",
            default=True,
            doc="decompress",
        ),
        ToolInput(
            "force",
            Boolean(optional=True),
            prefix="-f",
            doc="force overwrite of output file and compress links",
        ),
        ToolInput(
            "keep",
            Boolean(optional=True),
            prefix="-k",
            doc="keep (don't delete) input files",
        ),
        ToolInput(
            "list",
            Boolean(optional=True),
            prefix="-l",
            doc="list compressed file contents",
        ),
        ToolInput(
            "noName",
            Boolean(optional=True),
            prefix="-n",
            doc="do not save or restore the original name and time stamp",
        ),
        ToolInput(
            "name",
            Boolean(optional=True),
            prefix="-N",
            doc="save or restore the original name and time stamp",
        ),
        ToolInput(
            "quiet", Boolean(optional=True), prefix="-q", doc="suppress all warnings"
        ),
        ToolInput(
            "recursive",
            Boolean(optional=True),
            prefix="-r",
            doc="operate recursively on directories",
        ),
        ToolInput(
            "suffix",
            String(optional=True),
            prefix="-s",
            doc="use suffix SUF on compressed files",
        ),
        ToolInput(
            "test",
            Boolean(optional=True),
            prefix="-t",
            doc="test compressed file integrity",
        ),
        ToolInput("fast", Boolean(optional=True), prefix="-1", doc="compress faster"),
        ToolInput("best", Boolean(optional=True), prefix="-9", doc="compress better"),
        ToolInput(
            "rsyncable",
            Boolean(optional=True),
            prefix="--rsyncable",
            doc="Make rsync-friendly archive",
        ),
    ]

    def bind_metadata(self):
        return ToolMetadata(
            contributors=["Jiaan Yu"],
            dateCreated=datetime(2020, 6, 11),
            dateUpdated=datetime(2020, 6, 11),
            documentation="",
        )

    def tests(self):
        return [
            TTestCase(
                name="basic",
                input={
                    "file": "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/petermac_testdata/1000G_phase1.snps.high_confidence.hg38.BRCA1.vcf.gz",
                },
                output=[
                    TTestExpectedOutput(
                        tag="out",
                        preprocessor=TTestPreprocessor.FileSize,
                        operator=operator.eq,
                        expected_value=160525,
                    ),
                    TTestExpectedOutput(
                        tag="out",
                        preprocessor=TTestPreprocessor.LineCount,
                        operator=operator.eq,
                        expected_value=625,
                    ),
                ],
            )
        ]
