

"""
Each modification of this tool should duplicate this code
"""
from datetime import datetime
from typing import List, Dict, Any, Optional

from janis_core import PythonTool, File, Array, ToolMetadata
from janis_core.tool.tool import TOutput
from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput,
    TTestPreprocessor,
)

from .bioinformaticstool import BioinformaticsPythonTool


class ParseFastqcAdapters(BioinformaticsPythonTool):
    @staticmethod
    def code_block(
        read1_fastqc_datafile: File,
        read2_fastqc_datafile: File,
        adapters_lookup: File,
        contamination_lookup: File,
    ):
        """
        :param read1_fastqc_datafile: fastqc_datafile of read 1

        :param read2_fastqc_datafile: fastqc_datafile of read 2

        :param adapters_lookup: Specifies a file which contains the list of adapter sequences which will
            be explicity searched against the library. The file must contain sets of named adapters in
            the form name[tab]sequence. Lines prefixed with a hash will be ignored.

        :param contamination_lookup: Specifies a file which contains the list of universal adapter
            sequences which will be explicity searched against the library.

        :return: sequences in list
        """

        import mmap, re, csv
        from io import StringIO
        from sys import stderr

        def get_overrepresented_text(f):
            """
            Get the table "Overrepresented sequences" within the fastqc_data.txt
            """
            adapt_section_query = (
                br"(?s)>>Overrepresented sequences\t\S+\n(.*?)>>END_MODULE"
            )
            adapt_section_fail_query = (
                br"(?s)>>Overrepresented sequences\tfail+\n(.*?)>>END_MODULE"
            )
            # fastqc_datafile could be fairly large, so we'll use mmap, and then
            with open(f) as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as fp:
                overrepresented_sequences_match = re.search(adapt_section_query, fp)
                if overrepresented_sequences_match is None:
                    raise Exception(
                        f"Couldn't find query ('{adapt_section_query.decode('utf8')}') in {f}"
                    )
                elif re.search(adapt_section_fail_query, fp):
                    # parse the sequences when fail
                    overrepresented_sequences_match = (
                        overrepresented_sequences_match.groups()[0].decode("utf8")
                    )
                else:
                    overrepresented_sequences_match = None
                return overrepresented_sequences_match

        def parse_tsv_table(tbl: str, skip_headers):
            """
            Parse a TSV table from a string using csvreader
            """

            rd = csv.reader(StringIO(tbl), delimiter="\t", quotechar='"')
            ret = list(rd)
            if len(ret) == 0:
                return ret
            if skip_headers:
                ret.pop(0)  # discard headers
            return ret

        def get_adapt_map(adapter_file):
            """
            Helper method to parse the file 'adapters_lookup'/ 'contamination_lookup' with
            format: 'name[tab]sequence' into the dictionary: '{ name: sequence }'
            """
            adapter_map = {}
            with open(adapter_file) as fp:
                for row in fp:
                    st = row.strip()
                    if not st or st.startswith("#"):
                        continue

                    # In reality, the format is [\t+] (more than one tab)
                    # so we'll just split on a tab, and remove all the empty elements.
                    split = [f for f in st.split("\t") if bool(f) and len(f) > 0]

                    # Invalid format for line, so skip it.
                    if len(split) != 2:
                        print(
                            f"Skipping adapter line '{st}' as irregular elements ({len(split)})",
                            file=stderr,
                        )
                        continue
                    adapter_map[split[0]] = split[1]
            return adapter_map

        # Start doing the work
        # Look up overrepresented sequences
        i = 0
        for qc_file in [read1_fastqc_datafile, read2_fastqc_datafile]:
            overrepresented_ids = set()
            overrepresented_sequences = []
            text = get_overrepresented_text(qc_file)
            if text:
                overrepresented_ids = overrepresented_ids.union(
                    set(
                        a[3].split("(")[0].strip(" ")
                        for a in parse_tsv_table(text, skip_headers=True)
                    )
                )
                # print(overrepresented_ids)

            if overrepresented_ids:
                adapter_map = get_adapt_map(contamination_lookup)
                # print(overrepresented_ids)
                # print(adapter_map)
                for oid in overrepresented_ids:
                    if oid in adapter_map:
                        print(
                            f"Identified sequence '{oid}' as '{adapter_map.get(oid)}' in lookup",
                            file=stderr,
                        )
                        overrepresented_sequences.append(adapter_map.get(oid))
                    else:
                        print(
                            f"Couldn't find a corresponding sequence for '{oid}' in lookup map",
                            file=stderr,
                        )

            # Look up common adapter sequences
            adapter_sequences = []
            adapter_status = False
            for line in open(qc_file, "r"):
                if line.startswith(">>Adapter Content"):
                    adapter_qc_line = line
                    adapter_status = True
                elif adapter_status == False:
                    continue
                elif adapter_status == True and not line.startswith(">>END_MODULE"):
                    if line.startswith("#Position"):
                        adapter_names = line.strip("\n").split("\t")
                    else:
                        adapter_vals = line.strip("\n").split("\t")

            # Parse adapter id if adapter qc fails
            if "fail" in adapter_qc_line:
                adapter_map = get_adapt_map(adapters_lookup)
                for (aid, percentage) in zip(adapter_names[1:], adapter_vals[1:]):
                    if float(percentage) >= 10:
                        if aid in adapter_map:
                            sequence = adapter_map.get(aid)
                            print(
                                f"Identified adapter '{aid}' sequence '{sequence}' in lookup",
                                file=stderr,
                            )
                            adapter_sequences.append(sequence)
                        else:
                            print(
                                f"Couldn't find a corresponding sequence for '{aid}' in lookup map",
                                file=stderr,
                            )
            else:
                pass

            if i == 0:
                read1_sequences = list(overrepresented_sequences) + adapter_sequences
            else:
                read2_sequences = list(overrepresented_sequences) + adapter_sequences
            i += 1

        return {
            "out_R1_sequences": read1_sequences,
            "out_R2_sequences": read2_sequences,
        }

    def outputs(self) -> List[TOutput]:
        return [
            TOutput("out_R1_sequences", Array(str)),
            TOutput("out_R2_sequences", Array(str)),
        ]

    def id(self) -> str:
        return "ParseFastqcAdapters"

    def friendly_name(self):
        return "Parse FastQC Adapters"

    def version(self):
        return "v0.2.0"

    def tool_provider(self):
        return "Peter MacCallum Cancer Centre"

    def bind_metadata(self):
        self.metadata.documentation = (
            "Parse overrepresented region and lookup in adapter table"
        )
        self.metadata.contributors = ["Michael Franklin", "Jiaan Yu"]
        self.metadata.dateCreated = datetime(2020, 1, 7)
        self.metadata.dateUpdated = datetime(2021, 10, 6)
        self.metadata.version = "0.2.0"

    def tests(self):
        # New files needed to update to the space
        # Test not founding anything useful
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "read1_fastqc_datafile": f"{remote_dir}/NA12878-BRCA1_R1.fastqc_data.txt",
                    "read2_fastqc_datafile": f"{remote_dir}/NA12878-BRCA1_R2.fastqc_data.txt",
                    "adapters_lookup": f"{remote_dir}/adapter_list.txt",
                    "contamination_lookup": f"{remote_dir}/contaminant_list.txt",
                },
                output=[],
            ),
        ]
