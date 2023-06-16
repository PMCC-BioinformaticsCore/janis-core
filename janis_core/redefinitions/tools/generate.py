
# grep ^@SQ reference.dict | cut -f2,3 | sed 's/SN://' | sed 's/LN://'
import os
from datetime import datetime
from typing import List, Dict, Any, Optional
from janis_core import TOutput, OutputDocumentation, File
from janis_core.tool.test_classes import TTestCase

from ..types import FastaDict, TextFile, Bed
from .bioinformaticstool import BioinformaticsPythonTool

import janis_core as j



class GenerateMantaConfig(BioinformaticsPythonTool):
    @staticmethod
    def code_block(output_filename: str = "output.txt") -> Dict[str, Any]:
        """
        :param output_filename: Filename to output to
        """
        with open(output_filename, "w+") as out:
            out.write("\n")
            out.write("#\n")
            out.write(
                "# This section contains all configuration settings for the top-level manta workflow,\n"
            )
            out.write("#\n")
            out.write("[manta]\n")
            out.write("\n")
            out.write(
                "referenceFasta = /illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa\n"
            )
            out.write("\n")
            out.write(
                "# Run discovery and candidate reporting for all SVs/indels at or above this size\n"
            )
            out.write(
                "# Separate option (to provide different default) used for runs in RNA-mode\n"
            )
            out.write("minCandidateVariantSize = 8\n")
            out.write("rnaMinCandidateVariantSize = 1000\n")
            out.write("\n")
            out.write(
                "# Remove all edges from the graph unless they're supported by this many 'observations'.\n"
            )
            out.write(
                "# Note that one supporting read pair or split read usually equals one observation, but evidence is sometimes downweighted.\n"
            )
            out.write("minEdgeObservations = 3\n")
            out.write("\n")
            out.write(
                "# If both nodes of an edge have an edge count higher than this, then skip evaluation of the edge.\n"
            )
            out.write("# Set to 0 to turn this filtration off\n")
            out.write("graphNodeMaxEdgeCount = 10\n")
            out.write("\n")
            out.write(
                "# Run discovery and candidate reporting for all SVs/indels with at least this\n"
            )
            out.write("# many spanning support observations\n")
            out.write("minCandidateSpanningCount = 3\n")
            out.write("\n")
            out.write(
                "# After candidate identification, only score and report SVs/indels at or above this size:\n"
            )
            out.write("minScoredVariantSize = 50\n")
            out.write("\n")
            out.write(
                '# minimum VCF "QUAL" score for a variant to be included in the diploid vcf:\n'
            )
            out.write("minDiploidVariantScore = 10\n")
            out.write("\n")
            out.write(
                '# VCF "QUAL" score below which a variant is marked as filtered in the diploid vcf:\n'
            )
            out.write("minPassDiploidVariantScore = 20\n")
            out.write("\n")
            out.write(
                "# minimum genotype quality score below which single samples are filtered for a variant in the diploid vcf:\n"
            )
            out.write("minPassDiploidGTScore = 15\n")
            out.write("\n")
            out.write(
                "# somatic quality scores below this level are not included in the somatic vcf:\n"
            )
            out.write("minSomaticScore = 10\n")
            out.write("\n")
            out.write(
                "# somatic quality scores below this level are filtered in the somatic vcf:\n"
            )
            out.write("minPassSomaticScore = 30\n")
            out.write("\n")
            out.write(
                "# Remote read retrieval is used ot improve the assembly of putative insertions by retrieving any mate reads in remote\n"
            )
            out.write(
                "# locations with poor mapping quality, which pair to confidently mapping reads near the insertion locus. These reads\n"
            )
            out.write(
                "# can help to fully assemble longer insertions, under certain circumstances this feature can add a very large runtime\n"
            )
            out.write(
                "# burden. For instance, given the very high chimeric pair rates found in degraded FFPE samples, the runtime of the read\n"
            )
            out.write(
                "# retrieval process can be unpredicable. For this reason the feature is disabled by default for somatic variant calling.\n"
            )
            out.write(
                "# This feature can be enabled/disabled separately for germline and cancer calling below.\n"
            )
            out.write("#\n")
            out.write(
                '# Here "CancerCallingModes" includes tumor-normal subtraction and tumor-only calling. "GermlineCallingModes" includes\n'
            )
            out.write("# all other calling modes.\n")
            out.write("# custom set-up: https://github.com/Illumina/manta/issues/213\n")
            out.write(
                "enableRemoteReadRetrievalForInsertionsInGermlineCallingModes = 0\n"
            )
            out.write(
                "enableRemoteReadRetrievalForInsertionsInCancerCallingModes = 0\n"
            )
            out.write("\n")
            out.write(
                "# Set if an overlapping read pair will be considered as evidence\n"
            )
            out.write("# Set to 0 to skip overlapping read pairs\n")
            out.write("useOverlapPairEvidence = 0\n")
            out.write("\n")
            return {"out": output_filename}

    def outputs(self) -> List[TOutput]:
        return [
            TOutput(
                "out",
                File,
                doc=OutputDocumentation(doc="Custom Manta config file"),
            )
        ]

    def id(self) -> str:
        return "GenerateMantaConfig"

    def friendly_name(self) -> str:
        return "GenerateMantaConfig"

    def tool_provider(self):
        return "Peter MacCallum Cancer Centre"

    def version(self):
        return "v0.1.0"

    def bind_metadata(self):
        self.metadata.dateCreated = datetime(2021, 5, 27)
        self.metadata.dateUpdated = datetime(2021, 5, 27)
        self.metadata.contributors = ["Jiaan Yu"]
        self.metadata.documentation = """\
Generate custom manta config file.       
        """





class GenerateIntervalsByChromosome(BioinformaticsPythonTool):
    @staticmethod
    def code_block(
        reference: FastaDict,
        prefix="chr",
        allowed_contigs: Optional[List[str]] = None,
        max_size: Optional[int] = None,
        overlap: Optional[int] = 0,
        single_file=False,
    ) -> Dict[str, Any]:
        """
        This tool generates a list of BED intervals based on the reference input.
        It supports breaking up the intervals into subregions of length 'max_size',
        with an overlap (to consider indels spanning regions). Note the max_size is
        EXCLUSIVE of overlap, so the actual interval_size is (max_size + overlap).
        :param reference: FASTA reference with ^.dict reference
        :param prefix: contig prefix, default 'chr'
        :param allowed_contigs: Limits allowed_contigs to a list, this defaults of Human CHRs, 1-23,X,Y,Z
        :param max_size: Max size of interval, maybe 5000 for VarDict.
        :param overlap: Consider indels spanning regions, so choose
        :param single_file: Produce a SINGLE .bed file with all the listed regions
        """
        from re import sub

        if max_size is not None and overlap >= max_size:
            raise Exception(
                f"max_size ({max_size}) must be greater than overlap ({overlap})"
            )

        # Allowed contigs: use the standard human genome if none are provided
        # include M / MT for hg19 / hg39
        if allowed_contigs is None:
            allowed_contigs = list(
                map(lambda el: f"{prefix}{el}", [*range(23), "X", "Y", "M", "MT"])
            ) + [*range(23), "X", "Y", "M", "MT"]
        allowed_contigs = set(allowed_contigs)

        def contig_label(contig: str) -> str:
            """
            Add the prefix if the contig doesn't already contain it
            """
            if prefix and not contig.startswith(prefix):
                return prefix + contig
            return contig

        def get_contig_and_length_from_line(line):
            """
            Iterate through each col in TSV line, and find SN / LN
            :returns CONTIG, LENGTH (ordered tuple)
            """
            contig, length = None, None
            for col in line.split("\t"):
                if col.startswith("SN:") and len(col) > 3:
                    contig = col.strip("SN:")
                elif col.startswith("LN:") and len(col) > 3:
                    length = col.strip("LN:")
                if contig is not None and length is not None:
                    break

            return contig, length

        def prepare_regions(contig, length) -> List[List[str]]:
            """
            Split the region into INTERVALS for (max_size + overlap) if REQUIRED,
            else return a
            :param contig:
            :param length:
            :return:
            """
            length = int(length)
            label = contig_label(contig)

            # BASE case, the interval fits within the max_size, just return a single row
            if max_size is None or length < max_size:
                return [[contig, "1", str(length), str(label)]]

            # ELSE split into subregions
            subregions = []
            # BED regions start at 1
            start, counter, finish = 1, 1, None
            while finish is None or finish < length:
                finish = min((finish or 0) + max_size, length)
                subregions.append(
                    [str(contig), str(start), str(finish), f"{label}_{counter}"]
                )
                start = finish - overlap + 1
                counter += 1

            return subregions

        # Get the ^.dict from the end of the .fasta | .fa filename
        ref_dict = sub("\.fa(sta)?$", ".dict", reference)

        regions = []
        all_prepped_regions = []

        with open(ref_dict) as ref:
            for line in ref:
                if "GL" in line or "SN" not in line:
                    continue

                contig, length = get_contig_and_length_from_line(line)
                # If we couldn't find a contig or sequence length, skip
                if contig is None or length is None:
                    continue
                # Only proceed if the contig OR (prefix + contig) in allowed_contigs
                if not (
                    contig in allowed_contigs or prefix + contig in allowed_contigs
                ):
                    continue

                # Get total region list for contig, and considering max_size
                regions_to_write = prepare_regions(contig, length)
                prepped = [("\t".join(region) + "\n") for region in regions_to_write]
                all_prepped_regions.extend(prepped)

                if not single_file:
                    regions.append(f"{contig_label(contig)}.bed")
                    with open(regions[-1], "w") as f:
                        f.writelines(prepped)

        if single_file:
            regions.append("regions.bed")
            with open("regions.bed", "w") as f:
                f.writelines(all_prepped_regions)

        return {"out_regions": regions}

    def outputs(self) -> List[j.TOutput]:
        return [j.TOutput("out_regions", j.Array(Bed))]

    def id(self) -> str:
        return "GenerateIntervalsByChromosome"

    def friendly_name(self) -> Optional[str]:
        return "Generating genomic intervals by chromosome"

    def version(self):
        return "v0.1.0"

    def bind_metadata(self):
        meta: j.ToolMetadata = self.metadata

        meta.contributors = ["Michael Franklin"]
        meta.dateCreated = datetime(2020, 10, 19)
        meta.dateUpdated = datetime(2020, 10, 19)

    def tool_provider(self):
        return "Peter MacCallum Cancer Centre"







class GenerateGenomeFileForBedtoolsCoverage(BioinformaticsPythonTool):
    @staticmethod
    def code_block(
        reference: FastaDict, output_filename: str = "genome_file.txt"
    ) -> Dict[str, Any]:
        """
        :param reference: Reference file to generate genome for (must have ^.dict) pattern
        :param output_filename: Filename to output to
        """
        from re import sub

        ref_dict = sub("\.fa(sta)?$", ".dict", reference)

        with open(ref_dict) as inp, open(output_filename, "w+") as out:
            for line in inp:
                if not line.startswith("@SQ"):
                    continue
                pieces = line.split("\t")
                chrom = pieces[1].replace("SN:", "")
                length = pieces[2].replace("LN:", "")

                out.write(f"{chrom}\t{length}\n")

            return {"out": output_filename}

    def outputs(self) -> List[TOutput]:
        return [
            TOutput(
                "out",
                TextFile,
                doc=OutputDocumentation(doc="Genome file for BedToolsCoverage"),
            )
        ]

    def id(self) -> str:
        return "GenerateGenomeFileForBedtoolsCoverage"

    def friendly_name(self) -> str:
        return "Generate genome for BedtoolsCoverage"

    def tool_provider(self):
        return "Peter MacCallum Cancer Centre"

    def version(self):
        return "v0.1.0"

    def bind_metadata(self):
        self.metadata.dateCreated = datetime(2020, 7, 21)
        self.metadata.dateUpdated = datetime(2020, 6, 2)
        self.metadata.contributors = ["Michael Franklin", "Jiaan Yu"]
        self.metadata.documentation = """\
Generate --genome FILE for BedToolsCoverage      
        """

    def tests(self):
        remote_dir = "https://swift.rc.nectar.org.au/v1/AUTH_4df6e734a509497692be237549bbe9af/janis-test-data/bioinformatics/wgsgermline_data"
        return [
            TTestCase(
                name="basic",
                input={
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                },
                output=TextFile.basic_test("out", 15, "chr17\t83257441\n", 1),
            ),
            TTestCase(
                name="minimal",
                input={
                    "reference": f"{remote_dir}/Homo_sapiens_assembly38.chr17.fasta",
                },
                output=self.minimal_test(),
            ),
        ]
