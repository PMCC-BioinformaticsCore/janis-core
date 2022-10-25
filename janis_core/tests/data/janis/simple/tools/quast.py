# NOTE
# This is an automated translation of the 'quast' version '5.0.2' tool from a Galaxy XML tool wrapper.  
# Translation was performed by the gxtool2janis program (in-development)

import sys
sys.path.append('/home/grace/work/pp/gxtool2janis')

from janis_core.types.common_data_types import Boolean
from janis_core.types.common_data_types import String
from janis_core.types.common_data_types import Float
from janis_core.types.common_data_types import File
from janis_core.types.common_data_types import Int
from janis_core import CommandToolBuilder
from janis_core import WildcardSelector
from janis_core import ToolMetadata
from janis_core import ToolOutput
from janis_core import ToolInput
from janis_core import Array


metadata = ToolMetadata(
    short_documentation="Genome assembly Quality",
    keywords=[],
    contributors=['gxtool2janis'],
    dateCreated="2022-08-29 10:49:40",
    dateUpdated="2022-08-29 10:49:40",
    version="5.0.2",
    doi="None",
    citation="tool xml missing citation",
    documentationUrl=None,
    documentation="""

**What it does**

QUAST = QUality ASsessment Tool. The tool evaluates genome assemblies by computing various metrics.

If you have one or multiple genome assemblies, you can assess their quality with Quast. It works with or without reference genome. If you are new to Quast, start by reading its `manual page <http://quast.bioinf.spbau.ru/manual.html>`_.

**Using Quast without reference**

Without reference Quast can calculate a number of assembly related-metrics but cannot provide any information about potential misassemblies, inversions, translocations, etc. Suppose you have three assemblies produced by Unicycler corresponding to three different antibiotic treatments *car*, *pit*, and *cef* (these stand for carbenicillin, piperacillin, and cefsulodin, respectively). Evaluating them without reference will produce the following Quast outputs:
 
 * Quast report in HTML format
 * `Contig viewer <http://quast.bioinf.spbau.ru/manual.html#sec3.4>`_ (an HTML file)
 * `Quast report <http://quast.bioinf.spbau.ru/manual.html#sec3.1.1>`_ in Tab-delimited format
 * Quast log (a file technical information about Quast tool execution)

The **tab delimited Quast report** will contain the following information::

 Assembly                  pit_fna cef_fna car_fna
 # contigs (>= 0 bp)           100      91      94
 # contigs (>= 1000 bp)         62      58      61
 Total length (>= 0 bp)    6480635 6481216 6480271
 Total length (>= 1000 bp) 6466917 6468946 6467103
 # contigs                      71      66      70
 Largest contig             848753  848766  662053
 Total length              6473173 6474698 6473810
 GC (%)                      66.33   66.33   66.33
 N50                        270269  289027  254671
 N75                        136321  136321  146521
 L50                             7       7       8
 L75                            15      15      16
 # N's per 100 kbp            0.00    0.00    0.00

where values are defined as specified in `Quast manual <http://quast.bioinf.spbau.ru/manual.html#sec3.1.1>`_

**Quast report in HTML format** contains graphs in addition to the above metrics, while **Contig viewer** draws contigs ordered from longest to shortest. This ordering is suitable for comparing only largest contigs or number of contigs longer than a specific threshold. The viewer shows N50 and N75 with color and textual indication. If the reference genome is available or at least approximate genome length is known (see `--est-ref-size`), NG50 and NG75 are also shown. You can also tone down contigs shorter than a specified threshold using Icarus control panel:

.. image:: $PATH_TO_IMAGES/contig_view_noR.png
   :width: 558
   :height: 412

Also see `Plot description <http://quast.bioinf.spbau.ru/manual.html#sec3.2>`_ section of the manual. 

**Using Quast with reference**

Car, pit, and cef are in fact assemblies of *Pseudomonas aeruginosa* UCBPP-PA14, so we can use its genome as a reference (by supplying a Fasta file containing *P. aeruginosa* pa14 genome to **Reference genome** input box). The following outputs will be produced (note the alignment viewer):

 * Quast report in HTML format
 * `Contig viewer <http://quast.bioinf.spbau.ru/manual.html#sec3.4>`_ (an HTML file)
 * `Alignment viewer <http://quast.bioinf.spbau.ru/manual.html#sec3.4>`_ (an HTML file)
 * `Quast report <http://quast.bioinf.spbau.ru/manual.html#sec3.1.1>`_ in Tab-delimited format
 * Summary of `misassemblies <http://quast.bioinf.spbau.ru/manual.html#sec3.1.2>`_
 * Summary of `unaligned contigs <http://quast.bioinf.spbau.ru/manual.html#sec3.1.3>`_
 * Quast log (a file technical information about Quast tool execution)

With the reference Quast produces a much more comprehensive set of results::

 Assembly                  pit_fna cef_fna car_fna
 # contigs (>= 0 bp)           100      91      94
 # contigs (>= 1000 bp)         62      58      61
 Total length (>= 0 bp)    6480635 6481216 6480271
 Total length (>= 1000 bp) 6466917 6468946 6467103
 # contigs                      71      66      70
 Largest contig             848753  848766  662053
 Total length              6473173 6474698 6473810
 Reference length          6537648 6537648 6537648
 GC (%)                      66.33   66.33   66.33
 Reference GC (%)            66.29   66.29   66.29
 N50                        270269  289027  254671
 NG50                       270269  289027  254671
 N75                        136321  136321  146521
 NG75                       136321  136321  136321
 L50                             7       7       8
 LG50                            7       7       8
 L75                            15      15      16
 LG75                           15      15      17
 # misassemblies                 0       0       0
 # misassembled contigs          0       0       0
 Misassembled contigs length     0       0       0
 # local misassemblies           1       1       2
 # unaligned mis. contigs        0       0       0
 # unaligned contigs         0 + 0   0 + 0   0 + 0
                              part    part    part
 Unaligned length                0       0       0
 Genome fraction (%)        99.015  99.038  99.025
 Duplication ratio           1.000   1.000   1.000
 # N's per 100 kbp            0.00    0.00    0.00
 # mismatches per 100 kbp     3.82    3.63    3.49
 # indels per 100 kbp         1.19    1.13    1.13
 Largest alignment          848753  848766  662053
 Total aligned length      6473163 6474660 6473792
 NA50                       270269  289027  254671
 NGA50                      270269  289027  254671
 NA75                       136321  136321  146521
 NGA75                      136321  136321  136321
 LA50                            7       7       8
 LGA50                           7       7       8
 LA75                           15      15      16
 LGA75                          15      15      17 

where, again, values are defined as specified in `Quast manual <http://quast.bioinf.spbau.ru/manual.html#sec3.1.1>`_. You can see that this report includes a variety of data that can only be computer against a reference assembly. 

 Using reference also produces an **Alignment viewer**:

.. image:: $PATH_TO_IMAGES/Align_view.png
   :width: 515
   :height: 395

Alignment viewer highlights regions of interest as, in this case, missassemblies that can potentially point to genome rearrangements (see more `here <http://quast.bioinf.spbau.ru/manual.html#sec3.4>`_).

    
    """
)

inputs = [
	# Positionals
    ToolInput(
		'unknown1',
		File,
		doc="Unknown input. Accepts w.fastqc1.outTextFile. Please address.",
	),
    

	# Flags
	ToolInput(
		'circos',
		Boolean(optional=True),
		prefix='--circos',
		position=1,
		default="False",
	),
	ToolInput(
		'conservedGenesFinding',
		Boolean(optional=True),
		prefix='--conserved-genes-finding',
		position=1,
		default="False",
	),
	ToolInput(
		'eukaryote',
		Boolean(optional=True),
		prefix='--eukaryote',
		position=1,
		default="False",
	),
	ToolInput(
		'fragmented',
		Boolean(optional=True),
		prefix='--fragmented',
		position=1,
		default="False",
	),
	ToolInput(
		'fungus',
		Boolean(optional=True),
		prefix='--fungus',
		position=1,
		default="False",
	),
	ToolInput(
		'geneFinding',
		Boolean(optional=True),
		prefix='--gene-finding',
		position=1,
		default="False",
	),
	ToolInput(
		'glimmer',
		Boolean(optional=True),
		prefix='--glimmer',
		position=1,
		default="False",
	),
	ToolInput(
		'kMerStats',
		Boolean(optional=True),
		prefix='--k-mer-stats',
		position=1,
		default="False",
	),
	ToolInput(
		'large',
		Boolean(optional=True),
		prefix='--large',
		position=1,
		default="False",
	),
	ToolInput(
		'mgm',
		Boolean(optional=True),
		prefix='--mgm',
		position=1,
		default="False",
	),
	ToolInput(
		'rnaFinding',
		Boolean(optional=True),
		prefix='--rna-finding',
		position=1,
		default="False",
	),
	ToolInput(
		'skipUnalignedMisContigs',
		Boolean(optional=True),
		prefix='--skip-unaligned-mis-contigs',
		position=1,
		default="False",
	),
	ToolInput(
		'splitScaffolds',
		Boolean(optional=True),
		prefix='--split-scaffolds',
		position=1,
		default="False",
	),
	ToolInput(
		'strictNa',
		Boolean(optional=True),
		prefix='--strict-NA',
		position=1,
		default="False",
	),
	ToolInput(
		'testNoRef',
		Boolean(optional=True),
		prefix='--test-no-ref',
		position=1,
		default="False",
	),
	ToolInput(
		'useAllAlignments',
		Boolean(optional=True),
		prefix='--use-all-alignments',
		position=1,
		default="False",
	),

	# Options
	ToolInput(
		'ambiguityUsage',
		String,
		prefix='--ambiguity-usage',
		position=1,
		default="one",
		doc="How processing equally good alignments of a contig (probably repeats)?. possible values: all, none, one",
	),
	ToolInput(
		'contigThresholds',
		String(optional=True),
		prefix='--contig-thresholds',
		position=1,
		default="0,1000",
		doc="Comma-separated list of contig length thresholds (in bp). Used in # contigs ≥ x and total length (≥ x) metrics",
	),
	ToolInput(
		'estRefSize',
		Int(optional=True),
		prefix='--est-ref-size',
		position=1,
		default=None,
		doc="Estimated reference genome size (in bp) for computing NGx statistics",
	),
	ToolInput(
		'extensiveMisSize',
		Int,
		prefix='--extensive-mis-size',
		position=1,
		default=1000,
		doc="Lower threshold for the relocation size (gap or overlap size between left and right flanking sequence). Shorter relocations are considered as local misassemblies. It does not affect other types of extensive misassemblies (inversions and translocations). The default value is 1000 bp. Note that the threshold should be greater than maximum indel length which is 85 bp.",
	),
	ToolInput(
		'features',
		File(optional=True),
		prefix='--features',
		position=1,
		default=None,
		doc="Genomic feature positions in the reference genome. Gene coordinates for the reference genome",
	),
	ToolInput(
		'geneThresholds',
		String,
		prefix='--gene-thresholds',
		position=1,
		default="0,300,1500,3000",
		doc="Comma-separated list of thresholds (in bp) for gene lengths to find with a finding tool",
	),
	ToolInput(
		'kMerSize',
		Int,
		prefix='--k-mer-size',
		position=1,
		default=101,
		doc="Size of k",
	),
	ToolInput(
		'labels',
		String(optional=True),
		prefix='--labels',
		position=1,
		default=None,
	),
	ToolInput(
		'maxRefNum',
		Int,
		prefix='--max-ref-num',
		position=1,
		default=50,
		doc="Maximum number of reference genomes (per each assembly) to download after searching in the SILVA databa",
	),
	ToolInput(
		'minAlignment',
		Int,
		prefix='--min-alignment',
		position=1,
		default=65,
		doc="Minimum length of alignment. Alignments shorter than this value will be filtered. Note that all alignments shorter than 65 bp will be filtered regardless of this threshold.",
	),
	ToolInput(
		'minContig',
		Int,
		prefix='--min-contig',
		position=1,
		default=500,
		doc="Lower threshold for a contig length (in bp). Shorter contigs won't be taken into account",
	),
	ToolInput(
		'minIdentity',
		Float,
		prefix='--min-identity',
		position=1,
		default=95.0,
		doc="Minimum IDY% considered as proper alignment. Alignments with IDY% worse than this value will be filtered. ote that all alignments with IDY% less than 80.0% will be filtered regardless of this threshold. ",
	),
	ToolInput(
		'operons',
		File(optional=True),
		prefix='--operons',
		position=1,
		default=None,
		doc="Operon positions in the reference genome. Operon coordinates for the reference genome",
	),
	ToolInput(
		'optionO',
		String,
		prefix='-o',
		position=1,
		default="outputdir",
	),
	ToolInput(
		'optionR',
		Array(File(), optional=True),
		prefix='-r',
		position=1,
		default=None,
		doc="Reference genome",
	),
	ToolInput(
		'referencesList',
		String,
		prefix='--references-list',
		position=1,
		default=None,
		doc="Comma-separated list of reference genomes. MetaQUAST will search for these references in the NCBI database and will download the found ones",
	),
	ToolInput(
		'scaffoldGapMaxSize',
		Int,
		prefix='--scaffold-gap-max-size',
		position=1,
		default=1000,
		doc="Max allowed scaffold gap length difference for detecting corresponding type of misassemblies. Longer inconsistencies are considered as relocations and thus, counted as extensive misassemblies. The default value is 10000 bp. Note that the threshold make sense only if it is greater than extensive misassembly size",
	),
	ToolInput(
		'threads',
		Int,
		prefix='--threads',
		position=1,
		default=1,
	),
	ToolInput(
		'unalignedPartSize',
		Int,
		prefix='--unaligned-part-size',
		position=1,
		default=500,
		doc="Lower threshold for detecting partially unaligned contigs",
	),

]

outputs = [
		ToolOutput(
		'outQuastTabular',
		File,
		selector=WildcardSelector(r"outputdir/report.tsv"),
		doc="tabular report",
	),
		ToolOutput(
		'outReportHtml',
		File,
		selector=WildcardSelector(r"outputdir/report.html"),
		doc="HTML report",
	),
		ToolOutput(
		'outReportPdf',
		File,
		selector=WildcardSelector(r"outputdir/report.pdf"),
		doc="PDF report",
	),
		ToolOutput(
		'outLog',
		File,
		selector=WildcardSelector(r"outputdir/quast.log"),
		doc="Log",
	),
		ToolOutput(
		'outMisAss',
		File,
		selector=WildcardSelector(r"outputdir/contigs_reports/misassemblies_report.txt"),
		doc="Misassemblies",
	),
		ToolOutput(
		'outUnalign',
		File,
		selector=WildcardSelector(r"outputdir/contigs_reports/unaligned_report.tsv"),
		doc="Unaligned contigs",
	),
		ToolOutput(
		'outKmers',
		File,
		selector=WildcardSelector(r"outputdir/k_mer_stats/kmers_report.txt"),
		doc="K-mer-based metrics",
	),

]

quast = CommandToolBuilder(
    tool="quast",
    version="5.0.2",
    metadata=metadata,
    container="quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7",
    base_command=['quast', 'metaquast'],
    inputs=inputs,
    outputs=outputs
)


if __name__ == "__main__":
    quast().translate(
        "wdl", to_console=True
    )


