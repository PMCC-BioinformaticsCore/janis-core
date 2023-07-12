# NOTE
# This is an automated translation of the 'unicycler' version '0.4.8.0' tool from a Galaxy XML tool wrapper.  
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


metadata = ToolMetadata(
    short_documentation="",
    keywords=[],
    contributors=['gxtool2janis'],
    dateCreated="2022-08-29 10:49:40",
    dateUpdated="2022-08-29 10:49:40",
    version="0.4.8.0",
    doi="https://doi.org/10.1101/096412",
    citation="https://doi.org/10.1101/096412",
    documentationUrl=None,
    documentation="""

**Unicycler**

Unicycler is a hybrid assembly pipeline for bacterial genomes. It uses both Illumina reads and long reads (PacBio or Nanopore) to produce complete and accurate assemblies. It is written by `Ryan Wick`_ at the University of Melbourne's Centre for Systems Genomics. Much of the description below is lifted from Unicycler's `github page`_.

.. _`Ryan Wick`: https://github.com/rrwick
.. _`github page`: https://github.com/rrwick/Unicycler

-----

**Input data**

Unicycler accepts inputs short (Illumina) reads in FASTQ format. Galaxy places additional requirement of having these in FASTQ format with `Sanger encoding`_ of quality scores. Long reads (from Oxford Nanopore or PacBio) can be either in FASTQ of FASTA form.

.. _`Sanger encoding`: https://en.wikipedia.org/wiki/FASTQ_format#Quality

The input options are::

    -1 SHORT1, --short1 SHORT1
        FASTQ file of short reads (first reads in each pair)
    -2 SHORT2, --short2 SHORT2
        FASTQ file of short reads (second reads in each pair)
    -s SHORT_UNPAIRED, --short_unpaired SHORT_UNPAIRED
        FASTQ file of unpaired short reads
    -l LONG, --long LONG
        FASTQ or FASTA file of long reads, if all reads are available at start.

-----

**Bridging mode**

Unicycler can be run in three modes: conservative, normal (the default) and bold, set with the --mode option. Conservative mode is least likely to produce a complete assembly but has a very low risk of misassembly. Bold mode is most likely to produce a complete assembly but carries greater risk of misassembly. Normal mode is intermediate regarding both completeness and misassembly risk. See `description of modes`_ for more information.

.. _`description of modes`: https://github.com/rrwick/Unicycler#conservative-normal-and-bold

The available modes are::

    --mode {conservative,normal,bold}
        Bridging mode (default: normal)
        conservative = smaller contigs, lowest misassembly rate
        normal = moderate contig size and misassembly rate
        bold = longest contigs, higher misassembly rate

----

**Skip SPAdes error correction step**

Sequencing data contains a substantial number of sequencing errors that manifest themselves as deviations (bulges and non-connected components) within the assembly graph. One of the ways to improve the graph even constructing it is to minimize the amount sequencing errors by performing error correction. SPAdes, which is used by Unicycler for error correction and assembly, uses `BayesHammer`_ to correct the reads. Here is a brief summary of what it does:

 1. SPAdes (or rather BayesHammer) counts *k*-mers in reads and computed *k*-mer statistics that takes into account base quality values.
 2. `Hamming graph`_ is constructed for *k*-mers is which *k*-mers are nodes. In this graph edges connect nodes (*k*-mers) is they differ from each other by a number of nucleotides up to a certain threshold (the `Hamming distance`_). The graph is central to the error correction algorithm.
 3. At this step Bayesian subclustering of the graph produced in the previous step. For each *k*-mer we now know the center of its subcluster.
 4. Solid *k*-mers are derived from cluster centers and are assumed to be *error free*.
 5. Solid *k*-mers are mapped back to the reads and used to correct them.

This step takes considerable time, so if one need to quickly evaluate assemblies this step can be skipped. However, this is not recommended if one if trying to produce a final high quality assembly.

.. _`BayesHammer`: https://goo.gl/1iGkMe
.. _`Hamming graph`: https://en.wikipedia.org/wiki/Hamming_graph
.. _`Hamming distance`: https://en.wikipedia.org/wiki/Hamming_distance

This following option turns error correction on and off::

    --no_correct
        Skip SPAdes error correction step
        (default: conduct SPAdes error correction)

-----

**Do not rotate completed replicons to start at a standard gene**

Unicycler uses TBLASTN to search for dnaA or repA alleles in each completed replicon. If one is found, the sequence is rotated and/or flipped so that it begins with that gene encoded on the forward strand. This provides consistently oriented assemblies and reduces the risk that a gene will be split across the start and end of the sequence.

The following option turns rotation on and off::

    --no_rotate
        Do not rotate completed replicons
        to start at a standard gene
        (default: completed replicons are rotated)

**Do not use Pilon to polish the final assembly**

`Pilon`_ is a tool for improving overall quality of draft assemblies and finding variation among strains. Unicycler uses it for assembly *polishing*.

The following option turns pilon part of Unicycler pipeline on and off::

    --no_pilon
        Do not use Pilon to polish the
        final assembly (default: Pilon is used)

.. _`Pilon`: https://github.com/broadinstitute/pilon/wiki

------

**Expected number of linear sequences**

If you expect your sample to contain linear (non circular) sequences, set this option::

    --linear_seqs EXPECTED_LINEAR_SEQS
        The expected number of linear (i.e. non-circular)
        sequences in the underlying sequence

----

**SPAdes options**

This section provides control of SPAdes options::

    --min_kmer_frac MIN_KMER_FRAC
        Lowest k-mer size for SPAdes assembly,
        expressed as a fraction of the read length
        (default: 0.2)
    --max_kmer_frac MAX_KMER_FRAC
        Highest k-mer size for SPAdes assembly,
        expressed as a fraction of the read length
        (default: 0.95)
    --kmer_count KMER_COUNT
        Number of k-mer steps to use in
        SPAdes assembly (default: 10)
    --depth_filter DEPTH_FILTER
        Filter out contigs lower than this fraction
        of the chromosomal depth, if doing so does
        not result in graph dead ends (default: 0.25)

----

**Rotation options**

Unicycler attempts to rotate circular assemblies to make sure that they begin at a consistent starting gene. The following parameters control assembly rotation::

    --start_genes START_GENES
        FASTA file of genes for start point
        of rotated replicons
        (default: start_genes.fasta)
    --start_gene_id START_GENE_ID
        The minimum required BLAST percent identity
        for a start gene search
        (default: 90.0)
    --start_gene_cov START_GENE_COV
        The minimum required BLAST percent coverage
        for a start gene search
        (default: 95.0)

-----

**Graph cleaning options**

These options control the removal of small leftover sequences after bridging is complete::

    --min_component_size MIN_COMPONENT_SIZE
        Unbridged graph components smaller
        than this size (bp) will be removed
        from the final graph (default: 1000)
    --min_dead_end_size MIN_DEAD_END_SIZE
        Graph dead ends smaller than this size (bp)
        will be removed from the final graph
        (default: 1000)

-----

**Long read alignment options**

These options control the alignment of long reads to the assembly graph::

    --contamination CONTAMINATION
        FASTA file of known contamination in long reads
    --scores SCORES
        Comma-delimited string of alignment scores:
        match, mismatch, gap open, gap extend
        (default: 3,-6,-5,-2)
    --low_score LOW_SCORE
        Score threshold - alignments below this
        are considered poor
        (default: set threshold automatically)

-----

**Outputs**

Galaxy's wrapped for Unicycler produces two outputs:

 * final assembly in FASTA format
 * final assembly grapth in graph format

 While most will likely be interested in the FASTA dataset, the graph dataset is also quite useful and can be visualized using tools such as `Bandage`_.


.. _`Bandage`: https://github.com/rrwick/Bandage


    """
)

inputs = [
	# Positionals

	# Flags
	ToolInput(
		'largestComponent',
		Boolean(optional=True),
		prefix='--largest_component',
		position=1,
		default="False",
	),
	ToolInput(
		'noCorrect',
		Boolean(optional=True),
		prefix='--no_correct',
		position=1,
		default="False",
	),
	ToolInput(
		'noPilon',
		Boolean(optional=True),
		prefix='--no_pilon',
		position=1,
		default="False",
	),
	ToolInput(
		'noRotate',
		Boolean(optional=True),
		prefix='--no_rotate',
		position=1,
		default="False",
	),

	# Options
	ToolInput(
		'contamination',
		File(optional=True),
		prefix='--contamination',
		position=1,
		default=None,
		doc="FASTA file of known contamination in long reads, e.g. lambda, phiXm or puc18 spike-ins.",
	),
	ToolInput(
		'depthFilter',
		Float,
		prefix='--depth_filter',
		position=1,
		default=0.25,
		doc="Filter out contigs lower than this fraction of the chromosomal depth. It is done if does not result in graph dead ends",
	),
	ToolInput(
		'kmerCount',
		Int,
		prefix='--kmer_count',
		position=1,
		default=10,
		doc="Number of k-mer steps to use in SPAdes assembly",
	),
	ToolInput(
		'kmers',
		String(optional=True),
		prefix='--kmers',
		position=1,
		default=None,
		doc="Exact k-mers to use for SPAdes assembly, comma-separated",
	),
	ToolInput(
		'linearSeqs',
		Int,
		prefix='--linear_seqs',
		position=1,
		default=0,
		doc="The expected number of linear (i.e. non-circular) sequences in the assembly",
	),
	ToolInput(
		'lowScore',
		Int(optional=True),
		prefix='--low_score',
		position=1,
		default=None,
		doc="Score threshold - alignments below this are considered poor. default = set automatically",
	),
	ToolInput(
		'maxKmerFrac',
		Float,
		prefix='--max_kmer_frac',
		position=1,
		default=0.95,
		doc="Highest k-mer size for SPAdes assembly, expressed as a fraction of the read length",
	),
	ToolInput(
		'minAnchorSegLen',
		Int(optional=True),
		prefix='--min_anchor_seg_len',
		position=1,
		default=None,
		doc="Unicycler will not use segments shorter than this as scaffolding anchors",
	),
	ToolInput(
		'minComponentSize',
		Int,
		prefix='--min_component_size',
		position=1,
		default=1000,
		doc="Unbridged graph components smaller than this size will be removed from the final graph",
	),
	ToolInput(
		'minDeadEndSize',
		Int,
		prefix='--min_dead_end_size',
		position=1,
		default=1000,
		doc="Graph dead ends smaller than this size will be removed from the final graph",
	),
	ToolInput(
		'minFastaLength',
		Int,
		prefix='--min_fasta_length',
		position=1,
		default=100,
		doc="Exclude contigs from the FASTA file which are shorter than this length (bp)",
	),
	ToolInput(
		'minKmerFrac',
		Float,
		prefix='--min_kmer_frac',
		position=1,
		default=0.2,
		doc="Lowest k-mer size for SPAdes assembly, expressed as a fraction of the read length",
	),
	ToolInput(
		'minPolishSize',
		Int,
		prefix='--min_polish_size',
		position=1,
		default=1000,
		doc="Contigs shorter than this value (bp) will not be polished using Pilon",
	),
	ToolInput(
		'mode',
		String,
		prefix='--mode',
		position=1,
		default="normal",
		doc="Select Bridging mode. possible values: bold, conservative, normal",
	),
	ToolInput(
		'option11',
		File,
		prefix='-1',
		position=1,
		default=None,
		doc="Select first set of reads. Specify dataset with forward reads",
	),
	ToolInput(
		'option12',
		File,
		prefix='-2',
		position=1,
		default=None,
		doc="Select second set of reads. Specify dataset with reverse reads",
	),
	ToolInput(
		'optionL',
		File(optional=True),
		prefix='-l',
		position=1,
		default=None,
		doc="Select long reads. If there are no long reads, leave this empty",
	),
	ToolInput(
		'optionO',
		String,
		prefix='-o',
		position=1,
		default="./",
	),
	ToolInput(
		'optionS',
		File(optional=True),
		prefix='-s',
		position=1,
		default=None,
		doc="Select unpaired reads. Specify dataset with unpaired reads",
	),
	ToolInput(
		'optionT',
		Int,
		prefix='-t',
		position=1,
		default=4,
	),
	ToolInput(
		'pilonPath',
		String(optional=True),
		prefix='--pilon_path',
		position=1,
		default=None,
	),
	ToolInput(
		'scores',
		String,
		prefix='--scores',
		position=1,
		default="3,-6,-5,-2",
		doc="Comma-delimited string of alignment scores: match, mismatch, gap open, gap extend",
	),
	ToolInput(
		'startGeneCov',
		Float,
		prefix='--start_gene_cov',
		position=1,
		default=95,
		doc="The minimum required BLAST percent coverage for a start gene search",
	),
	ToolInput(
		'startGeneId',
		Float,
		prefix='--start_gene_id',
		position=1,
		default=90,
		doc="The minimum required BLAST percent identity for a start gene search",
	),
	ToolInput(
		'startGenes',
		File(optional=True),
		prefix='--start_genes',
		position=1,
		default=None,
		doc="FASTA file of genes for start point of rotated replicons",
	),
	ToolInput(
		'verbosity',
		Int,
		prefix='--verbosity',
		position=1,
		default=3,
	),

]

outputs = [
		ToolOutput(
		'outAssemblyGraph',
		File,
		selector=WildcardSelector(r"assembly.gfa"),
		doc="Final Assembly Graph",
	),
		ToolOutput(
		'outAssembly',
		File,
		selector=WildcardSelector(r"assembly.fasta"),
		doc="Final Assembly",
	),

]

unicycler = CommandToolBuilder(
    tool="unicycler",
    version="0.4.8.0",
    metadata=metadata,
    container="quay.io/biocontainers/unicycler:0.4.8--py39h98c8e45_5",
    base_command=['unicycler'],
    inputs=inputs,
    outputs=outputs
)


if __name__ == "__main__":
    unicycler().translate(
        "wdl", to_console=True
    )


