# NOTE
# This is an automated translation of the 'Unicycler training (imported from uploaded file)' version '3' workflow. 
# Translation was performed by the gxtool2janis program (in-development)

import sys
sys.path.append('/home/grace/work/pp/gxtool2janis')

from janis_core.types.common_data_types import File
from janis_core import WorkflowMetadata
from janis_core import WorkflowBuilder
from tools.unicycler import unicycler
from tools.multiqc import multiqc
from tools.fastqc import fastqc
from tools.prokka import prokka
from tools.quast import quast
from janis_core import Array

metadata = WorkflowMetadata(
    short_documentation="Unicycler Assembly",
    contributors=['gxtool2janis'],
    keywords=['assembly'],
    dateCreated="2022-08-29 10:49:40",
    dateUpdated="2022-08-29 10:49:40",
    version=3
)

w = WorkflowBuilder(
	"unicyclerTrainingImportedFromUploadedFile",
	version="3",
	doc="Unicycler Assembly",
	metadata=metadata
)

# ------
# INPUTS
# ------

"""
NOTE
Please provide values for these workflow inputs 
in inputs.yaml
"""

w.input("inForwardReads", File)
w.input("inReverseReads", File)
w.input("inLongReads", File)
w.input("fastqc1_adapters", File(optional=True))
w.input("fastqc1_contaminants", File(optional=True))
w.input("fastqc1_limits", File(optional=True))
w.input("fastqc2_adapters", File(optional=True))
w.input("fastqc2_contaminants", File(optional=True))
w.input("fastqc2_limits", File(optional=True))
w.input("prokka_proteins", File(optional=True))

# -------------
# STEP1: FastQC
# -------------

w.step(
	"fastqc1",
	fastqc(
		adapters=w.fastqc1_adapters,          # --adapters [File OPTIONAL] (WORKFLOW INPUT)
		contaminants=w.fastqc1_contaminants,  # --contaminants [File OPTIONAL] (WORKFLOW INPUT)
		limits=w.fastqc1_limits,              # --limits [File OPTIONAL] (WORKFLOW INPUT)
		inputFile=w.inForwardReads,           # [File] (WORKFLOW INPUT)
		extract=False,                        # --extract [Boolean OPTIONAL] (DEFAULT)
		kmers=7,                              # --kmers [Int] (DEFAULT)
		nogroup=False,                        # --nogroup [Boolean OPTIONAL] (DEFAULT)
		quiet=False,                          # --quiet [Boolean OPTIONAL] (DEFAULT)
		minLength=None,                       # --min_length [Int OPTIONAL] (DEFAULT)
		optionF=None,                         # -f [String OPTIONAL] (DEFAULT)
		outdir=None,                          # --outdir [String OPTIONAL] (DEFAULT)
	)
)

# -------------
# STEP2: FastQC
# -------------

w.step(
	"fastqc2",
	fastqc(
		adapters=w.fastqc2_adapters,          # --adapters [File OPTIONAL] (WORKFLOW INPUT)
		contaminants=w.fastqc2_contaminants,  # --contaminants [File OPTIONAL] (WORKFLOW INPUT)
		limits=w.fastqc2_limits,              # --limits [File OPTIONAL] (WORKFLOW INPUT)
		inputFile=w.inReverseReads,           # [File] (WORKFLOW INPUT)
		extract=False,                        # --extract [Boolean OPTIONAL] (DEFAULT)
		kmers=7,                              # --kmers [Int] (DEFAULT)
		nogroup=False,                        # --nogroup [Boolean OPTIONAL] (DEFAULT)
		quiet=False,                          # --quiet [Boolean OPTIONAL] (DEFAULT)
		minLength=None,                       # --min_length [Int OPTIONAL] (DEFAULT)
		optionF=None,                         # -f [String OPTIONAL] (DEFAULT)
		outdir=None,                          # --outdir [String OPTIONAL] (DEFAULT)
	)
)

# ---------------------------------------
# STEP3: Create assemblies with Unicycler
# ---------------------------------------

w.step(
	"unicycler",
	unicycler(
		option11=w.inForwardReads,  # -1 [File] (WORKFLOW INPUT)
		option12=w.inReverseReads,  # -2 [File] (WORKFLOW INPUT)
		optionL=w.inLongReads,      # -l [File OPTIONAL] (WORKFLOW INPUT)
		depthFilter=0.25,           # --depth_filter [Float] (DEFAULT)
		kmerCount=10,               # --kmer_count [Int] (DEFAULT)
		kmers="",                   # --kmers [String OPTIONAL]
		largestComponent=False,     # --largest_component [Boolean OPTIONAL] (DEFAULT)
		linearSeqs=0,               # --linear_seqs [Int] (DEFAULT)
		maxKmerFrac=0.95,           # --max_kmer_frac [Float] (DEFAULT)
		minComponentSize=1000,      # --min_component_size [Int] (DEFAULT)
		minDeadEndSize=1000,        # --min_dead_end_size [Int] (DEFAULT)
		minFastaLength=100,         # --min_fasta_length [Int] (DEFAULT)
		minKmerFrac=0.2,            # --min_kmer_frac [Float] (DEFAULT)
		minPolishSize=1000,         # --min_polish_size [Int] (DEFAULT)
		mode="normal",              # --mode [String] (DEFAULT)
		noCorrect=False,            # --no_correct [Boolean OPTIONAL] (DEFAULT)
		noPilon=False,              # --no_pilon [Boolean OPTIONAL] (DEFAULT)
		noRotate=False,             # --no_rotate [Boolean OPTIONAL] (DEFAULT)
		optionO="./",               # -o [String] (DEFAULT)
		optionT=4,                  # -t [Int] (DEFAULT)
		scores="",                  # --scores [String]
		startGeneCov=95.0,          # --start_gene_cov [Float]
		startGeneId=90.0,           # --start_gene_id [Float]
		verbosity=3,                # --verbosity [Int] (DEFAULT)
		contamination=None,         # --contamination [File OPTIONAL] (DEFAULT)
		lowScore=None,              # --low_score [Int OPTIONAL] (DEFAULT)
		minAnchorSegLen=None,       # --min_anchor_seg_len [Int OPTIONAL] (DEFAULT)
		optionS=None,               # -s [File OPTIONAL] (DEFAULT)
		pilonPath=None,             # --pilon_path [String OPTIONAL] (DEFAULT)
		startGenes=None,            # --start_genes [File OPTIONAL] (DEFAULT)
	)
)

# --------------
# STEP4: MultiQC
# --------------

w.step(
	"multiqc",
	multiqc(
		#UNKNOWN1=w.fastqc1.outTextFile,  # (CONNECTION)
		#UNKNOWN2=w.fastqc2.outTextFile,  # (CONNECTION)
		filename="report",                # --filename [String] (DEFAULT)
		comment=None,                     # --comment [String OPTIONAL] (DEFAULT)
		config=None,                      # --config [String OPTIONAL] (DEFAULT)
		title=None,                       # --title [String OPTIONAL] (DEFAULT)
	)
)

# ------------
# STEP5: Quast
# ------------

w.step(
	"quast",
	quast(
		#UNKNOWN1=w.unicycler.outAssembly,  # (CONNECTION)
		ambiguityUsage="one",               # --ambiguity-usage [String] (DEFAULT)
		circos=False,                       # --circos [Boolean OPTIONAL] (DEFAULT)
		conservedGenesFinding=False,        # --conserved-genes-finding [Boolean OPTIONAL] (DEFAULT)
		eukaryote=False,                    # --eukaryote [Boolean OPTIONAL] (DEFAULT)
		extensiveMisSize=1000,              # --extensive-mis-size [Int] (DEFAULT)
		fragmented=False,                   # --fragmented [Boolean OPTIONAL] (DEFAULT)
		fungus=False,                       # --fungus [Boolean OPTIONAL] (DEFAULT)
		geneFinding=False,                  # --gene-finding [Boolean OPTIONAL] (DEFAULT)
		geneThresholds="0,300,1500,3000",   # --gene-thresholds [String] (DEFAULT)
		glimmer=False,                      # --glimmer [Boolean OPTIONAL] (DEFAULT)
		kMerSize=101,                       # --k-mer-size [Int] (DEFAULT)
		kMerStats=False,                    # --k-mer-stats [Boolean OPTIONAL] (DEFAULT)
		large=False,                        # --large [Boolean OPTIONAL] (DEFAULT)
		maxRefNum=50,                       # --max-ref-num [Int] (DEFAULT)
		mgm=False,                          # --mgm [Boolean OPTIONAL] (DEFAULT)
		minAlignment=65,                    # --min-alignment [Int] (DEFAULT)
		minContig=500,                      # --min-contig [Int] (DEFAULT)
		minIdentity=95.0,                   # --min-identity [Float] (DEFAULT)
		optionO="outputdir",                # -o [String] (DEFAULT)
		referencesList="temp_ref_list_fp",  # --references-list [String]
		rnaFinding=False,                   # --rna-finding [Boolean OPTIONAL] (DEFAULT)
		scaffoldGapMaxSize=1000,            # --scaffold-gap-max-size [Int] (DEFAULT)
		skipUnalignedMisContigs=False,      # --skip-unaligned-mis-contigs [Boolean OPTIONAL] (DEFAULT)
		splitScaffolds=False,               # --split-scaffolds [Boolean OPTIONAL] (DEFAULT)
		strictNa=False,                     # --strict-NA [Boolean OPTIONAL] (DEFAULT)
		testNoRef=False,                    # --test-no-ref [Boolean OPTIONAL] (DEFAULT)
		threads=1,                          # --threads [Int] (DEFAULT)
		unalignedPartSize=500,              # --unaligned-part-size [Int] (DEFAULT)
		useAllAlignments=False,             # --use-all-alignments [Boolean OPTIONAL] (DEFAULT)
		contigThresholds=None,              # --contig-thresholds [String OPTIONAL]
		estRefSize=None,                    # --est-ref-size [Int OPTIONAL] (DEFAULT)
		features=None,                      # --features [File OPTIONAL] (DEFAULT)
		labels=None,                        # --labels [String OPTIONAL] (DEFAULT)
		operons=None,                       # --operons [File OPTIONAL] (DEFAULT)
		optionR=None,                       # -r [File ARRAY OPTIONAL] (DEFAULT)
	)
)

# -------------
# STEP6: Prokka
# -------------

w.step(
	"prokka",
	prokka(
		proteins=w.prokka_proteins,         # --proteins [File OPTIONAL] (WORKFLOW INPUT)
		inputFile=w.unicycler.outAssembly,  # [File] (CONNECTION)
		addgenes=False,                     # --addgenes [Boolean OPTIONAL] (DEFAULT)
		compliant=False,                    # --compliant [Boolean OPTIONAL] (DEFAULT)
		cpus=8,                             # --cpus [Int] (DEFAULT)
		evalue=1e-06,                       # --evalue [Float OPTIONAL] (DEFAULT)
		fast=False,                         # --fast [Boolean OPTIONAL] (DEFAULT)
		gcode=11,                           # --gcode [Int OPTIONAL] (DEFAULT)
		genus="Escherichia",                # --genus [String]
		gffver=3,                           # --gffver [Int] (DEFAULT)
		increment=10,                       # --increment [Int OPTIONAL]
		kingdom="Archaea",                  # --kingdom [String] (DEFAULT)
		locustag="PROKKA",                  # --locustag [String]
		metagenome=False,                   # --metagenome [Boolean OPTIONAL] (DEFAULT)
		mincontig=200,                      # --mincontig [Int OPTIONAL] (DEFAULT)
		norrna=False,                       # --norrna [Boolean OPTIONAL] (DEFAULT)
		outdir="outdir",                    # --outdir [String] (DEFAULT)
		prefix="prokka",                    # --prefix [String] (DEFAULT)
		quiet=False,                        # --quiet [Boolean OPTIONAL] (DEFAULT)
		rfam=False,                         # --rfam [Boolean OPTIONAL] (DEFAULT)
		species="Coli",                     # --species [String]
		strain="C-1",                       # --strain [String]
		usegenus=False,                     # --usegenus [Boolean OPTIONAL] (DEFAULT)
		centre=None,                        # --centre [String OPTIONAL] (DEFAULT)
		plasmid=None,                       # --plasmid [String OPTIONAL] (DEFAULT)
	)
)

# -------
# OUTPUTS
# -------

w.output(
	"fastqc1_outHtmlFile",
	File,
	source=(w.fastqc1, "outHtmlFile")
)

w.output(
	"fastqc1_outTextFile",
	File,
	source=(w.fastqc1, "outTextFile")
)

w.output(
	"fastqc2_outHtmlFile",
	File,
	source=(w.fastqc2, "outHtmlFile")
)

w.output(
	"fastqc2_outTextFile",
	File,
	source=(w.fastqc2, "outTextFile")
)

w.output(
	"unicycler_outAssembly",
	File,
	source=(w.unicycler, "outAssembly")
)

w.output(
	"multiqc_outStats",
	Array(File),
	source=(w.multiqc, "outStats")
)

w.output(
	"multiqc_outHtmlReport",
	File,
	source=(w.multiqc, "outHtmlReport")
)

w.output(
	"quast_outReportHtml",
	File,
	source=(w.quast, "outReportHtml")
)

w.output(
	"prokka_outGff",
	File,
	source=(w.prokka, "outGff")
)

w.output(
	"prokka_outGbk",
	File,
	source=(w.prokka, "outGbk")
)

w.output(
	"prokka_outFna",
	File,
	source=(w.prokka, "outFna")
)

w.output(
	"prokka_outFaa",
	File,
	source=(w.prokka, "outFaa")
)

w.output(
	"prokka_outFfn",
	File,
	source=(w.prokka, "outFfn")
)

w.output(
	"prokka_outSqn",
	File,
	source=(w.prokka, "outSqn")
)

w.output(
	"prokka_outFsa",
	File,
	source=(w.prokka, "outFsa")
)

w.output(
	"prokka_outTbl",
	File,
	source=(w.prokka, "outTbl")
)

w.output(
	"prokka_outErr",
	File,
	source=(w.prokka, "outErr")
)

w.output(
	"prokka_outTxt",
	File,
	source=(w.prokka, "outTxt")
)

w.output(
	"prokka_outLog",
	File,
	source=(w.prokka, "outLog")
)

if __name__ == "__main__":
    w.translate("cwl", **args)

