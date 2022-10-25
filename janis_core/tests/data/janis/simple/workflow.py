

from janis_core.types.common_data_types import File
from janis_core import WorkflowMetadata
from janis_core import WorkflowBuilder
from janis_core import Array

from janis_core.tests.data.janis.simple.tools.unicycler import unicycler
from janis_core.tests.data.janis.simple.tools.multiqc import multiqc
from janis_core.tests.data.janis.simple.tools.fastqc import fastqc
from janis_core.tests.data.janis.simple.tools.prokka import prokka
from janis_core.tests.data.janis.simple.tools.quast import quast

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
		kmers="",                   # --kmers [String OPTIONAL]
		scores="",                  # --scores [String]
		startGeneCov=95.0,          # --start_gene_cov [Float]
		startGeneId=90.0,           # --start_gene_id [Float]
	)
)

# --------------
# STEP4: MultiQC
# --------------

w.step(
	"multiqc",
	multiqc(
		unknown1=w.fastqc1.outTextFile,  # (CONNECTION)
		unknown2=w.fastqc2.outTextFile,  # (CONNECTION)
	)
)

# ------------
# STEP5: Quast
# ------------

w.step(
	"quast",
	quast(
		unknown1=w.unicycler.outAssembly,  # (CONNECTION)
		referencesList="temp_ref_list_fp",  # --references-list [String]
		contigThresholds=None,              # --contig-thresholds [String OPTIONAL]
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
		genus="Escherichia",                # --genus [String]
		increment=10,                       # --increment [Int OPTIONAL]
		locustag="PROKKA",                  # --locustag [String]
		species="Coli",                     # --species [String]
		strain="C-1",                       # --strain [String]
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

