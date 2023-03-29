

from janis_core import WorkflowMetadata
from janis_core import WorkflowBuilder
from janis_core import ScatterDescription
from janis_core import ScatterMethods
from janis_core import Array
from janis_core import File, Int

from janis_core.tests.testtools import FastqcTestTool
from janis_core.tests.testtools import UnicyclerTestTool
from janis_core.tests.testtools import CatTestTool



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
w.input("testInput", File)
w.input("fastqc1_adapters", File(optional=True))
w.input("fastqc1_contaminants", File(optional=True))
w.input("fastqc1_limits", File(optional=True))
w.input("fastqc2_adapters", File(optional=True))
w.input("fastqc2_contaminants", File(optional=True))
w.input("fastqc2_limits", File(optional=True))
w.input("min_len", Int(optional=True))

# -------------
# STEP1: FastQC
# -------------

w.step(
	"fastqc1",
	FastqcTestTool(
		adapters=w.fastqc1_adapters,          # --adapters [File OPTIONAL] (WORKFLOW INPUT)
		contaminants=w.fastqc1_contaminants,  # --contaminants [File OPTIONAL] (WORKFLOW INPUT)
		limits=w.fastqc1_limits,              # --limits [File OPTIONAL] (WORKFLOW INPUT)
		inputFile=w.inForwardReads,           # [File] (WORKFLOW INPUT)
	),
)

# -------------
# STEP2: FastQC
# -------------

w.step(
	"fastqc2",
	FastqcTestTool(
		adapters=w.fastqc2_adapters,          # --adapters [File OPTIONAL] (WORKFLOW INPUT)
		contaminants=w.fastqc2_contaminants,  # --contaminants [File OPTIONAL] (WORKFLOW INPUT)
		limits=w.fastqc2_limits,              # --limits [File OPTIONAL] (WORKFLOW INPUT)
		inputFile=w.inReverseReads,           # [File] (WORKFLOW INPUT)
	),
)

### TEST PURPOSES

w.step(
	"fastqc3",
	FastqcTestTool(
		inputFile=w.testInput,                  # [File] (WORKFLOW INPUT)
	),
)

w.step(
    "CatTestTool",
    CatTestTool(
        inp=w.fastqc3.outTextFile,
    )

)



# ---------------------------------------
# STEP3: Create assemblies with Unicycler
# ---------------------------------------

w.step(
	"unicycler",
	UnicyclerTestTool(
		option1=w.inForwardReads,  # -1 [File] (WORKFLOW INPUT)
		option2=w.inReverseReads,  # -2 [File] (WORKFLOW INPUT)
		optionL=w.inLongReads,      # -l [File OPTIONAL] (WORKFLOW INPUT)
		kmers="",                   # --kmers [String OPTIONAL]
		scores="",                  # --scores [String]
		startGeneCov=95.0,          # --start_gene_cov [Float]
		startGeneId=90.0,           # --start_gene_id [Float]
	),
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


if __name__ == "__main__":
    w.translate("cwl", **args)




# metadata = WorkflowMetadata(
#     short_documentation="Unicycler Assembly",
#     contributors=['gxtool2janis'],
#     keywords=['assembly'],
#     dateCreated="2022-08-29 10:49:40",
#     dateUpdated="2022-08-29 10:49:40",
#     version=3
# )

# w = WorkflowBuilder(
# 	"unicyclerTrainingImportedFromUploadedFile",
# 	version="3",
# 	doc="Unicycler Assembly",
# 	metadata=metadata
# )

# # ------
# # INPUTS
# # ------

# w.input("inForwardReads", Array(File))
# w.input("inReverseReads", Array(File))
# w.input("inLongReads", Array(File))
# w.input("testInput", File)
# w.input("fastqc1_adapters", File(optional=True))
# w.input("fastqc1_contaminants", File(optional=True))
# w.input("fastqc1_limits", File(optional=True))
# w.input("fastqc2_adapters", File(optional=True))
# w.input("fastqc2_contaminants", File(optional=True))
# w.input("fastqc2_limits", File(optional=True))

# # -------------
# # STEP1: FastQC
# # -------------

# w.step(
# 	"fastqc1",
# 	FastqcTestTool(
# 		adapters=w.fastqc1_adapters,          # --adapters [File OPTIONAL] (WORKFLOW INPUT)
# 		contaminants=w.fastqc1_contaminants,  # --contaminants [File OPTIONAL] (WORKFLOW INPUT)
# 		limits=w.fastqc1_limits,              # --limits [File OPTIONAL] (WORKFLOW INPUT)
# 		inputFile=w.inForwardReads,           # [File] (WORKFLOW INPUT)
# 	),
#     scatter="inputFile"
# )

# # -------------
# # STEP2: FastQC
# # -------------

# w.step(
# 	"fastqc2",
# 	FastqcTestTool(
# 		adapters=w.fastqc2_adapters,          # --adapters [File OPTIONAL] (WORKFLOW INPUT)
# 		contaminants=w.fastqc2_contaminants,  # --contaminants [File OPTIONAL] (WORKFLOW INPUT)
# 		limits=w.fastqc2_limits,              # --limits [File OPTIONAL] (WORKFLOW INPUT)
# 		inputFile=w.inReverseReads,           # [File] (WORKFLOW INPUT)
# 	),
#     scatter="inputFile"
# )

# ### TEST PURPOSES

# w.step(
# 	"fastqc3",
# 	FastqcTestTool(
# 		inputFile=w.testInput,                  # [File] (WORKFLOW INPUT)
# 	),
# )

# w.step(
#     "CatTestTool",
#     CatTestTool(
#         inp=w.fastqc3.outTextFile,
#     )

# )

# # ---------------------------------------
# # STEP3: Create assemblies with Unicycler
# # ---------------------------------------

# w.step(
# 	"unicycler",
# 	UnicyclerTestTool(
# 		option11=w.inForwardReads,  # -1 [File] (WORKFLOW INPUT)
# 		option12=w.inReverseReads,  # -2 [File] (WORKFLOW INPUT)
# 		optionL=w.inLongReads,      # -l [File OPTIONAL] (WORKFLOW INPUT)
# 		kmers="",                   # --kmers [String OPTIONAL]
# 		scores="",                  # --scores [String]
# 		startGeneCov=95.0,          # --start_gene_cov [Float]
# 		startGeneId=90.0,           # --start_gene_id [Float]
# 	),
#     scatter=ScatterDescription(fields=['option11', 'option12', 'optionL'], method=ScatterMethods.dot)
# )


# # -------
# # OUTPUTS
# # -------

# w.output(
# 	"fastqc1_outHtmlFile",
# 	Array(File),
# 	source=(w.fastqc1, "outHtmlFile")
# )

# w.output(
# 	"fastqc1_outTextFile",
# 	Array(File),
# 	source=(w.fastqc1, "outTextFile")
# )

# w.output(
# 	"fastqc2_outHtmlFile",
# 	Array(File),
# 	source=(w.fastqc2, "outHtmlFile")
# )

# w.output(
# 	"fastqc2_outTextFile",
# 	Array(File),
# 	source=(w.fastqc2, "outTextFile")
# )

# w.output(
# 	"unicycler_outAssembly",
# 	Array(File),
# 	source=(w.unicycler, "outAssembly")
# )


# if __name__ == "__main__":
#     w.translate("cwl", **args)

