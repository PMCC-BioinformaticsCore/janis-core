# NOTE
# This is an automated translation of the 'multiqc' version '1.7' tool from a Galaxy XML tool wrapper.  
# Translation was performed by the gxtool2janis program (in-development)

import sys
sys.path.append('/home/grace/work/pp/gxtool2janis')

from janis_core.types.common_data_types import String
from janis_core.types.common_data_types import File
from janis_core import CommandToolBuilder
from janis_core import WildcardSelector
from janis_core import ToolMetadata
from janis_core import ToolOutput
from janis_core import ToolInput
from janis_core import Array


metadata = ToolMetadata(
    short_documentation="aggregate results from bioinformatics analyses into a single report",
    keywords=[],
    contributors=['gxtool2janis'],
    dateCreated="2022-08-29 10:49:40",
    dateUpdated="2022-08-29 10:49:40",
    version="1.7",
    doi="https://doi.org/10.1093/bioinformatics/btw354",
    citation="https://doi.org/10.1093/bioinformatics/btw354",
    documentationUrl=None,
    documentation="""
**What it does**

`MultiQC <http://multiqc.info/>`_ aggregates results from bioinformatics analyses across many samples into a single report. It takes results of multiple analyses and creates a report that can be viewed as a single beautiful web-page. It's a general use tool, perfect for summarizing the output from numerous bioinformatics tools.

**Inputs**

MultiQC takes software output summaries/logs and creates a single report from them. You need to tell the tool which software was used to generate the report. This is done using the **Software name** dropdown. At present only the Galaxy tools found in the ToolShed produce logs that can used with MultiQC

----

The first integration of this tool was made by Cyril Monjeaud and Yvan Le Bras (`Enancio <http://enancio.fr/>`_ and Rennes GenOuest Bio-informatics Core Facility). It is now maintained by the `Intergalactic Utilities Commission <https://galaxyproject.org/iuc>`_.
    """
)

inputs = [
	# Positionals

	# Flags

	# Options
	ToolInput(
		'comment',
		String(optional=True),
		prefix='--comment',
		position=1,
		default=None,
		doc="Custom comment. It will be printed at the top of the report",
	),
	ToolInput(
		'config',
		String(optional=True),
		prefix='--config',
		position=1,
		default=None,
	),
	ToolInput(
		'filename',
		String,
		prefix='--filename',
		position=1,
		default="report",
	),
	ToolInput(
		'title',
		String(optional=True),
		prefix='--title',
		position=1,
		default=None,
		doc="Report title. It is printed as page header",
	),

]

outputs = [
		ToolOutput(
		'outStats',
		Array(File),
		selector=WildcardSelector(r"report_data/multiqc_.+\.txt"),
		doc="Stats",
	),
		ToolOutput(
		'outHtmlReport',
		File,
		selector=WildcardSelector(r"report.html"),
		doc="Webpage",
	),
		ToolOutput(
		'outLog',
		File,
		selector=WildcardSelector(r"report_data/multiqc.log"),
		doc="Log",
	),

]

multiqc = CommandToolBuilder(
    tool="multiqc",
    version="1.7",
    metadata=metadata,
    container="quay.io/biocontainers/multiqc:1.7--py_4",
    base_command=['multiqc', 'multiqc_WDir'],
    inputs=inputs,
    outputs=outputs
)


if __name__ == "__main__":
    multiqc().translate(
        "wdl", to_console=True
    )


