# NOTE
# This is an automated translation of the 'prokka' version '1.14.5' tool from a Galaxy XML tool wrapper.  
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
    short_documentation="Prokaryotic genome annotation",
    keywords=[],
    contributors=['gxtool2janis'],
    dateCreated="2022-08-29 10:49:40",
    dateUpdated="2022-08-29 10:49:40",
    version="1.14.5",
    doi="https://doi.org/10.1093/bioinformatics/btu135",
    citation="https://doi.org/10.1093/bioinformatics/btu135",
    documentationUrl=None,
    documentation="""
**What it does**

Prokka_ is a software tool to rapidly annotate bacterial, archaeal and viral genomes, and produce output files that require only minor tweaking to submit to GenBank/ENA/DDBJ.

.. _Prokka: http://github.com/tseemann/prokka

**Output files**

Prokka creates several output files, which are described in the **Additional outputs** section above.

**License and citation**

This Galaxy tool is Copyright © 2013 Lionel Guy, © 2013-2014 `CRS4 Srl.`_, © 2015-2016 `Earlham Institute`_, 2018 `Galaxy IUC` and is released under the `MIT license`_.

.. _CRS4 Srl.: http://www.crs4.it/
.. _Earlham Institute: http://earlham.ac.uk/
.. _MIT license: https://opensource.org/licenses/MIT

You can use this tool only if you agree to the license terms of: `Prokka`_.

.. _Prokka: http://github.com/tseemann/prokka
    """
)

inputs = [
	# Positionals
	ToolInput(
		'inputFile',
		File,
		position=2,
		default=None,
		doc="Contigs to annotate. FASTA format",
	),

	# Flags
	ToolInput(
		'addgenes',
		Boolean(optional=True),
		prefix='--addgenes',
		position=1,
		default="False",
	),
	ToolInput(
		'compliant',
		Boolean(optional=True),
		prefix='--compliant',
		position=1,
		default="False",
	),
	ToolInput(
		'fast',
		Boolean(optional=True),
		prefix='--fast',
		position=1,
		default="False",
	),
	ToolInput(
		'metagenome',
		Boolean(optional=True),
		prefix='--metagenome',
		position=1,
		default="False",
	),
	ToolInput(
		'norrna',
		Boolean(optional=True),
		prefix='--norrna',
		position=1,
		default="False",
	),
	ToolInput(
		'quiet',
		Boolean(optional=True),
		prefix='--quiet',
		position=1,
		default="False",
	),
	ToolInput(
		'rfam',
		Boolean(optional=True),
		prefix='--rfam',
		position=1,
		default="False",
	),
	ToolInput(
		'usegenus',
		Boolean(optional=True),
		prefix='--usegenus',
		position=1,
		default="False",
	),

	# Options
	ToolInput(
		'centre',
		String(optional=True),
		prefix='--centre',
		position=1,
		default=None,
		doc="Sequencing centre ID (--centre)",
	),
	ToolInput(
		'cpus',
		Int,
		prefix='--cpus',
		position=1,
		default=8,
	),
	ToolInput(
		'evalue',
		Float(optional=True),
		prefix='--evalue',
		position=1,
		default=1e-06,
		doc="Similarity e-value cut-off",
	),
	ToolInput(
		'gcode',
		Int(optional=True),
		prefix='--gcode',
		position=1,
		default=11,
		doc="Genetic code (transl_table)",
	),
	ToolInput(
		'genus',
		String,
		prefix='--genus',
		position=1,
		default=None,
		doc="Genus name (--genus). May be used to aid annotation, see --usegenus below",
	),
	ToolInput(
		'gffver',
		Int,
		prefix='--gffver',
		position=1,
		default=3,
	),
	ToolInput(
		'increment',
		Int(optional=True),
		prefix='--increment',
		position=1,
		default=1,
		doc="Locus tag counter increment (--increment)",
	),
	ToolInput(
		'kingdom',
		String,
		prefix='--kingdom',
		position=1,
		default="Archaea",
	),
	ToolInput(
		'locustag',
		String,
		prefix='--locustag',
		position=1,
		default=None,
		doc="Locus tag prefix (--locustag)",
	),
	ToolInput(
		'mincontig',
		Int(optional=True),
		prefix='--mincontig',
		position=1,
		default=200,
		doc="Minimum contig size (--mincontiglen). NCBI needs 200",
	),
	ToolInput(
		'outdir',
		String,
		prefix='--outdir',
		position=1,
		default="outdir",
	),
	ToolInput(
		'plasmid',
		String(optional=True),
		prefix='--plasmid',
		position=1,
		default=None,
		doc="Plasmid name or identifier (--plasmid)",
	),
	ToolInput(
		'prefix',
		String,
		prefix='--prefix',
		position=1,
		default="prokka",
	),
	ToolInput(
		'proteins',
		File(optional=True),
		prefix='--proteins',
		position=1,
		default=None,
		doc="Optional FASTA file of trusted proteins to first annotate from (--proteins)",
	),
	ToolInput(
		'species',
		String,
		prefix='--species',
		position=1,
		default=None,
		doc="Species name (--species)",
	),
	ToolInput(
		'strain',
		String,
		prefix='--strain',
		position=1,
		default=None,
		doc="Strain name (--strain)",
	),

]

outputs = [
		ToolOutput(
		'outGff',
		File,
		selector=WildcardSelector(r"outdir/prokka.gff"),
		doc="gff",
	),
		ToolOutput(
		'outGbk',
		File,
		selector=WildcardSelector(r"outdir/prokka.gbk"),
		doc="gbk",
	),
		ToolOutput(
		'outFna',
		File,
		selector=WildcardSelector(r"outdir/prokka.fna"),
		doc="fna",
	),
		ToolOutput(
		'outFaa',
		File,
		selector=WildcardSelector(r"outdir/prokka.faa"),
		doc="faa",
	),
		ToolOutput(
		'outFfn',
		File,
		selector=WildcardSelector(r"outdir/prokka.ffn"),
		doc="ffn",
	),
		ToolOutput(
		'outSqn',
		File,
		selector=WildcardSelector(r"outdir/prokka.sqn"),
		doc="sqn",
	),
		ToolOutput(
		'outFsa',
		File,
		selector=WildcardSelector(r"outdir/prokka.fsa"),
		doc="fsa",
	),
		ToolOutput(
		'outTbl',
		File,
		selector=WildcardSelector(r"outdir/prokka.tbl"),
		doc="tbl",
	),
		ToolOutput(
		'outTsv',
		File,
		selector=WildcardSelector(r"outdir/prokka.tsv"),
		doc="tsv",
	),
		ToolOutput(
		'outErr',
		File,
		selector=WildcardSelector(r"outdir/prokka.err"),
		doc="err",
	),
		ToolOutput(
		'outTxt',
		File,
		selector=WildcardSelector(r"outdir/prokka.txt"),
		doc="txt",
	),
		ToolOutput(
		'outLog',
		File,
		selector=WildcardSelector(r"outdir/prokka.log"),
		doc="log",
	),

]

prokka = CommandToolBuilder(
    tool="prokka",
    version="1.14.5",
    metadata=metadata,
    container="quay.io/biocontainers/prokka:1.14.5--pl526_1",
    base_command=['prokka'],
    inputs=inputs,
    outputs=outputs
)


if __name__ == "__main__":
    prokka().translate(
        "wdl", to_console=True
    )


