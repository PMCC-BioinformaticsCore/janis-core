
from janis_core import File


class Ab1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ab1"
		)


	@staticmethod
	def name():
		return "Ab1"


	def doc(self):
		return (
			"A binary sequence file in 'ab1' format with a '.ab1' file extension.  You must manually select this 'File Format' when uploading the file."
		)






class AccNos(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.accnos", alternate_extensions={".mothur.otulabels"}
		)


	@staticmethod
	def name():
		return "AccNos"


	def doc(self):
		return (
			""
		)






class Acedb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".acedb"
		)


	@staticmethod
	def name():
		return "Acedb"


	def doc(self):
		return (
			""
		)






class Affybatch(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".affybatch"
		)


	@staticmethod
	def name():
		return "Affybatch"


	def doc(self):
		return (
			""
		)






class AlignCheck(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.align.check"
		)


	@staticmethod
	def name():
		return "AlignCheck"


	def doc(self):
		return (
			""
		)






class AlignReport(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.align.report"
		)


	@staticmethod
	def name():
		return "AlignReport"


	def doc(self):
		return (
			""
		)






class AllegroLOD(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".allegro_fparam"
		)


	@staticmethod
	def name():
		return "AllegroLOD"


	def doc(self):
		return (
			""
		)






class Alohomora_maf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".alohomora_maf"
		)


	@staticmethod
	def name():
		return "Alohomora_maf"


	def doc(self):
		return (
			""
		)






class Alohomora_map(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".alohomora_map"
		)


	@staticmethod
	def name():
		return "Alohomora_map"


	def doc(self):
		return (
			""
		)






class Alohomora_ped(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".alohomora_ped"
		)


	@staticmethod
	def name():
		return "Alohomora_ped"


	def doc(self):
		return (
			""
		)






class Amos(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".afg"
		)


	@staticmethod
	def name():
		return "Amos"


	def doc(self):
		return (
			""
		)






class Analyze75(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".analyze75"
		)


	@staticmethod
	def name():
		return "Analyze75"


	def doc(self):
		return (
			""
		)






class Anndata(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".h5ad"
		)


	@staticmethod
	def name():
		return "Anndata"


	def doc(self):
		return (
			"An HDF5-based anndata File"
		)






class AnvioComposite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_composite"
		)


	@staticmethod
	def name():
		return "AnvioComposite"


	def doc(self):
		return (
			""
		)






class AnvioContigsDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_contigs_db"
		)


	@staticmethod
	def name():
		return "AnvioContigsDB"


	def doc(self):
		return (
			""
		)






class AnvioDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_db"
		)


	@staticmethod
	def name():
		return "AnvioDB"


	def doc(self):
		return (
			""
		)






class AnvioGenomesDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_genomes_db"
		)


	@staticmethod
	def name():
		return "AnvioGenomesDB"


	def doc(self):
		return (
			""
		)






class AnvioPanDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_pan_db"
		)


	@staticmethod
	def name():
		return "AnvioPanDB"


	def doc(self):
		return (
			""
		)






class AnvioProfileDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_profile_db"
		)


	@staticmethod
	def name():
		return "AnvioProfileDB"


	def doc(self):
		return (
			""
		)






class AnvioSamplesDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_samples_db"
		)


	@staticmethod
	def name():
		return "AnvioSamplesDB"


	def doc(self):
		return (
			""
		)






class AnvioStructureDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_structure_db"
		)


	@staticmethod
	def name():
		return "AnvioStructureDB"


	def doc(self):
		return (
			""
		)






class Anvio_classifier(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_classifier"
		)


	@staticmethod
	def name():
		return "Anvio_classifier"


	def doc(self):
		return (
			""
		)






class Anvio_cog_profile(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_cog_profile"
		)


	@staticmethod
	def name():
		return "Anvio_cog_profile"


	def doc(self):
		return (
			""
		)






class Anvio_pfam_profile(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_pfam_profile"
		)


	@staticmethod
	def name():
		return "Anvio_pfam_profile"


	def doc(self):
		return (
			""
		)






class Anvio_state(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_state"
		)


	@staticmethod
	def name():
		return "Anvio_state"


	def doc(self):
		return (
			""
		)






class Anvio_variability(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".anvio_variability"
		)


	@staticmethod
	def name():
		return "Anvio_variability"


	def doc(self):
		return (
			""
		)






class Arff(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".arff"
		)


	@staticmethod
	def name():
		return "Arff"


	def doc(self):
		return (
			""
		)






class Augustus(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".augustus"
		)


	@staticmethod
	def name():
		return "Augustus"


	def doc(self):
		return (
			""
		)






class Axes(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.axes"
		)


	@staticmethod
	def name():
		return "Axes"


	def doc(self):
		return (
			""
		)






class Axt(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".axt"
		)


	@staticmethod
	def name():
		return "Axt"


	def doc(self):
		return (
			"blastz pairwise alignment format.  Each alignment block in an axt file contains three lines: a summary line and 2 sequence lines.  Blocks are separated from one another by blank lines.  The summary line contains chromosomal position and size information about the alignment. It consists of 9 required fields."
		)






class BPF(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".par"
		)


	@staticmethod
	def name():
		return "BPF"


	def doc(self):
		return (
			""
		)






class BafTar(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".brukerbaf.d.tar"
		)


	@staticmethod
	def name():
		return "BafTar"


	def doc(self):
		return (
			""
		)






class Bai(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bai"
		)


	@staticmethod
	def name():
		return "Bai"


	def doc(self):
		return (
			""
		)






class Bam(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bam"
		)


	@staticmethod
	def name():
		return "Bam"


	def doc(self):
		return (
			"A binary file compressed in the BGZF format with a '.bam' file extension."
		)






class BamInputSorted(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".qname_input_sorted.bam"
		)


	@staticmethod
	def name():
		return "BamInputSorted"


	def doc(self):
		return (
			"A binary file compressed in the BGZF format with a '.bam' file extension and sorted based on the aligner output."
		)






class BamNative(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".unsorted.bam"
		)


	@staticmethod
	def name():
		return "BamNative"


	def doc(self):
		return (
			"A binary file compressed in the BGZF format with a '.bam' file extension."
		)






class BamQuerynameSorted(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".qname_sorted.bam"
		)


	@staticmethod
	def name():
		return "BamQuerynameSorted"


	def doc(self):
		return (
			"A binary file compressed in the BGZF format with a '.bam' file extension and sorted by queryname."
		)






class Bcf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bcf"
		)


	@staticmethod
	def name():
		return "Bcf"


	def doc(self):
		return (
			""
		)






class BcfUncompressed(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bcf_uncompressed"
		)


	@staticmethod
	def name():
		return "BcfUncompressed"


	def doc(self):
		return (
			""
		)






class Bcf_bgzip(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bcf_bgzip"
		)


	@staticmethod
	def name():
		return "Bcf_bgzip"


	def doc(self):
		return (
			""
		)






class Bed(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bed"
		)


	@staticmethod
	def name():
		return "Bed"


	def doc(self):
		return (
			"BED format provides a flexible way to define the data lines that are displayed in an annotation track. BED lines have three required columns and nine additional optional columns. The three required columns are chrom, chromStart and chromEnd."
		)






class Bed12(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bed12"
		)


	@staticmethod
	def name():
		return "Bed12"


	def doc(self):
		return (
			""
		)






class Bed6(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bed6"
		)


	@staticmethod
	def name():
		return "Bed6"


	def doc(self):
		return (
			""
		)






class BedGraph(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bedgraph"
		)


	@staticmethod
	def name():
		return "BedGraph"


	def doc(self):
		return (
			""
		)






class BedStrict(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bedstrict"
		)


	@staticmethod
	def name():
		return "BedStrict"


	def doc(self):
		return (
			""
		)






class Bgzip(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bgzip"
		)


	@staticmethod
	def name():
		return "Bgzip"


	def doc(self):
		return (
			""
		)






class Bif(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bif"
		)


	@staticmethod
	def name():
		return "Bif"


	def doc(self):
		return (
			""
		)






class BigBed(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bigbed"
		)


	@staticmethod
	def name():
		return "BigBed"


	def doc(self):
		return (
			""
		)






class BigWig(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bigwig"
		)


	@staticmethod
	def name():
		return "BigWig"


	def doc(self):
		return (
			""
		)






class Binary(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tpr", alternate_extensions={".binary"}
		)


	@staticmethod
	def name():
		return "Binary"


	def doc(self):
		return (
			""
		)






class Biom1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".biom1"
		)


	@staticmethod
	def name():
		return "Biom1"


	def doc(self):
		return (
			""
		)






class Biom2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".biom2"
		)


	@staticmethod
	def name():
		return "Biom2"


	def doc(self):
		return (
			""
		)






class BlastDomainDb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".blastdbd"
		)


	@staticmethod
	def name():
		return "BlastDomainDb"


	def doc(self):
		return (
			""
		)






class BlastDomainDb5(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".blastdbd5"
		)


	@staticmethod
	def name():
		return "BlastDomainDb5"


	def doc(self):
		return (
			""
		)






class BlastNucDb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".blastdbn"
		)


	@staticmethod
	def name():
		return "BlastNucDb"


	def doc(self):
		return (
			""
		)






class BlastNucDb5(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".blastdbn5"
		)


	@staticmethod
	def name():
		return "BlastNucDb5"


	def doc(self):
		return (
			""
		)






class BlastProtDb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".blastdbp"
		)


	@staticmethod
	def name():
		return "BlastProtDb"


	def doc(self):
		return (
			""
		)






class BlastProtDb5(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".blastdbp5"
		)


	@staticmethod
	def name():
		return "BlastProtDb5"


	def doc(self):
		return (
			""
		)






class BlastXml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".blastxml"
		)


	@staticmethod
	def name():
		return "BlastXml"


	def doc(self):
		return (
			""
		)






class BlibSQlite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".blib"
		)


	@staticmethod
	def name():
		return "BlibSQlite"


	def doc(self):
		return (
			""
		)






class Bmp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bmp"
		)


	@staticmethod
	def name():
		return "Bmp"


	def doc(self):
		return (
			""
		)






class BowtieBaseIndex(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bowtie_base_index"
		)


	@staticmethod
	def name():
		return "BowtieBaseIndex"


	def doc(self):
		return (
			""
		)






class BowtieColorIndex(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bowtie_color_index"
		)


	@staticmethod
	def name():
		return "BowtieColorIndex"


	def doc(self):
		return (
			""
		)






class Bref3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bref3"
		)


	@staticmethod
	def name():
		return "Bref3"


	def doc(self):
		return (
			"Bref3 format is a binary format for storing phased, non-missing genotypes for a list of samples. More information in https://faculty.washington.edu/browning/beagle/bref3.14May18.pdf"
		)






class Btf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".btf"
		)


	@staticmethod
	def name():
		return "Btf"


	def doc(self):
		return (
			""
		)






class Btwisted(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".btwisted"
		)


	@staticmethod
	def name():
		return "Btwisted"


	def doc(self):
		return (
			""
		)






class Bus(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".bus"
		)


	@staticmethod
	def name():
		return "Bus"


	def doc(self):
		return (
			""
		)






class CMAP(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cmap"
		)


	@staticmethod
	def name():
		return "CMAP"


	def doc(self):
		return (
			"The Bionano Genomics cmap format provides location information for label sites within a genome map or an in silico digestion of a reference or sequence data. A CMAP file contains two sections, header and the map information block"
		)






class CML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cml"
		)


	@staticmethod
	def name():
		return "CML"


	def doc(self):
		return (
			""
		)






class CRAM(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cram"
		)


	@staticmethod
	def name():
		return "CRAM"


	def doc(self):
		return (
			"CRAM is a file format for highly efficient and tunable reference-based compression of alignment data."
		)






class CSV(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".csv"
		)


	@staticmethod
	def name():
		return "CSV"


	def doc(self):
		return (
			""
		)






class Cai(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cai"
		)


	@staticmethod
	def name():
		return "Cai"


	def doc(self):
		return (
			""
		)






class Cat_db(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cat_db"
		)


	@staticmethod
	def name():
		return "Cat_db"


	def doc(self):
		return (
			""
		)






class Cel(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cel"
		)


	@staticmethod
	def name():
		return "Cel"


	def doc(self):
		return (
			""
		)






class Charge(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".charge"
		)


	@staticmethod
	def name():
		return "Charge"


	def doc(self):
		return (
			""
		)






class Checktrans(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".checktrans"
		)


	@staticmethod
	def name():
		return "Checktrans"


	def doc(self):
		return (
			""
		)






class Chips(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".chips"
		)


	@staticmethod
	def name():
		return "Chips"


	def doc(self):
		return (
			""
		)






class ChiraSQLite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".chira.sqlite"
		)


	@staticmethod
	def name():
		return "ChiraSQLite"


	def doc(self):
		return (
			""
		)






class ChromInfo(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".len"
		)


	@staticmethod
	def name():
		return "ChromInfo"


	def doc(self):
		return (
			""
		)






class ChromatinInteractions(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".chrint"
		)


	@staticmethod
	def name():
		return "ChromatinInteractions"


	def doc(self):
		return (
			""
		)






class CisML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cisml"
		)


	@staticmethod
	def name():
		return "CisML"


	def doc(self):
		return (
			""
		)






class Ckpt(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ckpt"
		)


	@staticmethod
	def name():
		return "Ckpt"


	def doc(self):
		return (
			""
		)






class Clustal(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".clustal"
		)


	@staticmethod
	def name():
		return "Clustal"


	def doc(self):
		return (
			""
		)






class Codata(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".codata"
		)


	@staticmethod
	def name():
		return "Codata"


	def doc(self):
		return (
			""
		)






class Codcmp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".codcmp"
		)


	@staticmethod
	def name():
		return "Codcmp"


	def doc(self):
		return (
			""
		)






class Coderet(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".coderet"
		)


	@staticmethod
	def name():
		return "Coderet"


	def doc(self):
		return (
			""
		)






class CompressedZipArchive(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".zip"
		)


	@staticmethod
	def name():
		return "CompressedZipArchive"


	def doc(self):
		return (
			""
		)






class Compseq(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".compseq"
		)


	@staticmethod
	def name():
		return "Compseq"


	def doc(self):
		return (
			""
		)






class ConnectivityTable(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ct"
		)


	@staticmethod
	def name():
		return "ConnectivityTable"


	def doc(self):
		return (
			""
		)






class ConsensusTaxonomy(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.cons.taxonomy"
		)


	@staticmethod
	def name():
		return "ConsensusTaxonomy"


	def doc(self):
		return (
			""
		)






class ConsensusXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".consensusxml"
		)


	@staticmethod
	def name():
		return "ConsensusXML"


	def doc(self):
		return (
			""
		)






class Cool(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cool"
		)


	@staticmethod
	def name():
		return "Cool"


	def doc(self):
		return (
			""
		)






class CountTable(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.count_table"
		)


	@staticmethod
	def name():
		return "CountTable"


	def doc(self):
		return (
			""
		)






class Cpgplot(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cpgplot"
		)


	@staticmethod
	def name():
		return "Cpgplot"


	def doc(self):
		return (
			""
		)






class Cpgreport(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cpgreport"
		)


	@staticmethod
	def name():
		return "Cpgreport"


	def doc(self):
		return (
			""
		)






class Cps(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cps"
		)


	@staticmethod
	def name():
		return "Cps"


	def doc(self):
		return (
			""
		)






class Cpt(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cpt"
		)


	@staticmethod
	def name():
		return "Cpt"


	def doc(self):
		return (
			""
		)






class CuffDiffSQlite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cuffdiff.sqlite"
		)


	@staticmethod
	def name():
		return "CuffDiffSQlite"


	def doc(self):
		return (
			""
		)






class Cusp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cusp"
		)


	@staticmethod
	def name():
		return "Cusp"


	def doc(self):
		return (
			""
		)






class CustomTrack(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".customtrack"
		)


	@staticmethod
	def name():
		return "CustomTrack"


	def doc(self):
		return (
			""
		)






class Cut(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cut"
		)


	@staticmethod
	def name():
		return "Cut"


	def doc(self):
		return (
			""
		)






class Cxb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cxb"
		)


	@staticmethod
	def name():
		return "Cxb"


	def doc(self):
		return (
			"Cuffquant output format"
		)






class D3_hierarchy(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".d3_hierarchy"
		)


	@staticmethod
	def name():
		return "D3_hierarchy"


	def doc(self):
		return (
			""
		)






class DAA(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".daa"
		)


	@staticmethod
	def name():
		return "DAA"


	def doc(self):
		return (
			""
		)






class DMND(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dmnd"
		)


	@staticmethod
	def name():
		return "DMND"


	def doc(self):
		return (
			""
		)






class DRF(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".drf"
		)


	@staticmethod
	def name():
		return "DRF"


	def doc(self):
		return (
			""
		)






class Dada2_dada(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dada2_dada"
		)


	@staticmethod
	def name():
		return "Dada2_dada"


	def doc(self):
		return (
			""
		)






class Dada2_errorrates(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dada2_errorrates"
		)


	@staticmethod
	def name():
		return "Dada2_errorrates"


	def doc(self):
		return (
			""
		)






class Dada2_mergepairs(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dada2_mergepairs"
		)


	@staticmethod
	def name():
		return "Dada2_mergepairs"


	def doc(self):
		return (
			""
		)






class Dada2_sequencetable(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dada2_sequencetable"
		)


	@staticmethod
	def name():
		return "Dada2_sequencetable"


	def doc(self):
		return (
			""
		)






class Dada2_uniques(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dada2_uniques"
		)


	@staticmethod
	def name():
		return "Dada2_uniques"


	def doc(self):
		return (
			""
		)






class Dan(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dan"
		)


	@staticmethod
	def name():
		return "Dan"


	def doc(self):
		return (
			""
		)






class Data(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".data"
		)


	@staticmethod
	def name():
		return "Data"


	def doc(self):
		return (
			""
		)






class DataIn(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".linkage_datain"
		)


	@staticmethod
	def name():
		return "DataIn"


	def doc(self):
		return (
			""
		)






class Data_manager_json(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".data_manager_json"
		)


	@staticmethod
	def name():
		return "Data_manager_json"


	def doc(self):
		return (
			""
		)






class Dbmotif(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dbmotif"
		)


	@staticmethod
	def name():
		return "Dbmotif"


	def doc(self):
		return (
			""
		)






class Dbnsfp_tabular(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dbnsfp.tabular"
		)


	@staticmethod
	def name():
		return "Dbnsfp_tabular"


	def doc(self):
		return (
			""
		)






class Dcd(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dcd"
		)


	@staticmethod
	def name():
		return "Dcd"


	def doc(self):
		return (
			""
		)






class Deeptools_compute_matrix_archive(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".deeptools_compute_matrix_archive"
		)


	@staticmethod
	def name():
		return "Deeptools_compute_matrix_archive"


	def doc(self):
		return (
			""
		)






class Deeptools_coverage_matrix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".deeptools_coverage_matrix"
		)


	@staticmethod
	def name():
		return "Deeptools_coverage_matrix"


	def doc(self):
		return (
			""
		)






class Diffseq(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".diffseq"
		)


	@staticmethod
	def name():
		return "Diffseq"


	def doc(self):
		return (
			""
		)






class Digest(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".digest"
		)


	@staticmethod
	def name():
		return "Digest"


	def doc(self):
		return (
			""
		)






class Directory(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".directory"
		)


	@staticmethod
	def name():
		return "Directory"


	def doc(self):
		return (
			""
		)






class DistanceMatrix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.dist"
		)


	@staticmethod
	def name():
		return "DistanceMatrix"


	def doc(self):
		return (
			""
		)






class DlibSQlite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dlib"
		)


	@staticmethod
	def name():
		return "DlibSQlite"


	def doc(self):
		return (
			""
		)






class DotBracket(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dbn"
		)


	@staticmethod
	def name():
		return "DotBracket"


	def doc(self):
		return (
			"Dot-Bracket format is a text-based format for storing both an RNA sequence and its corresponding 2D structure."
		)






class Dreg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dreg"
		)


	@staticmethod
	def name():
		return "Dreg"


	def doc(self):
		return (
			""
		)






class Dta(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dta"
		)


	@staticmethod
	def name():
		return "Dta"


	def doc(self):
		return (
			""
		)






class Dta2d(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dta2d"
		)


	@staticmethod
	def name():
		return "Dta2d"


	def doc(self):
		return (
			""
		)






class Dzi(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".dzi"
		)


	@staticmethod
	def name():
		return "Dzi"


	def doc(self):
		return (
			""
		)






class ENCODEPeak(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".encodepeak"
		)


	@staticmethod
	def name():
		return "ENCODEPeak"


	def doc(self):
		return (
			""
		)






class Edr(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".edr"
		)


	@staticmethod
	def name():
		return "Edr"


	def doc(self):
		return (
			""
		)






class Edta(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".edta"
		)


	@staticmethod
	def name():
		return "Edta"


	def doc(self):
		return (
			""
		)






class Eigenstratgeno(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".eigenstratgeno"
		)


	@staticmethod
	def name():
		return "Eigenstratgeno"


	def doc(self):
		return (
			""
		)






class Eigenstratpca(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".eigenstratpca"
		)


	@staticmethod
	def name():
		return "Eigenstratpca"


	def doc(self):
		return (
			""
		)






class Einverted(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".einverted"
		)


	@staticmethod
	def name():
		return "Einverted"


	def doc(self):
		return (
			""
		)






class Eland(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".eland"
		)


	@staticmethod
	def name():
		return "Eland"


	def doc(self):
		return (
			""
		)






class ElandMulti(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".elandmulti"
		)


	@staticmethod
	def name():
		return "ElandMulti"


	def doc(self):
		return (
			""
		)






class ElibSQlite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".elib"
		)


	@staticmethod
	def name():
		return "ElibSQlite"


	def doc(self):
		return (
			""
		)






class Embl(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".embl"
		)


	@staticmethod
	def name():
		return "Embl"


	def doc(self):
		return (
			""
		)






class Epestfind(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".epestfind"
		)


	@staticmethod
	def name():
		return "Epestfind"


	def doc(self):
		return (
			""
		)






class Eps(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".eps"
		)


	@staticmethod
	def name():
		return "Eps"


	def doc(self):
		return (
			""
		)






class Equicktandem(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".equicktandem"
		)


	@staticmethod
	def name():
		return "Equicktandem"


	def doc(self):
		return (
			""
		)






class Eset(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".eset"
		)


	@staticmethod
	def name():
		return "Eset"


	def doc(self):
		return (
			""
		)






class Est2genome(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".est2genome"
		)


	@staticmethod
	def name():
		return "Est2genome"


	def doc(self):
		return (
			""
		)






class Etandem(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".etandem"
		)


	@staticmethod
	def name():
		return "Etandem"


	def doc(self):
		return (
			""
		)






class Excel(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".excel"
		)


	@staticmethod
	def name():
		return "Excel"


	def doc(self):
		return (
			""
		)






class ExcelXls(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xls", alternate_extensions={".excel.xls"}
		)


	@staticmethod
	def name():
		return "ExcelXls"


	def doc(self):
		return (
			""
		)






class ExpressionJson(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".expression.json"
		)


	@staticmethod
	def name():
		return "ExpressionJson"


	def doc(self):
		return (
			""
		)






class FCS(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fcs"
		)


	@staticmethod
	def name():
		return "FCS"


	def doc(self):
		return (
			"A FCS binary sequence file with a '.fcs' file extension."
		)






class FPS(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fps"
		)


	@staticmethod
	def name():
		return "FPS"


	def doc(self):
		return (
			""
		)






class Fai(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fai"
		)


	@staticmethod
	def name():
		return "Fai"


	def doc(self):
		return (
			"A Fasta Index File is a text file consisting of lines each with five TAB-delimited columns : Name, Length, offset, linebases, Linewidth"
		)






class Fast5Archive(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fast5.tar"
		)


	@staticmethod
	def name():
		return "Fast5Archive"


	def doc(self):
		return (
			""
		)






class Fast5ArchiveBz2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fast5.tar.bz2"
		)


	@staticmethod
	def name():
		return "Fast5ArchiveBz2"


	def doc(self):
		return (
			""
		)






class Fast5ArchiveGz(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fast5.tar.gz"
		)


	@staticmethod
	def name():
		return "Fast5ArchiveGz"


	def doc(self):
		return (
			""
		)






class Fasta(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fasta", alternate_extensions={".fasta.gz"}
		)


	@staticmethod
	def name():
		return "Fasta"


	def doc(self):
		return (
			"A sequence in FASTA format consists of a single-line description, followed by lines of sequence data. The first character of the description line is a greater-than ('>') symbol in the first column. All lines should be shorter than 80 characters."
		)






class Fastg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fastg"
		)


	@staticmethod
	def name():
		return "Fastg"


	def doc(self):
		return (
			"fastg format faithfully represents genome assemblies in the face of allelic polymorphism and assembly uncertainty"
		)






class Fastq(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fastq", alternate_extensions={".fastq.gz", ".fastq.bz2.gz", ".fastq.bz2"}
		)


	@staticmethod
	def name():
		return "Fastq"


	def doc(self):
		return (
			"FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores."
		)






class FastqCSSanger(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fastqcssanger", alternate_extensions={".fastqcssanger.bz2", ".fastqcssanger.bz2.gz", ".fastqcssanger.gz"}
		)


	@staticmethod
	def name():
		return "FastqCSSanger"


	def doc(self):
		return (
			"sequence in in color space phred scored quality values 0:93 represented by ASCII 33:126"
		)






class FastqIllumina(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fastqillumina", alternate_extensions={".fastqillumina.bz2", ".fastqillumina.gz", ".fastqillumina.bz2.gz"}
		)


	@staticmethod
	def name():
		return "FastqIllumina"


	def doc(self):
		return (
			"Sanger variant of the FASTQ format: phred+64"
		)






class FastqSanger(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fastqsanger", alternate_extensions={".fastqsanger.bz2", ".fastqsanger.bz2.gz", ".fastqsanger.gz"}
		)


	@staticmethod
	def name():
		return "FastqSanger"


	def doc(self):
		return (
			"Sanger variant of the FASTQ format: phred+33"
		)






class FastqSolexa(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fastqsolexa", alternate_extensions={".fastqsolexa.gz", ".fastqsolexa.bz2", ".fastqsolexa.bz2.gz"}
		)


	@staticmethod
	def name():
		return "FastqSolexa"


	def doc(self):
		return (
			"Solexa variant of the FASTQ format: solexa+64"
		)






class Feattable(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".feattable"
		)


	@staticmethod
	def name():
		return "Feattable"


	def doc(self):
		return (
			""
		)






class FeatureLocationIndex(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fli"
		)


	@staticmethod
	def name():
		return "FeatureLocationIndex"


	def doc(self):
		return (
			""
		)






class FeatureXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".featurexml"
		)


	@staticmethod
	def name():
		return "FeatureXML"


	def doc(self):
		return (
			""
		)






class Ffdata(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ffdata"
		)


	@staticmethod
	def name():
		return "Ffdata"


	def doc(self):
		return (
			""
		)






class Ffindex(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ffindex"
		)


	@staticmethod
	def name():
		return "Ffindex"


	def doc(self):
		return (
			""
		)






class Fitch(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fitch"
		)


	@staticmethod
	def name():
		return "Fitch"


	def doc(self):
		return (
			""
		)






class Flowclr(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowclr"
		)


	@staticmethod
	def name():
		return "Flowclr"


	def doc(self):
		return (
			"A Flow text file containing population information with a '.flowclr' file extension."
		)






class Flowframe(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowframe"
		)


	@staticmethod
	def name():
		return "Flowframe"


	def doc(self):
		return (
			"Data saved from a R session containing just a flowFrame object"
		)






class Flowmfi(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowmfi"
		)


	@staticmethod
	def name():
		return "Flowmfi"


	def doc(self):
		return (
			"A Flow MFI file with a '.flowmfi' file extension."
		)






class Flowscore(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowscore"
		)


	@staticmethod
	def name():
		return "Flowscore"


	def doc(self):
		return (
			"A Flow Score file with a '.flowscore' file extension."
		)






class Flowset(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowset"
		)


	@staticmethod
	def name():
		return "Flowset"


	def doc(self):
		return (
			"Data saved from a R session containing just a flowSet object"
		)






class Flowstat1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowstat1"
		)


	@staticmethod
	def name():
		return "Flowstat1"


	def doc(self):
		return (
			"A Flow Stats file with a '.flowstat1' file extension."
		)






class Flowstat2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowstat2"
		)


	@staticmethod
	def name():
		return "Flowstat2"


	def doc(self):
		return (
			"A Flow Stats file with a '.flowstat2' file extension."
		)






class Flowstat3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowstat3"
		)


	@staticmethod
	def name():
		return "Flowstat3"


	def doc(self):
		return (
			"A Flow Stats file with a '.flowstat3' file extension."
		)






class Flowtext(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flowtext"
		)


	@staticmethod
	def name():
		return "Flowtext"


	def doc(self):
		return (
			"A Flow text file with a '.flowtext' file extension."
		)






class Flv(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".flv"
		)


	@staticmethod
	def name():
		return "Flv"


	def doc(self):
		return (
			""
		)






class Fped(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fped"
		)


	@staticmethod
	def name():
		return "Fped"


	def doc(self):
		return (
			""
		)






class Fphe(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fphe"
		)


	@staticmethod
	def name():
		return "Fphe"


	def doc(self):
		return (
			""
		)






class Freak(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".freak"
		)


	@staticmethod
	def name():
		return "Freak"


	def doc(self):
		return (
			""
		)






class Frequency(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.freq"
		)


	@staticmethod
	def name():
		return "Frequency"


	def doc(self):
		return (
			""
		)






class Fsom(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fsom"
		)


	@staticmethod
	def name():
		return "Fsom"


	def doc(self):
		return (
			"Data saved from a R session containing just a fSOM object"
		)






class Fuzznuc(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fuzznuc"
		)


	@staticmethod
	def name():
		return "Fuzznuc"


	def doc(self):
		return (
			""
		)






class Fuzzpro(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fuzzpro"
		)


	@staticmethod
	def name():
		return "Fuzzpro"


	def doc(self):
		return (
			""
		)






class Fuzztran(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fuzztran"
		)


	@staticmethod
	def name():
		return "Fuzztran"


	def doc(self):
		return (
			""
		)






class GAFASQLite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gafa.sqlite"
		)


	@staticmethod
	def name():
		return "GAFASQLite"


	def doc(self):
		return (
			""
		)






class Gal(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gal"
		)


	@staticmethod
	def name():
		return "Gal"


	def doc(self):
		return (
			""
		)






class Garnier(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".garnier"
		)


	@staticmethod
	def name():
		return "Garnier"


	def doc(self):
		return (
			""
		)






class Gatk_dbsnp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gatk_dbsnp"
		)


	@staticmethod
	def name():
		return "Gatk_dbsnp"


	def doc(self):
		return (
			""
		)






class Gatk_interval(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gatk_interval"
		)


	@staticmethod
	def name():
		return "Gatk_interval"


	def doc(self):
		return (
			""
		)






class Gatk_recal(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gatk_recal"
		)


	@staticmethod
	def name():
		return "Gatk_recal"


	def doc(self):
		return (
			""
		)






class Gatk_report(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gatk_report"
		)


	@staticmethod
	def name():
		return "Gatk_report"


	def doc(self):
		return (
			""
		)






class Gatk_tranche(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gatk_tranche"
		)


	@staticmethod
	def name():
		return "Gatk_tranche"


	def doc(self):
		return (
			""
		)






class Gcg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gcg"
		)


	@staticmethod
	def name():
		return "Gcg"


	def doc(self):
		return (
			""
		)






class Geecee(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".geecee"
		)


	@staticmethod
	def name():
		return "Geecee"


	def doc(self):
		return (
			""
		)






class GeminiSQLite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gemini.sqlite"
		)


	@staticmethod
	def name():
		return "GeminiSQLite"


	def doc(self):
		return (
			""
		)






class Genbank(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".genbank", alternate_extensions={".genbank.gz"}
		)


	@staticmethod
	def name():
		return "Genbank"


	def doc(self):
		return (
			""
		)






class GeneTrack(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".genetrack"
		)


	@staticmethod
	def name():
		return "GeneTrack"


	def doc(self):
		return (
			""
		)






class GenericAsn1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".asn1"
		)


	@staticmethod
	def name():
		return "GenericAsn1"


	def doc(self):
		return (
			""
		)






class GenericAsn1Binary(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".asn1-binary"
		)


	@staticmethod
	def name():
		return "GenericAsn1Binary"


	def doc(self):
		return (
			""
		)






class GenericXml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xml"
		)


	@staticmethod
	def name():
		return "GenericXml"


	def doc(self):
		return (
			""
		)






class GenomeGraphs(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gg"
		)


	@staticmethod
	def name():
		return "GenomeGraphs"


	def doc(self):
		return (
			""
		)






class GenotypeMatrix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".alohomora_gts"
		)


	@staticmethod
	def name():
		return "GenotypeMatrix"


	def doc(self):
		return (
			""
		)






class GeoJson(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".geojson"
		)


	@staticmethod
	def name():
		return "GeoJson"


	def doc(self):
		return (
			""
		)






class Gfa1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gfa1"
		)


	@staticmethod
	def name():
		return "Gfa1"


	def doc(self):
		return (
			""
		)






class Gfa2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gfa2"
		)


	@staticmethod
	def name():
		return "Gfa2"


	def doc(self):
		return (
			""
		)






class Gff(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gff"
		)


	@staticmethod
	def name():
		return "Gff"


	def doc(self):
		return (
			"GFF lines have nine required fields that must be tab-separated."
		)






class Gff3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gff3", alternate_extensions={".gff3.gz", ".gff3.bz2", ".gff3.bz2.gz"}
		)


	@staticmethod
	def name():
		return "Gff3"


	def doc(self):
		return (
			"The GFF3 format addresses the most common extensions to GFF, while preserving backward compatibility with previous formats."
		)






class Gif(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gif"
		)


	@staticmethod
	def name():
		return "Gif"


	def doc(self):
		return (
			""
		)






class Gifti(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gii", alternate_extensions={".gii.gz"}
		)


	@staticmethod
	def name():
		return "Gifti"


	def doc(self):
		return (
			""
		)






class Gmaj(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gmaj.zip"
		)


	@staticmethod
	def name():
		return "Gmaj"


	def doc(self):
		return (
			""
		)






class Gpr(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gpr"
		)


	@staticmethod
	def name():
		return "Gpr"


	def doc(self):
		return (
			""
		)






class Graph_dot(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".graph_dot"
		)


	@staticmethod
	def name():
		return "Graph_dot"


	def doc(self):
		return (
			""
		)






class Gro(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gro"
		)


	@staticmethod
	def name():
		return "Gro"


	def doc(self):
		return (
			""
		)






class Group(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.groups"
		)


	@staticmethod
	def name():
		return "Group"


	def doc(self):
		return (
			""
		)






class GroupAbund(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.shared"
		)


	@staticmethod
	def name():
		return "GroupAbund"


	def doc(self):
		return (
			""
		)






class Gtf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".gtf"
		)


	@staticmethod
	def name():
		return "Gtf"


	def doc(self):
		return (
			""
		)






class H5(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".h5"
		)


	@staticmethod
	def name():
		return "H5"


	def doc(self):
		return (
			""
		)






class H5MLM(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".h5mlm"
		)


	@staticmethod
	def name():
		return "H5MLM"


	def doc(self):
		return (
			""
		)






class HDT(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hdt"
		)


	@staticmethod
	def name():
		return "HDT"


	def doc(self):
		return (
			""
		)






class Hamamatsu(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vms"
		)


	@staticmethod
	def name():
		return "Hamamatsu"


	def doc(self):
		return (
			""
		)






class Hardklor(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hardklor"
		)


	@staticmethod
	def name():
		return "Hardklor"


	def doc(self):
		return (
			""
		)






class Helixturnhelix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".helixturnhelix"
		)


	@staticmethod
	def name():
		return "Helixturnhelix"


	def doc(self):
		return (
			""
		)






class Hennig86(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hennig86"
		)


	@staticmethod
	def name():
		return "Hennig86"


	def doc(self):
		return (
			""
		)






class Hep_root(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hep.root"
		)


	@staticmethod
	def name():
		return "Hep_root"


	def doc(self):
		return (
			"ROOT binary file."
		)






class HexrdEtaOmeNpz(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hexrd.eta_ome.npz"
		)


	@staticmethod
	def name():
		return "HexrdEtaOmeNpz"


	def doc(self):
		return (
			""
		)






class HexrdImagesNpz(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hexrd.images.npz"
		)


	@staticmethod
	def name():
		return "HexrdImagesNpz"


	def doc(self):
		return (
			""
		)






class HexrdMaterials(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hexrd.materials.h5"
		)


	@staticmethod
	def name():
		return "HexrdMaterials"


	def doc(self):
		return (
			""
		)






class Hhr(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hhr"
		)


	@staticmethod
	def name():
		return "Hhr"


	def doc(self):
		return (
			""
		)






class Hivtrace(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hivtrace"
		)


	@staticmethod
	def name():
		return "Hivtrace"


	def doc(self):
		return (
			""
		)






class Hmmer2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hmm2"
		)


	@staticmethod
	def name():
		return "Hmmer2"


	def doc(self):
		return (
			""
		)






class Hmmer3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hmm3"
		)


	@staticmethod
	def name():
		return "Hmmer3"


	def doc(self):
		return (
			""
		)






class Hmoment(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hmoment"
		)


	@staticmethod
	def name():
		return "Hmoment"


	def doc(self):
		return (
			""
		)






class Html(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".html"
		)


	@staticmethod
	def name():
		return "Html"


	def doc(self):
		return (
			""
		)






class Hyphy_results_json(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hyphy_results.json"
		)


	@staticmethod
	def name():
		return "Hyphy_results_json"


	def doc(self):
		return (
			""
		)






class ICM(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".icm"
		)


	@staticmethod
	def name():
		return "ICM"


	def doc(self):
		return (
			""
		)






class IQTree(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".iqtree"
		)


	@staticmethod
	def name():
		return "IQTree"


	def doc(self):
		return (
			""
		)






class IdXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".idxml"
		)


	@staticmethod
	def name():
		return "IdXML"


	def doc(self):
		return (
			""
		)






class Idat(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".idat"
		)


	@staticmethod
	def name():
		return "Idat"


	def doc(self):
		return (
			""
		)






class IdeasPre(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ideaspre"
		)


	@staticmethod
	def name():
		return "IdeasPre"


	def doc(self):
		return (
			""
		)






class IdpDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".idpdb"
		)


	@staticmethod
	def name():
		return "IdpDB"


	def doc(self):
		return (
			""
		)






class Ig(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ig"
		)


	@staticmethod
	def name():
		return "Ig"


	def doc(self):
		return (
			""
		)






class Im(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".im"
		)


	@staticmethod
	def name():
		return "Im"


	def doc(self):
		return (
			""
		)






class ImgtJson(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".imgt.json"
		)


	@staticmethod
	def name():
		return "ImgtJson"


	def doc(self):
		return (
			""
		)






class ImzML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".imzml"
		)


	@staticmethod
	def name():
		return "ImzML"


	def doc(self):
		return (
			""
		)






class InChI(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".inchi"
		)


	@staticmethod
	def name():
		return "InChI"


	def doc(self):
		return (
			""
		)






class InfernalCM(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".cm"
		)


	@staticmethod
	def name():
		return "InfernalCM"


	def doc(self):
		return (
			""
		)






class Intermine_tabular(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".intermine_tabular"
		)


	@staticmethod
	def name():
		return "Intermine_tabular"


	def doc(self):
		return (
			""
		)






class Interprophet_pepxml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".interprophet_pepxml"
		)


	@staticmethod
	def name():
		return "Interprophet_pepxml"


	def doc(self):
		return (
			""
		)






class Interval(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".interval"
		)


	@staticmethod
	def name():
		return "Interval"


	def doc(self):
		return (
			"File must start with definition line in the following format (columns may be in any order)."
		)






class Interval_index(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".interval_index"
		)


	@staticmethod
	def name():
		return "Interval_index"


	def doc(self):
		return (
			""
		)






class Ipynb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ipynb"
		)


	@staticmethod
	def name():
		return "Ipynb"


	def doc(self):
		return (
			""
		)






class IsaJson(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".isa-json"
		)


	@staticmethod
	def name():
		return "IsaJson"


	def doc(self):
		return (
			"ISA-JSON data type."
		)






class IsaTab(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".isa-tab"
		)


	@staticmethod
	def name():
		return "IsaTab"


	def doc(self):
		return (
			"ISA-Tab data type."
		)






class Isochore(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".isochore"
		)


	@staticmethod
	def name():
		return "Isochore"


	def doc(self):
		return (
			""
		)






class Itp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".itp"
		)


	@staticmethod
	def name():
		return "Itp"


	def doc(self):
		return (
			""
		)






class JP2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".jp2"
		)


	@staticmethod
	def name():
		return "JP2"


	def doc(self):
		return (
			""
		)






class Jackknifer(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".jackknifer"
		)


	@staticmethod
	def name():
		return "Jackknifer"


	def doc(self):
		return (
			""
		)






class Jackknifernon(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".jackknifernon"
		)


	@staticmethod
	def name():
		return "Jackknifernon"


	def doc(self):
		return (
			""
		)






class Jellyfish(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".jellyfish"
		)


	@staticmethod
	def name():
		return "Jellyfish"


	def doc(self):
		return (
			"Jellyfish database files are k-mer counts in binary format with a readable head. They are operated on and converted to human-readable text through jellyfish commands."
		)






class Jpg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".jpg"
		)


	@staticmethod
	def name():
		return "Jpg"


	def doc(self):
		return (
			""
		)






class Json(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".json"
		)


	@staticmethod
	def name():
		return "Json"


	def doc(self):
		return (
			""
		)






class Jsonld(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".jsonld"
		)


	@staticmethod
	def name():
		return "Jsonld"


	def doc(self):
		return (
			""
		)






class Kallisto_idx(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".kallisto.idx"
		)


	@staticmethod
	def name():
		return "Kallisto_idx"


	def doc(self):
		return (
			""
		)






class Kroenik(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".kroenik"
		)


	@staticmethod
	def name():
		return "Kroenik"


	def doc(self):
		return (
			""
		)






class Kronik(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".kronik"
		)


	@staticmethod
	def name():
		return "Kronik"


	def doc(self):
		return (
			""
		)






class Laj(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".laj"
		)


	@staticmethod
	def name():
		return "Laj"


	def doc(self):
		return (
			""
		)






class LaneMask(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.filter"
		)


	@staticmethod
	def name():
		return "LaneMask"


	def doc(self):
		return (
			""
		)






class LastDb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".lastdb"
		)


	@staticmethod
	def name():
		return "LastDb"


	def doc(self):
		return (
			""
		)






class Lav(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".lav"
		)


	@staticmethod
	def name():
		return "Lav"


	def doc(self):
		return (
			"Lav is the primary output format for BLASTZ.  The first line of a .lav file begins with #:lav.."
		)






class LineCount(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".linecount"
		)


	@staticmethod
	def name():
		return "LineCount"


	def doc(self):
		return (
			""
		)






class Linkage_pedin(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".linkage_pedin"
		)


	@staticmethod
	def name():
		return "Linkage_pedin"


	def doc(self):
		return (
			""
		)






class Loom(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".loom"
		)


	@staticmethod
	def name():
		return "Loom"


	def doc(self):
		return (
			"An HDF5-based Loom File"
		)






class LowerTriangleDistanceMatrix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.lower.dist"
		)


	@staticmethod
	def name():
		return "LowerTriangleDistanceMatrix"


	def doc(self):
		return (
			""
		)






class Lped(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".lped"
		)


	@staticmethod
	def name():
		return "Lped"


	def doc(self):
		return (
			""
		)






class MAlist(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".malist"
		)


	@staticmethod
	def name():
		return "MAlist"


	def doc(self):
		return (
			""
		)






class MCool(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mcool"
		)


	@staticmethod
	def name():
		return "MCool"


	def doc(self):
		return (
			""
		)






class MEMEXml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".memexml"
		)


	@staticmethod
	def name():
		return "MEMEXml"


	def doc(self):
		return (
			""
		)






class MOL(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mol"
		)


	@staticmethod
	def name():
		return "MOL"


	def doc(self):
		return (
			""
		)






class MOL2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mol2"
		)


	@staticmethod
	def name():
		return "MOL2"


	def doc(self):
		return (
			""
		)






class Maf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".maf"
		)


	@staticmethod
	def name():
		return "Maf"


	def doc(self):
		return (
			"TBA and multiz multiple alignment format.  The first line of a .maf file begins with ##maf. This word is followed by white-space-separated 'variable=value' pairs. There should be no white space surrounding the '='."
		)






class MafCustomTrack(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mafcustomtrack"
		)


	@staticmethod
	def name():
		return "MafCustomTrack"


	def doc(self):
		return (
			""
		)






class MarkerMap(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".linkage_map"
		)


	@staticmethod
	def name():
		return "MarkerMap"


	def doc(self):
		return (
			""
		)






class Markx0(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".markx0"
		)


	@staticmethod
	def name():
		return "Markx0"


	def doc(self):
		return (
			""
		)






class Markx1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".markx1"
		)


	@staticmethod
	def name():
		return "Markx1"


	def doc(self):
		return (
			""
		)






class Markx10(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".markx10"
		)


	@staticmethod
	def name():
		return "Markx10"


	def doc(self):
		return (
			""
		)






class Markx2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".markx2"
		)


	@staticmethod
	def name():
		return "Markx2"


	def doc(self):
		return (
			""
		)






class Markx3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".markx3"
		)


	@staticmethod
	def name():
		return "Markx3"


	def doc(self):
		return (
			""
		)






class MascotDat(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mascotdat"
		)


	@staticmethod
	def name():
		return "MascotDat"


	def doc(self):
		return (
			""
		)






class MascotXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mascotxml"
		)


	@staticmethod
	def name():
		return "MascotXML"


	def doc(self):
		return (
			""
		)






class MashSketch(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".msh"
		)


	@staticmethod
	def name():
		return "MashSketch"


	def doc(self):
		return (
			""
		)






class Maskinfo_asn1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".maskinfo-asn1"
		)


	@staticmethod
	def name():
		return "Maskinfo_asn1"


	def doc(self):
		return (
			""
		)






class Maskinfo_asn1_binary(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".maskinfo-asn1-binary"
		)


	@staticmethod
	def name():
		return "Maskinfo_asn1_binary"


	def doc(self):
		return (
			""
		)






class MassHunterTar(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".agilentmasshunter.d.tar"
		)


	@staticmethod
	def name():
		return "MassHunterTar"


	def doc(self):
		return (
			""
		)






class MassLynxTar(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".watersmasslynx.raw.tar"
		)


	@staticmethod
	def name():
		return "MassLynxTar"


	def doc(self):
		return (
			""
		)






class Match(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".match"
		)


	@staticmethod
	def name():
		return "Match"


	def doc(self):
		return (
			""
		)






class MatrixMarket(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mtx"
		)


	@staticmethod
	def name():
		return "MatrixMarket"


	def doc(self):
		return (
			""
		)






class MauveXmfa(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xmfa"
		)


	@staticmethod
	def name():
		return "MauveXmfa"


	def doc(self):
		return (
			""
		)






class Mdp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mdp"
		)


	@staticmethod
	def name():
		return "Mdp"


	def doc(self):
		return (
			""
		)






class Mega(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mega"
		)


	@staticmethod
	def name():
		return "Mega"


	def doc(self):
		return (
			""
		)






class Meganon(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".meganon"
		)


	@staticmethod
	def name():
		return "Meganon"


	def doc(self):
		return (
			""
		)






class MemePsp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".memepsp"
		)


	@staticmethod
	def name():
		return "MemePsp"


	def doc(self):
		return (
			"The MEME Position Specific Priors (PSP) format includes the name of the sequence for which a prior distribution corresponds."
		)






class Meryldb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".meryldb"
		)


	@staticmethod
	def name():
		return "Meryldb"


	def doc(self):
		return (
			"MerylDB is a tar.gz archive containing 64 binaries + 64 indexes."
		)






class Metacyto_clr_txt(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".metacyto_clr.txt"
		)


	@staticmethod
	def name():
		return "Metacyto_clr_txt"


	def doc(self):
		return (
			"List of clusters used in MetaCyto analyses"
		)






class Mgf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mgf"
		)


	@staticmethod
	def name():
		return "Mgf"


	def doc(self):
		return (
			""
		)






class Mirax(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mrxs"
		)


	@staticmethod
	def name():
		return "Mirax"


	def doc(self):
		return (
			""
		)






class Mkv(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mkv"
		)


	@staticmethod
	def name():
		return "Mkv"


	def doc(self):
		return (
			""
		)






class Mothur_align(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.align"
		)


	@staticmethod
	def name():
		return "Mothur_align"


	def doc(self):
		return (
			""
		)






class Mothur_design(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.design"
		)


	@staticmethod
	def name():
		return "Mothur_design"


	def doc(self):
		return (
			""
		)






class Mothur_filtered_masked_quan(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.filtered.masked.quan"
		)


	@staticmethod
	def name():
		return "Mothur_filtered_masked_quan"


	def doc(self):
		return (
			""
		)






class Mothur_filtered_quan(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.filtered.quan"
		)


	@staticmethod
	def name():
		return "Mothur_filtered_quan"


	def doc(self):
		return (
			""
		)






class Mothur_list(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.list"
		)


	@staticmethod
	def name():
		return "Mothur_list"


	def doc(self):
		return (
			""
		)






class Mothur_masked_quan(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.masked.quan"
		)


	@staticmethod
	def name():
		return "Mothur_masked_quan"


	def doc(self):
		return (
			""
		)






class Mothur_otu_corr(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.otu.corr"
		)


	@staticmethod
	def name():
		return "Mothur_otu_corr"


	def doc(self):
		return (
			""
		)






class Mothur_rabund(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.rabund"
		)


	@staticmethod
	def name():
		return "Mothur_rabund"


	def doc(self):
		return (
			""
		)






class Mothur_rdp_taxonomy(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.rdp.taxonomy"
		)


	@staticmethod
	def name():
		return "Mothur_rdp_taxonomy"


	def doc(self):
		return (
			""
		)






class Mothur_relabund(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.relabund"
		)


	@staticmethod
	def name():
		return "Mothur_relabund"


	def doc(self):
		return (
			""
		)






class Mothur_seq_taxonomy(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.seq.taxonomy"
		)


	@staticmethod
	def name():
		return "Mothur_seq_taxonomy"


	def doc(self):
		return (
			""
		)






class Mothur_tre(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.tre"
		)


	@staticmethod
	def name():
		return "Mothur_tre"


	def doc(self):
		return (
			""
		)






class Motif(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".motif"
		)


	@staticmethod
	def name():
		return "Motif"


	def doc(self):
		return (
			""
		)






class Mp3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mp3"
		)


	@staticmethod
	def name():
		return "Mp3"


	def doc(self):
		return (
			""
		)






class Mp4(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mp4"
		)


	@staticmethod
	def name():
		return "Mp4"


	def doc(self):
		return (
			""
		)






class Mpg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mpg"
		)


	@staticmethod
	def name():
		return "Mpg"


	def doc(self):
		return (
			""
		)






class Mrc2014(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mrc"
		)


	@staticmethod
	def name():
		return "Mrc2014"


	def doc(self):
		return (
			""
		)






class Mrm(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mrm"
		)


	@staticmethod
	def name():
		return "Mrm"


	def doc(self):
		return (
			""
		)






class Ms2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ms2"
		)


	@staticmethod
	def name():
		return "Ms2"


	def doc(self):
		return (
			""
		)






class Msf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".msf"
		)


	@staticmethod
	def name():
		return "Msf"


	def doc(self):
		return (
			""
		)






class Msp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".msp"
		)


	@staticmethod
	def name():
		return "Msp"


	def doc(self):
		return (
			""
		)






class Mz5(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mz5"
		)


	@staticmethod
	def name():
		return "Mz5"


	def doc(self):
		return (
			""
		)






class MzData(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mzdata"
		)


	@staticmethod
	def name():
		return "MzData"


	def doc(self):
		return (
			""
		)






class MzIdentML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mzid"
		)


	@staticmethod
	def name():
		return "MzIdentML"


	def doc(self):
		return (
			""
		)






class MzML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mzml"
		)


	@staticmethod
	def name():
		return "MzML"


	def doc(self):
		return (
			""
		)






class MzQuantML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mzq"
		)


	@staticmethod
	def name():
		return "MzQuantML"


	def doc(self):
		return (
			""
		)






class MzSQlite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mz.sqlite"
		)


	@staticmethod
	def name():
		return "MzSQlite"


	def doc(self):
		return (
			""
		)






class MzTab(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mztab"
		)


	@staticmethod
	def name():
		return "MzTab"


	def doc(self):
		return (
			""
		)






class MzTab2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mztab2"
		)


	@staticmethod
	def name():
		return "MzTab2"


	def doc(self):
		return (
			""
		)






class MzXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mzxml"
		)


	@staticmethod
	def name():
		return "MzXML"


	def doc(self):
		return (
			""
		)






class N3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".n3"
		)


	@staticmethod
	def name():
		return "N3"


	def doc(self):
		return (
			""
		)






class NTriples(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nt"
		)


	@staticmethod
	def name():
		return "NTriples"


	def doc(self):
		return (
			""
		)






class Names(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.names"
		)


	@staticmethod
	def name():
		return "Names"


	def doc(self):
		return (
			""
		)






class Nametable(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nametable"
		)


	@staticmethod
	def name():
		return "Nametable"


	def doc(self):
		return (
			""
		)






class Ncbi(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ncbi"
		)


	@staticmethod
	def name():
		return "Ncbi"


	def doc(self):
		return (
			""
		)






class NcbiTaxonomySQlite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ncbitaxonomy.sqlite"
		)


	@staticmethod
	def name():
		return "NcbiTaxonomySQlite"


	def doc(self):
		return (
			""
		)






class Ndpi(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ndpi"
		)


	@staticmethod
	def name():
		return "Ndpi"


	def doc(self):
		return (
			""
		)






class Ndx(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ndx"
		)


	@staticmethod
	def name():
		return "Ndx"


	def doc(self):
		return (
			""
		)






class Needle(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".needle"
		)


	@staticmethod
	def name():
		return "Needle"


	def doc(self):
		return (
			""
		)






class Neo4jDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".neostore"
		)


	@staticmethod
	def name():
		return "Neo4jDB"


	def doc(self):
		return (
			""
		)






class Neo4jDBzip(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".neostore.zip"
		)


	@staticmethod
	def name():
		return "Neo4jDBzip"


	def doc(self):
		return (
			""
		)






class NetCDF(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".netcdf"
		)


	@staticmethod
	def name():
		return "NetCDF"


	def doc(self):
		return (
			"Format used by netCDF software library for writing and reading chromatography-MS data files."
		)






class Newcpgreport(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".newcpgreport"
		)


	@staticmethod
	def name():
		return "Newcpgreport"


	def doc(self):
		return (
			""
		)






class Newcpgseek(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".newcpgseek"
		)


	@staticmethod
	def name():
		return "Newcpgseek"


	def doc(self):
		return (
			""
		)






class Newick(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".newick"
		)


	@staticmethod
	def name():
		return "Newick"


	def doc(self):
		return (
			""
		)






class Nexus(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nex", alternate_extensions={".nexus"}
		)


	@staticmethod
	def name():
		return "Nexus"


	def doc(self):
		return (
			""
		)






class Nexusnon(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nexusnon"
		)


	@staticmethod
	def name():
		return "Nexusnon"


	def doc(self):
		return (
			""
		)






class Nhdr(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nhdr"
		)


	@staticmethod
	def name():
		return "Nhdr"


	def doc(self):
		return (
			""
		)






class Nhx(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nhx"
		)


	@staticmethod
	def name():
		return "Nhx"


	def doc(self):
		return (
			""
		)






class Nifti1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nii1", alternate_extensions={".nii1.gz"}
		)


	@staticmethod
	def name():
		return "Nifti1"


	def doc(self):
		return (
			""
		)






class Nifti2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nii2", alternate_extensions={".nii2.gz"}
		)


	@staticmethod
	def name():
		return "Nifti2"


	def doc(self):
		return (
			""
		)






class NmrML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nmrml"
		)


	@staticmethod
	def name():
		return "NmrML"


	def doc(self):
		return (
			"nmrML is an open mark-up language for NMR data."
		)






class Noreturn(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".noreturn"
		)


	@staticmethod
	def name():
		return "Noreturn"


	def doc(self):
		return (
			""
		)






class Npz(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".npz"
		)


	@staticmethod
	def name():
		return "Npz"


	def doc(self):
		return (
			""
		)






class Nrrd(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".nrrd"
		)


	@staticmethod
	def name():
		return "Nrrd"


	def doc(self):
		return (
			""
		)






class OBFS(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".obfs"
		)


	@staticmethod
	def name():
		return "OBFS"


	def doc(self):
		return (
			""
		)






class OMETiff(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ome.tiff"
		)


	@staticmethod
	def name():
		return "OMETiff"


	def doc(self):
		return (
			""
		)






class OSW(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".osw"
		)


	@staticmethod
	def name():
		return "OSW"


	def doc(self):
		return (
			""
		)






class Obo(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".obo"
		)


	@staticmethod
	def name():
		return "Obo"


	def doc(self):
		return (
			""
		)






class Odgi(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".odgi"
		)


	@staticmethod
	def name():
		return "Odgi"


	def doc(self):
		return (
			"Genomic variation graphs self index used by odgi."
		)






class Oligos(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.oligos"
		)


	@staticmethod
	def name():
		return "Oligos"


	def doc(self):
		return (
			""
		)






class Onnx(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".onnx"
		)


	@staticmethod
	def name():
		return "Onnx"


	def doc(self):
		return (
			"ONNX (Open neural network exchange) is data format for storing and sharing machine learning and deep learning models."
		)






class Otu(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.otu"
		)


	@staticmethod
	def name():
		return "Otu"


	def doc(self):
		return (
			""
		)






class Owl(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".owl"
		)


	@staticmethod
	def name():
		return "Owl"


	def doc(self):
		return (
			""
		)






class OxliCountGraph(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".oxlicg"
		)


	@staticmethod
	def name():
		return "OxliCountGraph"


	def doc(self):
		return (
			""
		)






class OxliGraphLabels(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".oxligl"
		)


	@staticmethod
	def name():
		return "OxliGraphLabels"


	def doc(self):
		return (
			""
		)






class OxliNodeGraph(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".oxling"
		)


	@staticmethod
	def name():
		return "OxliNodeGraph"


	def doc(self):
		return (
			""
		)






class OxliStopTags(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".oxlist"
		)


	@staticmethod
	def name():
		return "OxliStopTags"


	def doc(self):
		return (
			""
		)






class OxliSubset(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".oxliss"
		)


	@staticmethod
	def name():
		return "OxliSubset"


	def doc(self):
		return (
			""
		)






class OxliTagSet(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".oxlits"
		)


	@staticmethod
	def name():
		return "OxliTagSet"


	def doc(self):
		return (
			""
		)






class PDB(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pdb"
		)


	@staticmethod
	def name():
		return "PDB"


	def doc(self):
		return (
			""
		)






class PDBQT(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pdbqt"
		)


	@staticmethod
	def name():
		return "PDBQT"


	def doc(self):
		return (
			""
		)






class PEFF(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".peff"
		)


	@staticmethod
	def name():
		return "PEFF"


	def doc(self):
		return (
			""
		)






class PHAR(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".phar"
		)


	@staticmethod
	def name():
		return "PHAR"


	def doc(self):
		return (
			""
		)






class PQP(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pqp"
		)


	@staticmethod
	def name():
		return "PQP"


	def doc(self):
		return (
			""
		)






class PQR(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pqr"
		)


	@staticmethod
	def name():
		return "PQR"


	def doc(self):
		return (
			""
		)






class PSMS(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".psms"
		)


	@staticmethod
	def name():
		return "PSMS"


	def doc(self):
		return (
			""
		)






class Paf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".paf", alternate_extensions={".paf.gz"}
		)


	@staticmethod
	def name():
		return "Paf"


	def doc(self):
		return (
			""
		)






class Pair(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pair"
		)


	@staticmethod
	def name():
		return "Pair"


	def doc(self):
		return (
			""
		)






class PairwiseDistanceMatrix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.pair.dist"
		)


	@staticmethod
	def name():
		return "PairwiseDistanceMatrix"


	def doc(self):
		return (
			""
		)






class Palindrome(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".palindrome"
		)


	@staticmethod
	def name():
		return "Palindrome"


	def doc(self):
		return (
			""
		)






class Paramxml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".paramxml"
		)


	@staticmethod
	def name():
		return "Paramxml"


	def doc(self):
		return (
			""
		)






class Parquet(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".parquet"
		)


	@staticmethod
	def name():
		return "Parquet"


	def doc(self):
		return (
			""
		)






class Pbed(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pbed"
		)


	@staticmethod
	def name():
		return "Pbed"


	def doc(self):
		return (
			""
		)






class Pbm(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pbm"
		)


	@staticmethod
	def name():
		return "Pbm"


	def doc(self):
		return (
			""
		)






class Pcd(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pcd"
		)


	@staticmethod
	def name():
		return "Pcd"


	def doc(self):
		return (
			""
		)






class Pcx(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pcx"
		)


	@staticmethod
	def name():
		return "Pcx"


	def doc(self):
		return (
			""
		)






class Pdf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pdf"
		)


	@staticmethod
	def name():
		return "Pdf"


	def doc(self):
		return (
			""
		)






class PepList(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".peplist"
		)


	@staticmethod
	def name():
		return "PepList"


	def doc(self):
		return (
			""
		)






class PepXml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pepxml"
		)


	@staticmethod
	def name():
		return "PepXml"


	def doc(self):
		return (
			""
		)






class PepXmlReport(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pepxml.tsv"
		)


	@staticmethod
	def name():
		return "PepXmlReport"


	def doc(self):
		return (
			""
		)






class Pepcoil(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pepcoil"
		)


	@staticmethod
	def name():
		return "Pepcoil"


	def doc(self):
		return (
			""
		)






class Pepinfo(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pepinfo"
		)


	@staticmethod
	def name():
		return "Pepinfo"


	def doc(self):
		return (
			""
		)






class Pepstats(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pepstats"
		)


	@staticmethod
	def name():
		return "Pepstats"


	def doc(self):
		return (
			""
		)






class Peptideprophet_pepxml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".peptideprophet_pepxml"
		)


	@staticmethod
	def name():
		return "Peptideprophet_pepxml"


	def doc(self):
		return (
			""
		)






class Peptideshaker_archive(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".peptideshaker_archive"
		)


	@staticmethod
	def name():
		return "Peptideshaker_archive"


	def doc(self):
		return (
			""
		)






class Percin(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".percin"
		)


	@staticmethod
	def name():
		return "Percin"


	def doc(self):
		return (
			""
		)






class Percout(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".percout"
		)


	@staticmethod
	def name():
		return "Percout"


	def doc(self):
		return (
			""
		)






class Pgm(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pgm"
		)


	@staticmethod
	def name():
		return "Pgm"


	def doc(self):
		return (
			""
		)






class Pheno(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pheno"
		)


	@staticmethod
	def name():
		return "Pheno"


	def doc(self):
		return (
			""
		)






class Phylip(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".phylip"
		)


	@staticmethod
	def name():
		return "Phylip"


	def doc(self):
		return (
			""
		)






class Phylipnon(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".phylipnon"
		)


	@staticmethod
	def name():
		return "Phylipnon"


	def doc(self):
		return (
			""
		)






class Phyloxml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".phyloxml"
		)


	@staticmethod
	def name():
		return "Phyloxml"


	def doc(self):
		return (
			""
		)






class Picard_interval_list(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".picard_interval_list"
		)


	@staticmethod
	def name():
		return "Picard_interval_list"


	def doc(self):
		return (
			""
		)






class Pileup(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pileup"
		)


	@staticmethod
	def name():
		return "Pileup"


	def doc(self):
		return (
			""
		)






class Pir(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pir"
		)


	@staticmethod
	def name():
		return "Pir"


	def doc(self):
		return (
			""
		)






class PlantTribesKsComponents(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ptkscmp"
		)


	@staticmethod
	def name():
		return "PlantTribesKsComponents"


	def doc(self):
		return (
			""
		)






class PlyAscii(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".plyascii"
		)


	@staticmethod
	def name():
		return "PlyAscii"


	def doc(self):
		return (
			""
		)






class PlyBinary(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".plybinary"
		)


	@staticmethod
	def name():
		return "PlyBinary"


	def doc(self):
		return (
			""
		)






class Png(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".png"
		)


	@staticmethod
	def name():
		return "Png"


	def doc(self):
		return (
			""
		)






class Polydot(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".polydot"
		)


	@staticmethod
	def name():
		return "Polydot"


	def doc(self):
		return (
			""
		)






class PostgresqlArchive(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".postgresql"
		)


	@staticmethod
	def name():
		return "PostgresqlArchive"


	def doc(self):
		return (
			""
		)






class Pphe(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pphe"
		)


	@staticmethod
	def name():
		return "Pphe"


	def doc(self):
		return (
			""
		)






class Ppm(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ppm"
		)


	@staticmethod
	def name():
		return "Ppm"


	def doc(self):
		return (
			""
		)






class Preg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".preg"
		)


	@staticmethod
	def name():
		return "Preg"


	def doc(self):
		return (
			""
		)






class Pretext(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pretext"
		)


	@staticmethod
	def name():
		return "Pretext"


	def doc(self):
		return (
			""
		)






class Prettyseq(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".prettyseq"
		)


	@staticmethod
	def name():
		return "Prettyseq"


	def doc(self):
		return (
			""
		)






class Primersearch(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".primersearch"
		)


	@staticmethod
	def name():
		return "Primersearch"


	def doc(self):
		return (
			""
		)






class ProBam(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".probam"
		)


	@staticmethod
	def name():
		return "ProBam"


	def doc(self):
		return (
			""
		)






class ProBed(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".probed"
		)


	@staticmethod
	def name():
		return "ProBed"


	def doc(self):
		return (
			""
		)






class ProtXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".protxml"
		)


	@staticmethod
	def name():
		return "ProtXML"


	def doc(self):
		return (
			""
		)






class ProtXmlReport(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".protxml.tsv"
		)


	@staticmethod
	def name():
		return "ProtXmlReport"


	def doc(self):
		return (
			""
		)






class Protobuf2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".protobuf2"
		)


	@staticmethod
	def name():
		return "Protobuf2"


	def doc(self):
		return (
			"Protocol Buffers (Protobuf) is data format for serializing structured data."
		)






class Protobuf3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".protobuf3"
		)


	@staticmethod
	def name():
		return "Protobuf3"


	def doc(self):
		return (
			"Protocol Buffers (Protobuf) is data format for serializing structured data."
		)






class Psd(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".psd"
		)


	@staticmethod
	def name():
		return "Psd"


	def doc(self):
		return (
			""
		)






class Pssm_asn1(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".pssm-asn1"
		)


	@staticmethod
	def name():
		return "Pssm_asn1"


	def doc(self):
		return (
			""
		)






class QCML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".qcml"
		)


	@staticmethod
	def name():
		return "QCML"


	def doc(self):
		return (
			"Quality control data in XML format (https://code.google.com/p/qcml/)."
		)






class QualityScore(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".qual"
		)


	@staticmethod
	def name():
		return "QualityScore"


	def doc(self):
		return (
			""
		)






class QualityScore454(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".qual454"
		)


	@staticmethod
	def name():
		return "QualityScore454"


	def doc(self):
		return (
			""
		)






class QualityScoreIllumina(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".qualillumina"
		)


	@staticmethod
	def name():
		return "QualityScoreIllumina"


	def doc(self):
		return (
			""
		)






class QualityScoreSOLiD(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".qualsolid"
		)


	@staticmethod
	def name():
		return "QualityScoreSOLiD"


	def doc(self):
		return (
			""
		)






class QualityScoreSolexa(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".qualsolexa"
		)


	@staticmethod
	def name():
		return "QualityScoreSolexa"


	def doc(self):
		return (
			""
		)






class Quantile(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.quan"
		)


	@staticmethod
	def name():
		return "Quantile"


	def doc(self):
		return (
			""
		)






class RData(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata"
		)


	@staticmethod
	def name():
		return "RData"


	def doc(self):
		return (
			"Stored data from an R session"
		)






class RMA6(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rma6"
		)


	@staticmethod
	def name():
		return "RMA6"


	def doc(self):
		return (
			""
		)






class RNADotPlotMatrix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rna_eps"
		)


	@staticmethod
	def name():
		return "RNADotPlotMatrix"


	def doc(self):
		return (
			""
		)






class Rast(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rast"
		)


	@staticmethod
	def name():
		return "Rast"


	def doc(self):
		return (
			""
		)






class Raw_pepxml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".raw_pepxml"
		)


	@staticmethod
	def name():
		return "Raw_pepxml"


	def doc(self):
		return (
			""
		)






class Rdata_camera_negative(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.camera.negative"
		)


	@staticmethod
	def name():
		return "Rdata_camera_negative"


	def doc(self):
		return (
			""
		)






class Rdata_camera_positive(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.camera.positive"
		)


	@staticmethod
	def name():
		return "Rdata_camera_positive"


	def doc(self):
		return (
			""
		)






class Rdata_camera_quick(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.camera.quick"
		)


	@staticmethod
	def name():
		return "Rdata_camera_quick"


	def doc(self):
		return (
			""
		)






class Rdata_eset(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.eset"
		)


	@staticmethod
	def name():
		return "Rdata_eset"


	def doc(self):
		return (
			"Stored RDS from a ExpressionSet Object"
		)






class Rdata_msnbase_raw(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.msnbase.raw"
		)


	@staticmethod
	def name():
		return "Rdata_msnbase_raw"


	def doc(self):
		return (
			""
		)






class Rdata_sce(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.sce"
		)


	@staticmethod
	def name():
		return "Rdata_sce"


	def doc(self):
		return (
			"Stored RDS from a SingleCellObject"
		)






class Rdata_xcms_fillpeaks(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.xcms.fillpeaks"
		)


	@staticmethod
	def name():
		return "Rdata_xcms_fillpeaks"


	def doc(self):
		return (
			""
		)






class Rdata_xcms_findchrompeaks(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.xcms.findchrompeaks"
		)


	@staticmethod
	def name():
		return "Rdata_xcms_findchrompeaks"


	def doc(self):
		return (
			""
		)






class Rdata_xcms_group(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.xcms.group"
		)


	@staticmethod
	def name():
		return "Rdata_xcms_group"


	def doc(self):
		return (
			""
		)






class Rdata_xcms_raw(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.xcms.raw"
		)


	@staticmethod
	def name():
		return "Rdata_xcms_raw"


	def doc(self):
		return (
			""
		)






class Rdata_xcms_retcor(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdata.xcms.retcor"
		)


	@staticmethod
	def name():
		return "Rdata_xcms_retcor"


	def doc(self):
		return (
			""
		)






class Rdf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdf"
		)


	@staticmethod
	def name():
		return "Rdf"


	def doc(self):
		return (
			""
		)






class Rdock_as(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rdock_as"
		)


	@staticmethod
	def name():
		return "Rdock_as"


	def doc(self):
		return (
			"rDock active site format"
		)






class RefTaxonomy(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.ref.taxonomy"
		)


	@staticmethod
	def name():
		return "RefTaxonomy"


	def doc(self):
		return (
			""
		)






class Regions(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".regions"
		)


	@staticmethod
	def name():
		return "Regions"


	def doc(self):
		return (
			""
		)






class RexpBase(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rexpbase"
		)


	@staticmethod
	def name():
		return "RexpBase"


	def doc(self):
		return (
			""
		)






class Rgb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rgb"
		)


	@staticmethod
	def name():
		return "Rgb"


	def doc(self):
		return (
			""
		)






class Rgenetics(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".rgenetics"
		)


	@staticmethod
	def name():
		return "Rgenetics"


	def doc(self):
		return (
			""
		)






class Roadmaps(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".roadmaps"
		)


	@staticmethod
	def name():
		return "Roadmaps"


	def doc(self):
		return (
			""
		)






class SDF(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sdf"
		)


	@staticmethod
	def name():
		return "SDF"


	def doc(self):
		return (
			""
		)






class SMILES(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".smi"
		)


	@staticmethod
	def name():
		return "SMILES"


	def doc(self):
		return (
			""
		)






class SNPMatrix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".snpmatrix"
		)


	@staticmethod
	def name():
		return "SNPMatrix"


	def doc(self):
		return (
			""
		)






class SPLib(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".splib"
		)


	@staticmethod
	def name():
		return "SPLib"


	def doc(self):
		return (
			""
		)






class SPLibNoIndex(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".splib_noindex"
		)


	@staticmethod
	def name():
		return "SPLibNoIndex"


	def doc(self):
		return (
			""
		)






class SQlite(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sqlite"
		)


	@staticmethod
	def name():
		return "SQlite"


	def doc(self):
		return (
			""
		)






class SQmass(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sqmass"
		)


	@staticmethod
	def name():
		return "SQmass"


	def doc(self):
		return (
			""
		)






class STL(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".stl"
		)


	@staticmethod
	def name():
		return "STL"


	def doc(self):
		return (
			""
		)






class Sabund(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.sabund"
		)


	@staticmethod
	def name():
		return "Sabund"


	def doc(self):
		return (
			""
		)






class Sakura(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".svslide"
		)


	@staticmethod
	def name():
		return "Sakura"


	def doc(self):
		return (
			""
		)






class Sam(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sam"
		)


	@staticmethod
	def name():
		return "Sam"


	def doc(self):
		return (
			""
		)






class Sbml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sbml"
		)


	@staticmethod
	def name():
		return "Sbml"


	def doc(self):
		return (
			""
		)






class ScIdx(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".scidx"
		)


	@staticmethod
	def name():
		return "ScIdx"


	def doc(self):
		return (
			""
		)






class Scf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".scf"
		)


	@staticmethod
	def name():
		return "Scf"


	def doc(self):
		return (
			"A binary sequence file in 'scf' format with a '.scf' file extension.  You must manually select this 'File Format' when uploading the file."
		)






class Scn(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".scn"
		)


	@staticmethod
	def name():
		return "Scn"


	def doc(self):
		return (
			""
		)






class Score(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".score"
		)


	@staticmethod
	def name():
		return "Score"


	def doc(self):
		return (
			""
		)






class SearchGuiArchive(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".searchgui_archive"
		)


	@staticmethod
	def name():
		return "SearchGuiArchive"


	def doc(self):
		return (
			""
		)






class SecondaryStructureMap(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.map"
		)


	@staticmethod
	def name():
		return "SecondaryStructureMap"


	def doc(self):
		return (
			""
		)






class Seqtable(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".seqtable"
		)


	@staticmethod
	def name():
		return "Seqtable"


	def doc(self):
		return (
			""
		)






class SequenceSplitLocations(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".fqtoc"
		)


	@staticmethod
	def name():
		return "SequenceSplitLocations"


	def doc(self):
		return (
			""
		)






class Sequences(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sequences"
		)


	@staticmethod
	def name():
		return "Sequences"


	def doc(self):
		return (
			""
		)






class Sf3(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sf3"
		)


	@staticmethod
	def name():
		return "Sf3"


	def doc(self):
		return (
			""
		)






class Sff(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sff"
		)


	@staticmethod
	def name():
		return "Sff"


	def doc(self):
		return (
			"A binary file in 'Standard Flowgram Format' with a '.sff' file extension."
		)






class SffFlow(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.sff.flow"
		)


	@staticmethod
	def name():
		return "SffFlow"


	def doc(self):
		return (
			""
		)






class Shapefile(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".shp"
		)


	@staticmethod
	def name():
		return "Shapefile"


	def doc(self):
		return (
			""
		)






class Showfeat(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".showfeat"
		)


	@staticmethod
	def name():
		return "Showfeat"


	def doc(self):
		return (
			""
		)






class Showorf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".showorf"
		)


	@staticmethod
	def name():
		return "Showorf"


	def doc(self):
		return (
			""
		)






class Sif(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sif"
		)


	@staticmethod
	def name():
		return "Sif"


	def doc(self):
		return (
			""
		)






class Simple(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".simple"
		)


	@staticmethod
	def name():
		return "Simple"


	def doc(self):
		return (
			""
		)






class Sirius_ms(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sirius.ms"
		)


	@staticmethod
	def name():
		return "Sirius_ms"


	def doc(self):
		return (
			""
		)






class Sixpack(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sixpack"
		)


	@staticmethod
	def name():
		return "Sixpack"


	def doc(self):
		return (
			""
		)






class Smat(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".smat"
		)


	@staticmethod
	def name():
		return "Smat"


	def doc(self):
		return (
			""
		)






class SnapHmm(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".snaphmm"
		)


	@staticmethod
	def name():
		return "SnapHmm"


	def doc(self):
		return (
			""
		)






class SnpEffDb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".snpeffdb"
		)


	@staticmethod
	def name():
		return "SnpEffDb"


	def doc(self):
		return (
			""
		)






class SnpSiftDbNSFP(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".snpsiftdbnsfp"
		)


	@staticmethod
	def name():
		return "SnpSiftDbNSFP"


	def doc(self):
		return (
			""
		)






class Snptest(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".snptest"
		)


	@staticmethod
	def name():
		return "Snptest"


	def doc(self):
		return (
			""
		)






class Source_c(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".source.c"
		)


	@staticmethod
	def name():
		return "Source_c"


	def doc(self):
		return (
			"C source file"
		)






class Source_cpp(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".source.cpp"
		)


	@staticmethod
	def name():
		return "Source_cpp"


	def doc(self):
		return (
			"C++ source file"
		)






class Source_cs(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".source.cs"
		)


	@staticmethod
	def name():
		return "Source_cs"


	def doc(self):
		return (
			"C# source file"
		)






class Source_go(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".source.go"
		)


	@staticmethod
	def name():
		return "Source_go"


	def doc(self):
		return (
			"Go source file"
		)






class Source_h(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".source.h"
		)


	@staticmethod
	def name():
		return "Source_h"


	def doc(self):
		return (
			"C or cpp header file"
		)






class Source_py(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".source.py"
		)


	@staticmethod
	def name():
		return "Source_py"


	def doc(self):
		return (
			"Python source file"
		)






class Source_rs(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".source.rs"
		)


	@staticmethod
	def name():
		return "Source_rs"


	def doc(self):
		return (
			"Rust source file"
		)






class SpalnNuclDb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".spalndbnp"
		)


	@staticmethod
	def name():
		return "SpalnNuclDb"


	def doc(self):
		return (
			""
		)






class SpalnProtDb(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".spalndba"
		)


	@staticmethod
	def name():
		return "SpalnProtDb"


	def doc(self):
		return (
			""
		)






class SquareDistanceMatrix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mldist", alternate_extensions={".mothur.square.dist"}
		)


	@staticmethod
	def name():
		return "SquareDistanceMatrix"


	def doc(self):
		return (
			""
		)






class Sra(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sra"
		)


	@staticmethod
	def name():
		return "Sra"


	def doc(self):
		return (
			"A binary file archive format from the NCBI Sequence Read Archive with a '.sra' file extension."
		)






class SraManifest(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".sra_manifest.tabular"
		)


	@staticmethod
	def name():
		return "SraManifest"


	def doc(self):
		return (
			""
		)






class Srs(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".srs"
		)


	@staticmethod
	def name():
		return "Srs"


	def doc(self):
		return (
			""
		)






class Srspair(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".srspair"
		)


	@staticmethod
	def name():
		return "Srspair"


	def doc(self):
		return (
			""
		)






class Staden(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".staden"
		)


	@staticmethod
	def name():
		return "Staden"


	def doc(self):
		return (
			""
		)






class Star(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".star"
		)


	@staticmethod
	def name():
		return "Star"


	def doc(self):
		return (
			""
		)






class Stockholm_1_0(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".stockholm"
		)


	@staticmethod
	def name():
		return "Stockholm_1_0"


	def doc(self):
		return (
			""
		)






class Strider(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".strider"
		)


	@staticmethod
	def name():
		return "Strider"


	def doc(self):
		return (
			""
		)






class Summary(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.summary"
		)


	@staticmethod
	def name():
		return "Summary"


	def doc(self):
		return (
			""
		)






class Supermatcher(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".supermatcher"
		)


	@staticmethod
	def name():
		return "Supermatcher"


	def doc(self):
		return (
			""
		)






class Svg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".svg"
		)


	@staticmethod
	def name():
		return "Svg"


	def doc(self):
		return (
			""
		)






class Svs(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".svs"
		)


	@staticmethod
	def name():
		return "Svs"


	def doc(self):
		return (
			""
		)






class Swiss(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".swiss"
		)


	@staticmethod
	def name():
		return "Swiss"


	def doc(self):
		return (
			""
		)






class Syco(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".syco"
		)


	@staticmethod
	def name():
		return "Syco"


	def doc(self):
		return (
			""
		)






class TSV(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tsv"
		)


	@staticmethod
	def name():
		return "TSV"


	def doc(self):
		return (
			""
		)






class Tabix(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tabix"
		)


	@staticmethod
	def name():
		return "Tabix"


	def doc(self):
		return (
			""
		)






class Table(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".table"
		)


	@staticmethod
	def name():
		return "Table"


	def doc(self):
		return (
			""
		)






class Tabular(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tabular", alternate_extensions={".allegro_descent.gz", ".allegro_ihaplo.gz", ".allegro_ihaplo", ".allegro_descent", ".tabular.gz"}
		)


	@staticmethod
	def name():
		return "Tabular"


	def doc(self):
		return (
			", Any data in tab delimited format (tabular)."
		)






class Tagseq(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tagseq"
		)


	@staticmethod
	def name():
		return "Tagseq"


	def doc(self):
		return (
			""
		)






class TandemXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tandem"
		)


	@staticmethod
	def name():
		return "TandemXML"


	def doc(self):
		return (
			""
		)






class Tar(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tar"
		)


	@staticmethod
	def name():
		return "Tar"


	def doc(self):
		return (
			""
		)






class Taxonomy(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".taxonomy"
		)


	@staticmethod
	def name():
		return "Taxonomy"


	def doc(self):
		return (
			""
		)






class TaxonomySummary(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".mothur.tax.summary"
		)


	@staticmethod
	def name():
		return "TaxonomySummary"


	def doc(self):
		return (
			""
		)






class Tck(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tck"
		)


	@staticmethod
	def name():
		return "Tck"


	def doc(self):
		return (
			""
		)






class TdfTar(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".brukertdf.d.tar"
		)


	@staticmethod
	def name():
		return "TdfTar"


	def doc(self):
		return (
			""
		)






class Text(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".txt"
		)


	@staticmethod
	def name():
		return "Text"


	def doc(self):
		return (
			"Any text file."
		)






class TextGrid(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".textgrid"
		)


	@staticmethod
	def name():
		return "TextGrid"


	def doc(self):
		return (
			""
		)






class Textsearch(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".textsearch"
		)


	@staticmethod
	def name():
		return "Textsearch"


	def doc(self):
		return (
			""
		)






class Tf2(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tf2"
		)


	@staticmethod
	def name():
		return "Tf2"


	def doc(self):
		return (
			""
		)






class Tf8(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tf8"
		)


	@staticmethod
	def name():
		return "Tf8"


	def doc(self):
		return (
			""
		)






class Tgz(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tgz"
		)


	@staticmethod
	def name():
		return "Tgz"


	def doc(self):
		return (
			""
		)






class ThermoRAW(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".thermo.raw"
		)


	@staticmethod
	def name():
		return "ThermoRAW"


	def doc(self):
		return (
			""
		)






class Tif(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tif"
		)


	@staticmethod
	def name():
		return "Tif"


	def doc(self):
		return (
			""
		)






class Tiff(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".tiff"
		)


	@staticmethod
	def name():
		return "Tiff"


	def doc(self):
		return (
			""
		)






class Toml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".toml"
		)


	@staticmethod
	def name():
		return "Toml"


	def doc(self):
		return (
			""
		)






class Toolshed_gz(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".toolshed.gz"
		)


	@staticmethod
	def name():
		return "Toolshed_gz"


	def doc(self):
		return (
			""
		)






class Top(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".top"
		)


	@staticmethod
	def name():
		return "Top"


	def doc(self):
		return (
			""
		)






class TraML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".traml"
		)


	@staticmethod
	def name():
		return "TraML"


	def doc(self):
		return (
			""
		)






class TrafoXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".trafoxml"
		)


	@staticmethod
	def name():
		return "TrafoXML"


	def doc(self):
		return (
			""
		)






class Triples(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".triples"
		)


	@staticmethod
	def name():
		return "Triples"


	def doc(self):
		return (
			""
		)






class Trk(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".trk"
		)


	@staticmethod
	def name():
		return "Trk"


	def doc(self):
		return (
			""
		)






class Trr(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".trr"
		)


	@staticmethod
	def name():
		return "Trr"


	def doc(self):
		return (
			""
		)






class Turtle(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ttl"
		)


	@staticmethod
	def name():
		return "Turtle"


	def doc(self):
		return (
			""
		)






class TwoBit(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".twobit"
		)


	@staticmethod
	def name():
		return "TwoBit"


	def doc(self):
		return (
			""
		)






class UCSCTrackHub(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".trackhub"
		)


	@staticmethod
	def name():
		return "UCSCTrackHub"


	def doc(self):
		return (
			""
		)






class UniProtXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".uniprotxml"
		)


	@staticmethod
	def name():
		return "UniProtXML"


	def doc(self):
		return (
			""
		)






class Vcf(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vcf"
		)


	@staticmethod
	def name():
		return "Vcf"


	def doc(self):
		return (
			""
		)






class VcfGz(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vcf_bgzip"
		)


	@staticmethod
	def name():
		return "VcfGz"


	def doc(self):
		return (
			""
		)






class Vectorstrip(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vectorstrip"
		)


	@staticmethod
	def name():
		return "Vectorstrip"


	def doc(self):
		return (
			""
		)






class Vel(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vel"
		)


	@staticmethod
	def name():
		return "Vel"


	def doc(self):
		return (
			""
		)






class Velvet(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".velvet"
		)


	@staticmethod
	def name():
		return "Velvet"


	def doc(self):
		return (
			""
		)






class Vg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vg"
		)


	@staticmethod
	def name():
		return "Vg"


	def doc(self):
		return (
			"Genomic variation graphs."
		)






class Vmu(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vmu"
		)


	@staticmethod
	def name():
		return "Vmu"


	def doc(self):
		return (
			""
		)






class VtkAscii(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vtkascii"
		)


	@staticmethod
	def name():
		return "VtkAscii"


	def doc(self):
		return (
			""
		)






class VtkBinary(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".vtkbinary"
		)


	@staticmethod
	def name():
		return "VtkBinary"


	def doc(self):
		return (
			""
		)






class Wav(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".wav"
		)


	@staticmethod
	def name():
		return "Wav"


	def doc(self):
		return (
			""
		)






class Wiff(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".wiff"
		)


	@staticmethod
	def name():
		return "Wiff"


	def doc(self):
		return (
			""
		)






class WiffTar(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".wiff.tar"
		)


	@staticmethod
	def name():
		return "WiffTar"


	def doc(self):
		return (
			""
		)






class Wiggle(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".wig"
		)


	@staticmethod
	def name():
		return "Wiggle"


	def doc(self):
		return (
			"The wiggle format is line-oriented.  Wiggle data is preceded by a track definition line, which adds a number of options for controlling the default display of this track."
		)






class Wobble(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".wobble"
		)


	@staticmethod
	def name():
		return "Wobble"


	def doc(self):
		return (
			""
		)






class Wordcount(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".wordcount"
		)


	@staticmethod
	def name():
		return "Wordcount"


	def doc(self):
		return (
			""
		)






class XHunterAslFormat(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".hlf"
		)


	@staticmethod
	def name():
		return "XHunterAslFormat"


	def doc(self):
		return (
			""
		)






class Xbm(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xbm"
		)


	@staticmethod
	def name():
		return "Xbm"


	def doc(self):
		return (
			""
		)






class Xg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xg"
		)


	@staticmethod
	def name():
		return "Xg"


	def doc(self):
		return (
			"Genomic variation graphs with vg index."
		)






class Xgmml(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xgmml"
		)


	@staticmethod
	def name():
		return "Xgmml"


	def doc(self):
		return (
			""
		)






class Xlsx(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xlsx"
		)


	@staticmethod
	def name():
		return "Xlsx"


	def doc(self):
		return (
			""
		)






class Xpm(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xpm"
		)


	@staticmethod
	def name():
		return "Xpm"


	def doc(self):
		return (
			""
		)






class XquestSpecXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".spec.xml"
		)


	@staticmethod
	def name():
		return "XquestSpecXML"


	def doc(self):
		return (
			""
		)






class XquestXML(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xquest.xml"
		)


	@staticmethod
	def name():
		return "XquestXML"


	def doc(self):
		return (
			""
		)






class Xtc(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xtc"
		)


	@staticmethod
	def name():
		return "Xtc"


	def doc(self):
		return (
			""
		)






class Xvg(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".xvg"
		)


	@staticmethod
	def name():
		return "Xvg"


	def doc(self):
		return (
			""
		)






class YepTar(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".agilentbrukeryep.d.tar"
		)


	@staticmethod
	def name():
		return "YepTar"


	def doc(self):
		return (
			""
		)






class csFasta(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".csfasta"
		)


	@staticmethod
	def name():
		return "csFasta"


	def doc(self):
		return (
			""
		)






class grd(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".grd"
		)


	@staticmethod
	def name():
		return "grd"


	def doc(self):
		return (
			""
		)






class grdtgz(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".grd.tgz"
		)


	@staticmethod
	def name():
		return "grdtgz"


	def doc(self):
		return (
			""
		)






class ldIndep(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".ldindep"
		)


	@staticmethod
	def name():
		return "ldIndep"


	def doc(self):
		return (
			""
		)






class mStats(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".metacyto_stats.txt"
		)


	@staticmethod
	def name():
		return "mStats"


	def doc(self):
		return (
			"Table of statistics generated by a MetaCyto analysis"
		)






class mSummary(File):
	def __init__(self, optional=False):
		super().__init__(
			optional, extension=".metacyto_summary.txt"
		)


	@staticmethod
	def name():
		return "mSummary"


	def doc(self):
		return (
			"Summary table generated by MetaCyto preprocessing of FCS files"
		)


