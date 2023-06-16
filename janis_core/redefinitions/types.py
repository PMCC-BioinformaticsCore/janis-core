
import operator
import os.path
import zipfile
import subprocess

from typing import Any, Dict, List, Optional

from janis_core import File, Array, Logger
from janis_core.tool.test_classes import TTestExpectedOutput, TTestPreprocessor


### HTML ###

class HtmlFile(File):
    def __init__(self, optional=False, extension=".html"):
        super().__init__(optional, extension=extension)

    @staticmethod
    def name():
        return "HtmlFile"

    def doc(self):
        return "A HTML file"



### ZIPFILES ###

class Gunzipped(File):
    def __init__(self, inner_type=File(), optional=False, extension=".gz"):
        super().__init__(optional, extension=extension)
        self.inner_type = inner_type

    def id(self):
        inner = f"Gzipped<{self.inner_type.name()}>"
        if self.optional:
            return f"Optional<{inner}>"
        return inner

    @staticmethod
    def name():
        return "Gzip"

    def doc(self):
        return "A gzipped file"
    

class FileTabix(Gunzipped):
    @staticmethod
    def secondary_files():
        return [".tbi"]


class ZipFile(File):
    def __init__(self, optional=False, extension=".zip"):
        super().__init__(optional, extension=extension)

    @staticmethod
    def name():
        return "Zip"

    def doc(self):
        return "A zip archive, ending with .zip"

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_size: int,
        md5: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_size,
            ),
        ]
        if md5 is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.FileMd5,
                    operator=operator.eq,
                    expected_value=md5,
                ),
            ]
        return outcome
    

### TEXTFILES ###

class TextFile(File):
    def __init__(self, optional=False, extension=".txt"):
        super().__init__(optional, extension=extension)

    @staticmethod
    def name():
        return "TextFile"

    def doc(self):
        return "A textfile, ending with .txt"

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_size: int,
        min_required_content: Optional[str] = None,
        line_count: Optional[int] = None,
        md5: Optional[str] = None,
        expected_file_path: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_size,
            ),
        ]

        if min_required_content is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.FileContent,
                    operator=operator.contains,
                    expected_value=min_required_content,
                ),
            ]

        if line_count is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.LineCount,
                    operator=operator.eq,
                    expected_value=line_count,
                ),
            ]

        if md5 is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.FileMd5,
                    operator=operator.eq,
                    expected_value=md5,
                ),
            ]

        if expected_file_path is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.FileContent,
                    operator=operator.eq,
                    expected_file=expected_file_path,
                ),
            ]

        return outcome


class Tsv(TextFile):
    def __init__(self, optional=False, extension=".tsv"):
        super().__init__(optional, extension=extension)

    @staticmethod
    def name():
        return "tsv"

    def doc(self):
        return "A tab separated file"
    

class Csv(TextFile):
    def __init__(self, optional=False, extension=".csv"):
        super().__init__(optional, extension=extension)

    @staticmethod
    def name():
        return "csv"

    def doc(self):
        return "A comma separated file"



### VCF ###

class Vcf(File):
    def __init__(self, optional=False):
        super().__init__(optional=optional, extension=".vcf")

    @staticmethod
    def name():
        return "VCF"

    def doc(self):
        return """
    Variant Call Format:

    The Variant Call Format (VCF) specifies the format of a text file 
    used in bioinformatics for storing gene sequence variations. 

    Documentation: https://samtools.github.io/hts-specs/VCFv4.3.pdf
    """.strip()

    @classmethod
    def md5_without_header(cls, file_path: str, headers_to_remove: List[str]) -> str:
        """
        Compute md5 of a vcf file with unwanted headers removed

        :param file_path: path to the file
        :type file_path: str
        :param headers_to_remove: headers to be removed before computing md5
        :type headers_to_remove: List[str]
        :return: md5
        :rtype: str
        """
        with open(file_path, "r") as f:
            hash_md5 = hashlib.md5()
            while True:
                line = f.readline()
                if not line:
                    break
                if all(("##" + header) not in line for header in headers_to_remove):
                    hash_md5.update(line.encode())
        return hash_md5.hexdigest()

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_size: int,
        line_count: Optional[int] = None,
        headers_to_remove: Optional[List[str]] = None,
        md5_value: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_size,
            ),
        ]

        if line_count is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.LineCount,
                    operator=operator.eq,
                    expected_value=line_count,
                ),
            ]

        if md5_value is not None:
            if headers_to_remove is not None:
                outcome += [
                    TTestExpectedOutput(
                        tag=tag,
                        preprocessor=Vcf.md5_without_header,
                        operator=operator.eq,
                        expected_value=md5_value,
                        preprocessor_params={"headers_to_remove": headers_to_remove},
                    ),
                ]
            else:
                outcome += [
                    TTestExpectedOutput(
                        tag=tag,
                        preprocessor=TTestPreprocessor.FileMd5,
                        operator=operator.eq,
                        expected_value=md5_value,
                    ),
                ]

        return outcome


class CompressedVcf(Gunzipped):
    def __init__(self, optional=False):
        super().__init__(inner_type=Vcf, optional=optional, extension=".vcf.gz")

    @staticmethod
    def name():
        return "CompressedVCF"

    def doc(self):
        return ".vcf.gz"

    @classmethod
    def LineCount(cls, file_path: str):
        count = 0
        with gzip.open(file_path, "rt") as f:
            while f.readline():
                count += 1
        return count

    @classmethod
    def md5_without_header(cls, file_path: str, headers_to_remove: List[str]) -> str:
        """
        Compute md5 of a vcf.gz file with unwanted headers removed

        :param file_path: path to the file
        :type file_path: str
        :param headers_to_remove: headers to be removed before computing md5
        :type headers_to_remove: List[str]
        :return: md5
        :rtype: str
        """
        with gzip.open(file_path, "rt") as f:
            hash_md5 = hashlib.md5()
            while True:
                line = f.readline()
                if not line:
                    break
                if all(("##" + header) not in line for header in headers_to_remove):
                    hash_md5.update(line.encode())
        return hash_md5.hexdigest()

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_size: int,
        line_count: Optional[int] = None,
        headers_to_remove: Optional[List[str]] = None,
        md5_value: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_size,
            ),
        ]

        if line_count is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=CompressedVcf.LineCount,
                    operator=operator.eq,
                    expected_value=line_count,
                ),
            ]

        if md5_value is not None:
            if headers_to_remove is not None:
                outcome += [
                    TTestExpectedOutput(
                        tag=tag,
                        preprocessor=CompressedVcf.md5_without_header,
                        operator=operator.eq,
                        expected_value=md5_value,
                        preprocessor_params={"headers_to_remove": headers_to_remove},
                    ),
                ]
            else:
                outcome += [
                    TTestExpectedOutput(
                        tag=tag,
                        preprocessor=TTestPreprocessor.FileMd5,
                        operator=operator.eq,
                        expected_value=md5_value,
                    ),
                ]

        return outcome


class VcfTabix(CompressedVcf, FileTabix):
    @staticmethod
    def name():
        return "CompressedIndexedVCF"

    def doc(self):
        return ".vcf.gz with .vcf.gz.tbi file"

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_vcf_size: int,
        min_tbi_size: int,
        line_count: Optional[int] = None,
        headers_to_remove: Optional[List[str]] = None,
        vcf_md5: Optional[str] = None,
        tbi_md5: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = super().basic_test(
            tag, min_vcf_size, line_count, headers_to_remove, vcf_md5
        ) + [
            TTestExpectedOutput(
                tag=tag,
                suffix_secondary_file=".tbi",
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_tbi_size,
            ),
        ]
        if tbi_md5 is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.FileMd5,
                    operator=operator.eq,
                    expected_value=tbi_md5,
                ),
            ]
        return outcome






### FASTA ###

class Fasta(File):
    def __init__(self, optional=False):
        super().__init__(
            optional, extension=".fasta", alternate_extensions={".fa", ".fna"}
        )

    @staticmethod
    def name():
        return "Fasta"

    def can_receive_from(self, other, source_has_default=False):

        if isinstance(other, Fasta):
            if other.optional and not self.optional:
                return False
        elif not super().can_receive_from(other, source_has_default):
            return False

        if not self.secondary_files():
            return True

        return set(self.secondary_files()).issubset(set(other.secondary_files() or []))


class FastaGz(File):
    def __init__(self, optional=False):
        super().__init__(optional, extension=".fa.gz")

    @staticmethod
    def name():
        return "FastaGz"

    def can_receive_from(self, other, source_has_default=False):

        if isinstance(other, FastaGz):
            if other.optional and not self.optional:
                return False
        elif not super().can_receive_from(other, source_has_default):
            return False

        if not self.secondary_files():
            return True

        return set(self.secondary_files()).issubset(set(other.secondary_files() or []))
    

class FastaFai(Fasta):
    @staticmethod
    def name():
        return "FastaFai"

    @staticmethod
    def secondary_files():
        return [".fai"]


class FastaBwa(Fasta):
    @staticmethod
    def name():
        return "FastaBwa"

    @staticmethod
    def secondary_files():
        return [".amb", ".ann", ".bwt", ".pac", ".sa"]


class FastaDict(Fasta):
    @staticmethod
    def name():
        return "FastDict"

    @staticmethod
    def secondary_files():
        return ["^.dict"]


class FastaWithIndexes(Fasta):
    @staticmethod
    def name():
        return "FastaWithIndexes"

    @staticmethod
    def secondary_files():
        return [
            *FastaFai.secondary_files(),
            *FastaBwa.secondary_files(),
            *FastaDict.secondary_files(),
        ]
    
FastaWithDict = FastaWithIndexes


### FASTQ ###


class Fastq(File):
    def __init__(self, optional=False):
        super().__init__(
            optional=optional, extension=".fastq", alternate_extensions={".fq"}
        )

    @staticmethod
    def name():
        return "Fastq"

    def doc(self):
        return (
            "FASTQ files are text files containing sequence data with quality score, there are different types"
            "with no standard: https://www.drive5.com/usearch/manual/fastq_files.html"
        )


class FastqGz(File):
    def __init__(self, optional=False):
        super().__init__(
            optional=optional, extension=".fastq.gz", alternate_extensions={".fq.gz"}
        )

    @staticmethod
    def name():
        return "FastqGz"

    def doc(self):
        return (
            "FastqGz files are compressed sequence data with quality score, there are different types"
            "with no standard: https://en.wikipedia.org/wiki/FASTQ_format"
        )

    @classmethod
    def basic_test(cls, tag: str, min_size: int) -> List[TTestExpectedOutput]:
        return [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_size,
            ),
        ]
    

class FastqPairedEnd(Array):
    def __init__(self, optional=False):
        super().__init__(Fastq(optional=False), optional=optional)

    def id(self):
        if self.optional:
            return f"Optional<{self.name()}>"
        return self.name()

    @staticmethod
    def name():
        return "FastqPair"

    def doc(self):
        return "Paired end Fastq files "

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool):
        if not super().validate_value(meta, allow_null_if_not_optional):
            return False
        return len(meta) == 2

    def invalid_value_hint(self, meta):
        prev = super().invalid_value_hint(meta)
        hints = []
        if prev:
            hints.append(prev)
        if len(meta) != 2:
            hints.append(f"There must be exactly 2 (found {len(meta)}) fastq files")
        return ", ".join(hints)

    def is_paired(self):
        return True
    

class FastqGzPairedEnd(Array):
    def __init__(self, optional=False):
        super().__init__(FastqGz, optional=optional)

    @staticmethod
    def name():
        return "FastqGzPair"

    def id(self):
        if self.optional:
            return f"Optional<{self.name()}>"
        return self.name()

    def doc(self):
        return "Paired end FastqGz files"

    def validate_value(self, meta: Any, allow_null_if_not_optional: bool):
        super_is_valid = super().validate_value(meta, allow_null_if_not_optional)
        if not super_is_valid or meta is None:
            return super_is_valid

        return len(meta) == 2

    def invalid_value_hint(self, meta):
        prev = super().invalid_value_hint(meta)
        hints = []
        if prev:
            hints.append(prev)

        if meta is not None and len(meta) != 2:
            hints.append(f"There must be exactly 2 (found {len(meta)}) fastq files")
        return ", ".join(hints)

    def is_paired(self):
        return True

    @classmethod
    def ge(cls, file_paths: str, expected_sizes: List[int]):
        """

        :param file_paths: a string containing all file paths, separated by |
        :type file_paths: str
        :param expected_sizes: expected minimum sizes of all files
        :type expected_sizes: List[int]
        :return: a boolean value indicating if all files are bigger than or equal to their expected minimum sizes
        """
        files = file_paths.split("|")
        if len(files) != len(expected_sizes):
            return "Number of expected values don't match number of outputs"
        for file in files:
            unzipped = zipfile.ZipFile(file)
            if "R1" in unzipped.namelist()[0]:
                if os.path.getsize(file) < expected_sizes[0]:
                    return False
            else:
                if os.path.getsize(file) < expected_sizes[1]:
                    return False
        return True

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_total_size: int,
        min_first_size: Optional[int] = None,
        min_second_size: Optional[int] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.ListSize,
                operator=operator.eq,
                expected_value=2,
            ),
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.ListOfFilesTotalSize,
                operator=operator.ge,
                expected_value=min_total_size,
            ),
        ]

        # An example of how FastqGzPairedEnd.ge is used; can be removed if deemed unnecessary
        if min_first_size is not None and min_second_size is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.Value,
                    operator=FastqGzPairedEnd.ge,
                    expected_value=[min_first_size, min_second_size],
                )
            ]
        return outcome

FastqGzPair = FastqGzPairedEnd
FastqPair = FastqPairedEnd
    

### BAM ###

class Sam(File):
    def __init__(self, optional=False):
        super().__init__(optional=optional, extension=".sam")

    @staticmethod
    def name():
        return "SAM"

    def doc(self):
        return "Tab-delimited text file that contains sequence alignment data"


class Cram(File):
    def __init__(self, optional=False):
        super().__init__(optional, extension=".cram")

    @staticmethod
    def name():
        return "CRAM"

    def doc(self):
        return "A binary version of a SAM file, https://samtools.github.io/hts-specs/CRAMv3.pdf"


class CramCrai(Cram):
    @staticmethod
    def name():
        return "CramPair"

    @staticmethod
    def secondary_files():
        return [".crai"]

    def doc(self):
        return "A Cram and Crai as the secondary"


class Bam(File):
    def __init__(self, optional=False):
        super().__init__(optional, extension=".bam")

    @staticmethod
    def name():
        return "BAM"

    def doc(self):
        return "A binary version of a SAM file, http://software.broadinstitute.org/software/igv/bam"

    @classmethod
    def flagstat(cls, file_path: str):
        command = ["samtools", "flagstat", file_path]
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        if result.stderr:
            raise Exception(result.stderr)

        return result.stdout

    @classmethod
    def equal(cls, file_path_1: str, file_path_2: str):
        flagstat1 = cls.flagstat(file_path_1)
        flagstat2 = cls.flagstat(file_path_2)

        return flagstat1 == flagstat2

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_bam_size: int,
        flagstat: Optional[str] = None,
        md5: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = [
            TTestExpectedOutput(
                tag=tag,
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_bam_size,
            )
        ]
        if flagstat is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=Bam.flagstat,
                    operator=operator.eq,
                    expected_file=flagstat,
                )
            ]
        if md5 is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    preprocessor=TTestPreprocessor.FileMd5,
                    operator=operator.eq,
                    expected_value=md5,
                ),
            ]
        return outcome


class BamBai(Bam):
    @staticmethod
    def name():
        return "IndexedBam"

    @staticmethod
    def secondary_files():
        return [".bai"]

    def doc(self):
        return "A Bam and bai as the secondary"

    @classmethod
    def basic_test(
        cls,
        tag: str,
        min_bam_size: int,
        min_bai_size: int,
        flagstat: Optional[str] = None,
        bam_md5: Optional[str] = None,
        bai_md5: Optional[str] = None,
    ) -> List[TTestExpectedOutput]:
        outcome = super().basic_test(tag, min_bam_size, flagstat, bam_md5) + [
            TTestExpectedOutput(
                tag=tag,
                suffix_secondary_file=".bai",
                preprocessor=TTestPreprocessor.FileSize,
                operator=operator.ge,
                expected_value=min_bai_size,
            ),
        ]
        if bai_md5 is not None:
            outcome += [
                TTestExpectedOutput(
                    tag=tag,
                    suffix_secondary_file=".bai",
                    preprocessor=TTestPreprocessor.FileMd5,
                    operator=operator.eq,
                    expected_value=bai_md5,
                ),
            ]
        return outcome



### BED ###


class Bed(File):
    def __init__(self, optional=False):
        super().__init__(optional=optional, extension=".bed")

    @staticmethod
    def name():
        return "bed"
    
class BedGz(Gunzipped):
    def __init__(self, optional=False):
        super().__init__(inner_type=Bed, optional=optional, extension=".bed.gz")

    @staticmethod
    def name():
        return "BedGz"

class BedTabix(FileTabix, BedGz):
    @staticmethod
    def name():
        return "BedTABIX"