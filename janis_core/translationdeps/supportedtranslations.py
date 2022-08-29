from enum import Enum


class SupportedIngestion(Enum):
    JANIS = "janis"
    GALAXY = "galaxy"
    CWL = "cwl"
    #NEXTFLOW = "nf"    # future
    #WDL = "wdl"        # future

    def __str__(self) -> str:
        return self.value
    
    @staticmethod
    def all() -> list[str]:
        return [
            str(SupportedIngestion.JANIS),
            str(SupportedIngestion.GALAXY),
            str(SupportedIngestion.CWL),
        ]


class SupportedTranslation(Enum):
    CWL = "cwl"
    WDL = "wdl"
    NEXTFLOW = "nf"
    JANIS = "janis"

    def __str__(self):
        return self.value

    def __eq__(self, other):
        return str(self) == str(other)

    def get_translator(self):
        if self == SupportedTranslation.CWL:
            from ..translations.cwl import CwlTranslator

            return CwlTranslator()
        elif self == SupportedTranslation.WDL:
            from ..translations.wdl import WdlTranslator

            return WdlTranslator()
        elif self == SupportedTranslation.NEXTFLOW:
            from ..translations.nextflow import NextflowTranslator

            return NextflowTranslator()

        elif self == SupportedTranslation.JANIS:
            from ..translations.janis import JanisTranslator

            return JanisTranslator()

    @staticmethod
    def all():
        return [
            SupportedTranslation.CWL,
            SupportedTranslation.WDL,
            SupportedTranslation.JANIS,
            SupportedTranslation.NEXTFLOW
        ]
