
from enum import Enum
from typing import Any

class SupportedTranslation(Enum):
    CWL = "cwl"
    WDL = "wdl"
    NEXTFLOW = "nextflow"
    JANIS = "janis"

    def __str__(self):
        return self.value

    def __eq__(self, other):
        return str(self) == str(other)
    
    @staticmethod
    def from_str(label):
        for st in SupportedTranslation:
            if st.value == label:
                return st
        raise RuntimeError(f"Could not find SupportedTranslation for {label}")

    def get_translator(self) -> Any:
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
        
        else:
            raise NotImplementedError

    @staticmethod
    def all():
        return [
            SupportedTranslation.CWL,
            SupportedTranslation.WDL,
            SupportedTranslation.JANIS,
            SupportedTranslation.NEXTFLOW
        ]
