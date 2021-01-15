from enum import Enum


class SupportedTranslation(Enum):
    CWL = "cwl"
    WDL = "wdl"
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

        elif self == SupportedTranslation.JANIS:
            from ..translations.janis import JanisTranslator

            return JanisTranslator()

    @staticmethod
    def all():
        return [
            SupportedTranslation.CWL,
            SupportedTranslation.WDL,
            SupportedTranslation.JANIS,
        ]
