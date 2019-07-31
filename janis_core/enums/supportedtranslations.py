from enum import Enum

SupportedTranslation = str


class SupportedTranslations(Enum):
    CWL = "cwl"
    WDL = "wdl"

    def __str__(self):
        return self.value

    def __eq__(self, other):
        return str(self) == str(other)

    @staticmethod
    def all():
        return [SupportedTranslations.CWL, SupportedTranslations.WDL]
