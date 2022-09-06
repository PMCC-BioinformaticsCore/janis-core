
from enum import Enum


class SupportedIngestion(Enum):
    JANIS = "janis"
    GALAXY = "galaxy"
    CWL = "cwl"
    WDL = "wdl"        
    #NEXTFLOW = "nf"    # future

    @staticmethod
    def all() -> list[str]:
        """return the value of each enum item"""
        return [str(item.value) for item in SupportedIngestion]