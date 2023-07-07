

from dataclasses import dataclass


@dataclass
class XMLScript:
    varname: str
    filename: str
    contents: str