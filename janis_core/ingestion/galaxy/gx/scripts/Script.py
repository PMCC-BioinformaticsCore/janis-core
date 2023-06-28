

from dataclasses import dataclass


@dataclass
class Script:
    varname: str
    filename: str
    contents: str