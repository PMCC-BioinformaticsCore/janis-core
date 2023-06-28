

from dataclasses import dataclass


@dataclass
class Configfile:
    varname: str
    contents: str