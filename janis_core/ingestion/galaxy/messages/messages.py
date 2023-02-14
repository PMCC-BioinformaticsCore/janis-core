


from dataclasses import dataclass
from enum import Enum, auto
from typing import Tuple



class MsgCategory(Enum):
    WARNING     = auto()
    FIX         = auto()
    CHECK       = auto()
    ATTENTION   = auto()


@dataclass
class Message:
    category: MsgCategory
    title: str
    contents: str
    file_reference: str
    line_reference: Tuple[int, int]
