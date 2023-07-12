

from dataclasses import dataclass
from typing import Any, Optional

@dataclass
class ValidCheck:
    outname: str
    ctype: str
    value: Optional[Any] = None
    reffile: Optional[str] = None