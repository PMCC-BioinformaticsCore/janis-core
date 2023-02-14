



from typing import Optional
from dataclasses import dataclass


@dataclass
class JanisDatatype:
    format: str
    source: str
    classname: str
    extensions: Optional[str]
    import_path: str