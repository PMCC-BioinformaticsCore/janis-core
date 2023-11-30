
from enum import Enum


class ErrorCategory(Enum):
    DATATYPES       = ('WARNING','DATATYPES')
    PLUMBING        = ('WARNING','PLUMBING')
    METADATA        = ('WARNING','METADATA')
    EXPERIMENTAL    = ('WARNING','EXPERIMENTAL')
    FALLBACKS       = ('ERROR','FALLBACKS')
    SCRIPTING       = ('ERROR','SCRIPTING')
    FATAL           = ('ERROR','FATAL')

    @staticmethod
    def from_str(category: str):
        if category == "DATATYPES":
            return ErrorCategory.DATATYPES
        if category == "PLUMBING":
            return ErrorCategory.PLUMBING
        if category == "METADATA":
            return ErrorCategory.METADATA
        if category == "EXPERIMENTAL":
            return ErrorCategory.EXPERIMENTAL
        if category == "FALLBACKS":
            return ErrorCategory.FALLBACKS
        if category == "SCRIPTING":
            return ErrorCategory.SCRIPTING
        if category == "FATAL":
            return ErrorCategory.FATAL
        raise ValueError(f"Unknown ErrorCategory: {category}")
    