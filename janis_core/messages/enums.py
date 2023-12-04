
from enum import Enum, auto


class FormatCategory(Enum):
    MAIN        = auto()
    INPUT       = auto()
    ARGUMENT    = auto()
    STEP        = auto()
    OUTPUT      = auto()

    def __str__(self) -> str:
        return self.name
    
    @staticmethod
    def from_str(category: str):
        if category == "MAIN":
            return FormatCategory.MAIN
        if category == "INPUT":
            return FormatCategory.INPUT
        if category == "STEP":
            return FormatCategory.STEP
        if category == "OUTPUT":
            return FormatCategory.OUTPUT
        raise ValueError(f"Unknown FormatCategory: {category}")


class ErrorCategory(Enum):
    FATAL           = ('ERROR','FATAL')
    FALLBACKS       = ('ERROR','FALLBACKS')
    DATATYPES       = ('WARNING','DATATYPES')
    PLUMBING        = ('WARNING','PLUMBING')
    METADATA        = ('WARNING','METADATA')
    EXPERIMENTAL    = ('WARNING','EXPERIMENTAL')
    SCRIPTING       = ('ERROR','SCRIPTING')

    @staticmethod
    def from_str(category: str):
        # this is awful
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
    