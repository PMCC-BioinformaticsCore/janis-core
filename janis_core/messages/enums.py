
from enum import Enum

class ErrorLevel(Enum):
    INFO    = 'INFO'
    WARNING = 'WARNING'
    ERROR   = 'ERROR'

class ErrorCategory(Enum):
    SCRIPT          = 'SCRIPT'
    DATATYPE        = 'DATATYPE'
    PLUMBING        = 'PLUMBING'
    VERSION         = 'VERSION'
    FALLBACK        = 'FALLBACK'
    EXPERIMENTAL    = 'EXPERIMENTAL'
    FATAL           = 'FATAL'