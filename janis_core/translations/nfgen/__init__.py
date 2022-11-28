"""
    Nextflow modelling:
        Using DSL2

    Video: https://youtu.be/I-hunuzsh6A
    DSL2: https://www.nextflow.io/docs/latest/dsl2.html

"""

from .common import NFFile, Import, ImportItem, Function
from .workflow import Workflow, WorkflowTake, WorkflowEmit
from .directives import *  # so bad. 

from .process import (
    Process,
    ProcessInput,
    ProcessOutput,
    ProcessScriptType,
    # InputProcessQualifier,
    # OutputProcessQualifier,
    # TupleElementForOutput
)

from . import settings
from . import channels
from . import params
from . import ordering
from . import process

from .call import get_args
from .channels import ChannelOperation
from .config import generate_config
from .unwrap import unwrap_expression, translate_string_formatter
from .formatting import format_process_call
from .nfgen_utils import to_groovy
from .register import register_params_channels






