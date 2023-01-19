"""
    Nextflow modelling:
        Using DSL2

    Video: https://youtu.be/I-hunuzsh6A
    DSL2: https://www.nextflow.io/docs/latest/dsl2.html

"""

from .common import NFFile, Import, ImportItem, Function
from .workflow import Workflow, WorkflowTake, WorkflowEmit
from .directives import (
    DiskDirective,
    CpusDirective,
    TimeDirective,
    DebugDirective,
    CacheDirective,
    MemoryDirective,
    ProcessDirective,
    PublishDirDirective,
    ContainerDirective
)

from .process import (
    Process,
    ProcessScriptType,
)
from .process.outputs import ProcessOutput
from .process.inputs import (
    ProcessInput,
    ValProcessInput,
    PathProcessInput,
    TupleProcessInput
)

from . import settings
from . import channels
from . import params
from . import ordering
from . import process
from . import naming
from .plumbing import call

from .scope import Scope
from .channels import ChannelOperation
from .channels import Channel
from .config import generate_config
from .unwrap import unwrap_expression
from .plumbing import format_process_call
from .nfgen_utils import to_groovy
from .register import register_params_channels






