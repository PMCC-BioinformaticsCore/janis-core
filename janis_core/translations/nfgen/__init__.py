"""
    Nextflow modelling:
        Using DSL2

    Video: https://youtu.be/I-hunuzsh6A
    DSL2: https://www.nextflow.io/docs/latest/dsl2.html

"""

from .common import Channel, NFFile, Import
from .process import (
    Process,
    ProcessInput,
    ProcessOutput,
    ProcessScriptType,
    InputProcessQualifier,
    OutputProcessQualifier,
    TupleElementForOutput
)

from .workflow import Workflow, WorkflowInput, WorkflowOutput, WorkflowPublish

from .directives import *
