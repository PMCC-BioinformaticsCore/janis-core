"""
    Nextflow modelling:
        Using DSL2

    Video: https://youtu.be/I-hunuzsh6A
    DSL2: https://www.nextflow.io/docs/latest/dsl2.html

"""
CONFIG_FILENAME = "nextflow.config"

from .common import NFFile, Import, ImportItem, Function
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
from .channels import channel_factory, ChannelDeclaration, ChannelDeclarationBlock
from .params import ParamDeclaration, ParamDeclarationBlock
from .directives import *
