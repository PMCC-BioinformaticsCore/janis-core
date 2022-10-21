"""
    Nextflow modelling:
        Using DSL2

    Video: https://youtu.be/I-hunuzsh6A
    DSL2: https://www.nextflow.io/docs/latest/dsl2.html

"""
CONFIG_FILENAME = "nextflow.config"

from .common import NFFile, Import, ImportItem, Function
from .workflow import Workflow, WorkflowInput, WorkflowOutput, WorkflowPublish
from .channels import channel_factory, ChannelDeclaration, ChannelDeclarationBlock
from .params import ParamDeclaration, ParamDeclarationBlock
from .directives import *  # so bad. 
from .process import (
    Process,
    ProcessInput,
    ProcessOutput,
    ProcessScriptType,
    InputProcessQualifier,
    OutputProcessQualifier,
    TupleElementForOutput
)

from .script import gen_script_for_cmdtool
from .prescript import gen_prescript_for_cmdtool
from .unwrap import unwrap_expression, translate_string_formatter
from .formatting import format_process_call
from .utils import wrap_value
