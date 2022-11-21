"""
    Nextflow modelling:
        Using DSL2

    Video: https://youtu.be/I-hunuzsh6A
    DSL2: https://www.nextflow.io/docs/latest/dsl2.html

"""

from .common import NFFile, Import, ImportItem, Function
from .workflow import Workflow, WorkflowInput, WorkflowOutput, WorkflowPublish
from .directives import *  # so bad. 
from .process.process import (
    Process,
    ProcessInput,
    ProcessOutput,
    ProcessScriptType,
    InputProcessQualifier,
    OutputProcessQualifier,
    TupleElementForOutput
)

from . import settings
from . import channels
from .channels import ChannelOperation
from . import params
from . import ordering
from .config import generate_config
from .script import gen_script_for_cmdtool
from .unwrap import unwrap_expression, unwrap_source, translate_string_formatter
from .formatting import format_process_call
from .utils import to_groovy_str
from .register import register_workflow_inputs
from .process import gen_script_for_cmdtool





