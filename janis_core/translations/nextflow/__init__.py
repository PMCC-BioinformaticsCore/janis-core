"""
    Nextflow modelling:
        Using DSL2

    Video: https://youtu.be/I-hunuzsh6A
    DSL2: https://www.nextflow.io/docs/latest/dsl2.html

"""

from .model import process
from .model import files
from .model import workflow

from . import trace
from . import channels
from . import params
from . import ordering
from . import workflow
from . import process
from . import naming
from . import data_sources
from . import parsing
from . import preprocessing

from .plumbing import gen_task_call
from .plumbing import is_datatype_mismatch
from .plumbing import gen_datatype_mismatch_plumbing
from .plumbing import call

from .main import NextflowTranslator

from .scope import Scope
from .channels import ChannelOperation
from .channels import Channel
from .parsing.config import generate_config
from .unwrap import unwrap_expression
from .nfgen_utils import to_groovy





