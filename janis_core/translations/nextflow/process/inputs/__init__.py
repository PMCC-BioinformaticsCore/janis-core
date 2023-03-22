

from .factory import create_nextflow_process_inputs

from .data_sources import update
from .data_sources import process_inputs
from .data_sources import param_inputs
from .data_sources import internal_inputs

from .model import ProcessInput
from .model import PathProcessInput
from .model import ValProcessInput
from .model import TupleProcessInput