
 
from .directives import gen_directives_for_process
from .inputs import create_nextflow_process_inputs
from .outputs import create_nextflow_process_outputs
from .script import gen_script_for_cmdtool
from .process import gen_functions_for_process
from .process import gen_imports_for_process
from .process import gen_process_from_cmdtool
from .process import gen_process_from_codetool