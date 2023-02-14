




# params
PARAM_VAR = "%PARAM%"
LIST_OF_FILES_PARAM = "%LIST_OF_FILES_PARAM%"
LIST_OF_FILE_PAIRS_PARAM = "%LIST_OF_FILE_PAIRS_PARAM%"
PYTHON_CODE_FILE_SYMBOL = "code_file"

# filenames
OUTPUT_METADATA_FILENAME = "janis.outputs.metadata"
LIB_FILENAME = "lib.nf"
FINAL_STEP_NAME = "janis_outputs"
TOOL_STDOUT_FILENAME = "janisstdout"
CONFIG_FILENAME = "nextflow.config"
NO_FILE_PATH_PREFIX = f"JANIS_NO_FILE"
PYTHON_CODE_OUTPUT_FILENAME_PREFIX = "out_"
PYTHON_SHEBANG = "#!/usr/bin/env python"

# text case 
NF_INDENT = '    '
NF_OUTDIR_CASE = 'snake'
NF_PARAM_CASE = 'snake'
NF_CHANNEL_CASE = 'snake'
NF_PROCESS_CASE = 'snake_caps'
NF_PROCESS_INPUT_CASE = 'snake'
NF_MAIN_NAME = 'main'

# process translation
MINIMAL_PROCESS = True
JANIS_ASSISTANT = False

# comments
RENDER_COMMENTS = True

# translation mode
MODE = 'workflow'

# resources
WITH_RESOURCE_OVERRIDES = False

# directories
BASE_OUTDIR = ''
PROCESS_OUTDIR = 'modules'
SUBWORKFLOW_OUTDIR = 'subworkflows'
CODE_FILES_OUTDIR = 'templates'

# containers
WITH_CONTAINER = True
ALLOW_EMPTY_CONTAINER = False
CONTAINER_OVERRIDE = None