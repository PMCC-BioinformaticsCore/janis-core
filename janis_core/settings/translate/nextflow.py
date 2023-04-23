




# params
# PARAM_VAR = "%PARAM%"
# LIST_OF_FILES_PARAM = "%LIST_OF_FILES_PARAM%"
# LIST_OF_FILE_PAIRS_PARAM = "%LIST_OF_FILE_PAIRS_PARAM%"
PYTHON_CODE_FILE_SYMBOL = "code_file"

# filenames
CONFIG_FILENAME = 'nextflow.config'
MAIN_WORKFLOW_NAME = 'main.nf'
PYTHON_CODE_OUTPUT_FILENAME_PREFIX = "out_"
PYTHON_SHEBANG = "#!/usr/bin/env python"
# (deprecated)
# OUTPUT_METADATA_FILENAME = "janis.outputs.metadata"
# LIB_FILENAME = "lib.nf"
# FINAL_STEP_NAME = "janis_outputs"
# TOOL_STDOUT_FILENAME = "janisstdout"
# NO_FILE_PATH_PREFIX = f"JANIS_NO_FILE"

# directories
BASE_OUTDIR = ''
PROCESS_OUTDIR = 'modules'
SUBWORKFLOW_OUTDIR = 'subworkflows'
TEMPLATES_OUTDIR = 'templates'

# text case 
NF_INDENT = '    '
NF_OUTDIR_CASE = 'snake'
NF_PARAM_CASE = 'snake'
NF_CHANNEL_CASE = 'snake'
NF_PROCESS_CASE = 'snake_caps'
NF_PROCESS_INPUT_CASE = 'snake'
NF_FILE_CASE = 'snake'
NF_MAIN_NAME = 'main'

# process translation
MINIMAL_PROCESS = True
# JANIS_ASSISTANT = False

# translation mode
MODE = 'workflow'


