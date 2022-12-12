

# from typing import Optional
# from typing import Any


# params
PARAM_VAR = "%PARAM%"
LIST_OF_FILES_PARAM = "%LIST_OF_FILES_PARAM%"
LIST_OF_FILE_PAIRS_PARAM = "%LIST_OF_FILE_PAIRS_PARAM%"
PYTHON_CODE_FILE_PATH_PARAM = "%PYTHON_CODE_FILE_PATH%"

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
NF_MAIN_NAME = 'main'

# process translation
MINIMAL_PROCESS = True
JANIS_ASSISTANT = False

# translation mode
MODE = 'workflow'


# def configure(key: Optional[str]=None, val: Optional[Any]=None, dict_config: Optional[dict[str, Any]]=None) -> None:
#     if key is not None and val is not None:
#         configure_keyval(key, val)
#     elif dict_config is not None:
#         configure_dict(dict_config)
#     else:
#         raise RuntimeError

# def configure_keyval(key: str, val: Any) -> None:
#     globals()[key] = val

# def configure_dict(dict_config: dict[str, Any]) -> None:
#     for key, val in dict_config.items():
#         globals()[key] = val
