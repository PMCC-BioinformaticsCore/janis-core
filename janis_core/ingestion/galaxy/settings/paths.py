
import os

package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

USER_DATA_DIR = f'{package_dir}/data'
GALAXY_CONFIG = f'{USER_DATA_DIR}/galaxy_config.yaml'
GALAXY_DATATYPES_YAML = f'{USER_DATA_DIR}/galaxy_janis_types.yaml'
CONTAINER_CACHE = f'{USER_DATA_DIR}/container_url_cache.json'
WRAPPER_CACHE = f'{USER_DATA_DIR}/wrappers.json'
DOWNLOADED_WRAPPERS_DIR = f'{USER_DATA_DIR}/downloads/tool_wrappers'


# LOGGING_CONFIG = f'{USER_DATA_DIR}/logging_config.yaml'
# DEFAULT_OUTDIR = 'parsed'
# DEFAULT_TOOL_OUTDIR = 'parsed/tools'
# DEFAULT_WORKFLOW_OUTDIR = 'parsed/workflows'

# COMMON_SUBFOLDERS: list[str] = [
#     'logs',
# ]
# WORKFLOW_SUBFOLDERS: list[str] = [
#     'subworkflows',
#     'tools',
#     'tools/wrappers',
#     'tools/scripts',
#     'tools/untranslated'
# ]


# _outdir: str = DEFAULT_OUTDIR
# _tool_subdir: Optional[str] = None

# def set_outdir(value: str) -> None:
#     global _outdir
#     _outdir = value

# def set_tool_subdir(value: Optional[str]) -> None:
#     global _tool_subdir
#     _tool_subdir = value


# # runtime files
# def outdir() -> str:
#     return _outdir

# def common_subfolders() -> list[str]:
#     return [f'{_outdir}/{sf}' for sf in COMMON_SUBFOLDERS]

# def workflow_subfolders() -> list[str]:
#     return [f'{_outdir}/{sf}' for sf in WORKFLOW_SUBFOLDERS]

# def janis_log() -> str:
#     return f'{_outdir}/logs/janis.log'

# def message_log() -> str:
#     return f'{_outdir}/logs/messages.log'


# # workflow files
# def tool(tool_id: str) -> str:
#     if _tool_subdir:
#         return f'{_outdir}/{_tool_subdir}/{tool_id}.py'
#     else:
#         return f'{_outdir}/{tool_id}.py'

# def workflow() -> str:
#     return f'{_outdir}/workflow.py'

# def inputs(file_format: str='yaml') -> str:
#     return f'{_outdir}/inputs.{file_format}'

# def config(file_format: str='yaml') -> str:
#     return f'{_outdir}/config.{file_format}'

# def subworkflow(tool_id: str) -> str:
#     return f'{_outdir}/subworkflows/{tool_id}.py'

# def untranslated(tool_id: str) -> str: 
#     return f'{_outdir}/tools/untranslated/{tool_id}.txt'

# def configfile(tool_id: str, configfile_name: str) -> str:
#     return f'{_outdir}/tools/scripts/{tool_id}_{configfile_name}'

# def script(filename: str) -> str:
#     return f'{_outdir}/tools/scripts/{filename}'

# def wrapper(tool_id: str, revision: str) -> str:
#     return f'{_outdir}/tools/wrappers/{tool_id}-{revision}'

