
# from logging import getLogger, config
# from typing import Optional
# import yaml
# import sys 

import warnings

# from janis_core.ingestion.galaxy.settings import paths
from janis_core.ingestion.galaxy import settings

# -------------
# configuration
# -------------

warnings.filterwarnings("ignore")


# logging 
# def configure_logging() -> None:
#     config_path = f'{paths.USER_DATA_DIR}/{paths.LOGGING_CONFIG}'
#     with open(config_path, "r") as fp:
#         the_dict = yaml.safe_load(fp)
#         the_dict['handlers']['janis_log']['filename'] = paths.janis_log()
#         the_dict['handlers']['message_log']['filename'] = paths.message_log()
#     config.dictConfig(the_dict)


# -------
# logging
# -------

# messages
def msg_parsing_workflow(path: str):
    print(f'\nparsing workflow {path.split("/")[-1]}\n')

def msg_parsing_tool():
    filename = settings.tool.tool_path.rsplit('/', 1)[-1]
    print(f'parsing tool {filename}')

def msg_downloading_tool(url: str):
    print(f'downloading wrapper from {url}')


# # debug
# # (just runtime data and info messages which don't impact program)
# def runtime_data(data: str):
#     logger = getLogger('tool')
#     logger.debug(str(data))

# def evaluation_failed():
#     logger = getLogger('tool')
#     logger.debug('cheetah evaluation failed')


# # info
# # (things that MAY cause issues or just metric collecting)

# def has_preprocessing():
#     logger = getLogger('tool')
#     logger.info('preprocessing')

# def has_postprocessing():
#     logger = getLogger('tool')
#     logger.info('postprocessing')

# def has_backtick_statement():
#     logger = getLogger('tool')
#     logger.info('backtick statement')

# def has_multiline_str():
#     logger = getLogger('tool')
#     logger.info('multiline string')

# def has_repeat():
#     logger = getLogger('tool')
#     logger.info('repeat param')

# def has_cheetah_loop(text: Optional[str]=None):
#     logger = getLogger('tool')
#     logger.info(f'cheetah loop')

# def has_cheetah_function():
#     logger = getLogger('tool')
#     logger.info('cheetah function')

# def no_base_cmd():
#     logger = getLogger('tool')
#     logger.info('no base cmd')

# def has_configfile():
#     logger = getLogger('tool')
#     logger.info('configfile')

# def no_container():
#     logger = getLogger('tool')
#     logger.info('no container')

# def no_ga4gh_data():
#     logger = getLogger('tool')
#     logger.info('no ga4gh data')

# def container_version_mismatch():
#     logger = getLogger('tool')
#     logger.info('container version mismatch')

# def color_param_ignored():
#     logger = getLogger('tool')
#     logger.info('ignored unsupported param type: color')

# def uncertain_output():
#     logger = getLogger('tool')
#     logger.info('uncertain outputs')

# def unlinked_input_connection():
#     logger = getLogger('workflow')
#     logger.info('unknown input')


# # warning
# # (things that will probably cause issues and require human editing in the workflow)

# def no_inputs():
#     logger = getLogger('tool')
#     logger.warning('no inputs')

# def no_outputs():
#     logger = getLogger('tool')
#     logger.warning('no outputs')

# def zero_length_tag():
#     logger = getLogger('tool')
#     logger.warning('zero length tag')

# # error
# # (unsupported features)

# # critical
# # (program failed - uncaught exceptions)

# def tool_exception():
#     logger = getLogger('tool')
#     logger.critical('exception')
#     sys.exit(1)

# def no_close_quotation():
#     logger = getLogger('tool')
#     logger.critical('no closing quotation')
#     sys.exit(1)

# def workflow_exception():
#     logger = getLogger('workflow')
#     logger.critical('exception')
#     sys.exit(1)
