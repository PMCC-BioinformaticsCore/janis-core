

from logging import getLogger, config
from .logfile import LogFile

import os
import warnings
import yaml
from pathlib import Path

# -------------
# configuration
# -------------

warnings.filterwarnings("ignore")

THIS_DIR = (os.path.dirname(os.path.realpath(__file__)))
CONFIG_FILE_PATH = f'{THIS_DIR}/config.yaml'
MESSAGE_LOG_PATH = f'./messages.log'


# logging 
def configure_logging() -> None:
    # delete previous log
    path = Path(MESSAGE_LOG_PATH)
    if path.exists():
        path.unlink()
    # set up logging package conf
    with open(CONFIG_FILE_PATH, "r") as fp:
        the_dict = yaml.safe_load(fp)
        the_dict['handlers']['messages']['filename'] = MESSAGE_LOG_PATH
    config.dictConfig(the_dict)


# console

def info_ingesting_workflow(spec: str, name: str) -> None:
    logger = getLogger('console')
    logger.info(f'ingesting {spec} workflow - "{name}"')

def info_ingesting_tool(spec: str, name: str) -> None:
    logger = getLogger('console')
    logger.info(f'ingesting {spec} tool - "{name}"')


logfile = LogFile(MESSAGE_LOG_PATH)

# message log
def log_message(uuid: str, msg: str) -> None:
    # add a message to the logfile (saves to file & adds to in-memory)
    global logfile
    logfile.add('INFO', uuid, msg)

def get_messages(uuid: str) -> list[str]:
    # get messages corresponding to this uuid from logfile
    global logfile
    loglines = logfile.messages[uuid]
    messages = [l.message for l in loglines]
    return messages

