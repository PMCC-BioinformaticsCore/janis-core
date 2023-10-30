

from logging import getLogger, config
from .logfile import LogFile
from .logfile import LogLine


from typing import Optional
import os
import warnings
import yaml
from pathlib import Path
from enum import Enum
"""
This file was built for ingest / translation message logging. 

If there are any errors / unsupported features encountered during 
ingestion, you can use the functions here to log them as messages to disk.

During translation, you can then retrieve those error messages and display them 
to the user. 

the uuid field allows you to keep track of which messages belong to which 
janis entity in an ingested workflow / tool. 

for example, if we created an Edge during ingestion where we couldn't identify the 
data source, the we can provide uuid=Edge.uuid while logging this message. 
In the translate unit, you could then use get_messages(uuid=Edge.uuid()) to pull that 
message. 

this is the same as giving janis entities a 'error_messages: list[str] = []' field,
but decided against this.

you can also provide uuid=None, which would make sense for generalised issues
not related to a specific entity. 

currently, these janis entities have uuids:
    janis_core.graph.steptaginput.Edge  
    (issues around connecting data source to data sink)

    janis_core.graph.steptaginput.StepTagInput  
    (issues around step inputs more generally - scatter, step input doesn't 
    map to a tool input etc)


"""

# -------------
# configuration
# -------------

class ErrorCategory(Enum):
    SCRIPT      = 'script'
    DATATYPE    = 'datatype'
    DATAFLOW    = 'dataflow'
    VERSION     = 'version'
    FALLBACK    = 'fallback'
    FATAL       = 'fatal'

warnings.filterwarnings("ignore")

WORK_DIR = (os.getcwd())  # the working directory where janis translate is run from
PACKAGE_DIR = (os.path.dirname(os.path.realpath(__file__)))  #the directory where this file is located
CONFIG_FILE_PATH = f'{PACKAGE_DIR}/config.yaml'
MESSAGE_LOG_PATH = f'{WORK_DIR}/.janis/messages.log'

def configure_logging() -> None:
    # ensure we create the folder
    if not os.path.exists(os.path.dirname(MESSAGE_LOG_PATH)):
        os.mkdir(os.path.dirname(MESSAGE_LOG_PATH))

    # delete previous log
    path = Path(MESSAGE_LOG_PATH)
    if path.exists():
        path.unlink()

    # set up logging package conf
    with open(CONFIG_FILE_PATH, "r") as fp:
        the_dict = yaml.safe_load(fp)
        the_dict['handlers']['messages']['filename'] = MESSAGE_LOG_PATH
    config.dictConfig(the_dict)


# ----------
# to console
# ----------
def info_ingesting_workflow(spec: str, name: str) -> None:
    logger = getLogger('console')
    logger.info(f'ingesting {spec} workflow - "{name}"')

def info_ingesting_tool(spec: str, name: str) -> None:
    logger = getLogger('console')
    logger.info(f'ingesting {spec} tool - "{name}"')

# -------
# to file
# -------

def _log_message(level: str, category: Optional[ErrorCategory], uuid: Optional[str], msg: str) -> None:
    logfile = LogFile(MESSAGE_LOG_PATH)
    # if no uuid provided, consider this a general message provided during ingestion / translation. 
    # these messages can be shown to the user at the top of the main parsed file (ie the main workflow / tool), 
    # or you could generate a file in the output folder for the user to show this info. 
    if uuid is None:
        uuid = 'general'
    if category is not None:
        cat = category.value        
    logfile.add(level, cat, uuid, msg)

def log_info(uuid: Optional[str], msg: str) -> None:
    """logs a message which is considered INFO to the logfile"""
    _log_message('INFO', None, uuid, msg)

def log_warning(uuid: str, msg: str) -> None:
    """logs a message which is considered a WARNING to the logfile"""
    _log_message('WARNING', None, uuid, msg)

def log_error(uuid: Optional[str], category: ErrorCategory, msg: str) -> None:
    """logs a message which is considered an ERROR to the logfile"""
    _log_message('ERROR', category, uuid, msg)

def get_messages(uuid: Optional[str]=None, level: Optional[str]=None) -> list[LogLine]:
    logfile = LogFile(MESSAGE_LOG_PATH)
    loglines = logfile.lines
    
    # filter for specific entity
    if uuid is not None:
        loglines = [x for x in loglines if x.uuid == uuid]
    
    # filter for specific type of log message if 'level' supplied
    if level is not None:
        valid_levels = ['INFO', 'WARNING', 'ERROR']
        if level not in valid_levels:
            raise RuntimeError
        loglines = [x for x in loglines if x.level == level]

    return loglines

