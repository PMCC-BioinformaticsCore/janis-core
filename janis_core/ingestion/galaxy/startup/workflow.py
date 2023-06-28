

from janis_core.ingestion.galaxy.logs import logging

from typing import Any
from janis_core.ingestion.galaxy.runtime.exceptions import InputError
from janis_core.ingestion.galaxy.utils import galaxy as utils


def workflow_setup(args: dict[str, Any]) -> None:
    path = args['infile']
    logging.msg_parsing_workflow(path)
    _validate_workflow_settings(path)

### VALIDATION ###

def _validate_workflow_settings(path: str) -> None:
    if not _valid_workflow(path):
        raise InputError('please check workflow file path')

def _valid_workflow(path: str) -> bool:
    if utils.is_galaxy_workflow(path):
        return True
    return False

