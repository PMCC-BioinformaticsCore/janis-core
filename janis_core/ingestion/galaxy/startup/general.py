
import os
import shutil
from typing import Any, Optional

from janis_core.ingestion.galaxy.logs import logging
from janis_core.ingestion.galaxy import settings
from janis_core.ingestion.galaxy import fileio
from janis_core.ingestion.galaxy import paths
from janis_core.ingestion.galaxy import datatypes


def general_setup(args: dict[str, Any]) -> None:
    settings.general.set_command(args['command']) 

    parent_outdir = format_parent_outdir(args)
    tool_outdir = format_tool_subdir(args)

    paths.set_outdir(parent_outdir)
    paths.set_tool_subdir(tool_outdir)

    setup_output_folder(args)
    setup_data_folder()

    logging.configure_logging()


def setup_data_folder() -> None:
    # get package data folder path
    package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    package_data_dir = f'{package_dir}/data'

    # create user data folder if not exists
    if not os.path.exists(paths.USER_DATA_DIR):
        os.mkdir(paths.USER_DATA_DIR)
    
    # create user wrapper download folder if not exists
    if not os.path.exists(f'{paths.USER_DATA_DIR}/{paths.DOWNLOADED_WRAPPERS_DIR}'):
        os.mkdir(f'{paths.USER_DATA_DIR}/{paths.DOWNLOADED_WRAPPERS_DIR}')

    # move files from package to user if they don't exist
    data_files = [
        paths.LOGGING_CONFIG,
        paths.GALAXY_CONFIG,
        paths.GALAXY_DATATYPES_YAML,
        paths.CONTAINER_CACHE,
        paths.WRAPPER_CACHE,
    ]
    for f in data_files:
        package_path = f'{package_data_dir}/{f}'
        user_path = f'{paths.USER_DATA_DIR}/{f}'
        if not os.path.exists(user_path):
            shutil.copyfile(package_path, user_path)

def format_tool_subdir(args: dict[str, Any]) -> Optional[str]:
    if args['command'] == 'tool':
        return None
    else:
        return 'tools'

def format_parent_outdir(args: dict[str, Any]) -> str:
    # user specified outdir
    if args['outdir']:
        return args['outdir']
    # auto formatted outdir
    else:
        parent_dir = auto_parent_dir(args)
        project_name = auto_project_name(args)
        return f'{parent_dir}/{project_name}'
    
def auto_parent_dir(args: dict[str, Any]) -> str:
    assert(args['outdir'] is None)
    match args['command']:
        case 'tool':
            return paths.DEFAULT_TOOL_OUTDIR
        case 'workflow':
            return paths.DEFAULT_WORKFLOW_OUTDIR
        case _:
            raise RuntimeError()

def auto_project_name(args: dict[str, Any]) -> str:
    match args['command']:
        case 'tool':
            return auto_tool_project_name(args)
        case 'workflow':
            return auto_workflow_project_name(args)
        case _:
            raise RuntimeError()

def auto_tool_project_name(args: dict[str, Any]) -> str:
    if args['infile']:
        return args['infile'].rsplit('/', 1)[-1].split('.', 1)[0]
    elif args['remote']:
        return args['remote'].split(',')[2]
    else:
        raise RuntimeError()

def auto_workflow_project_name(args: dict[str, Any]) -> str:
    return args['infile'].rsplit('/', 1)[-1].split('.', 1)[0]

def setup_output_folder(args: dict[str, Any]) -> None:
    # create outdir
    fileio.init_folder(paths.outdir())
    # create common subfolders
    for subfolder in paths.common_subfolders():
        fileio.init_folder(subfolder)
    # create workflow specific subfolders
    if args['command'] == 'workflow':
        for subfolder in paths.workflow_subfolders():
            fileio.init_folder(subfolder)


