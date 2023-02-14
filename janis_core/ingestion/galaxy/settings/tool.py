

import os
from typing import Any, Optional

from janis_core.ingestion.galaxy.gx.wrappers import Wrapper
from janis_core.ingestion.galaxy.gx.wrappers import fetch_wrapper
from janis_core.ingestion.galaxy.utils.galaxy import get_xml_id

tool_path: str
tool_id: str
owner: Optional[str] = None
repo: Optional[str] = None
revision: Optional[str] = None


# properties
def xml_basename() -> str:
    return tool_path.rsplit('/', 1)[-1].split('.')[0]

def xml_dir() -> str:
    if '/' not in tool_path:
        return '.'
    return tool_path.rsplit('/', 1)[0]

def logfile_path() -> str:
    return f'{xml_basename()}.log'

def set(from_args: Optional[dict[str, Any]]=None, from_wrapper: Optional[Wrapper]=None) -> None:
    if not from_args and not from_wrapper:
        raise RuntimeError('supply either args or wrapper to update')
    if from_args:
        update_args(from_args)
    elif from_wrapper:
        update_wrapper(from_wrapper)

def update_args(args: dict[str, Any]) -> None:
    global tool_path
    global tool_id
    global owner
    global repo
    global revision

    if args['infile']:
        rel_path = args['infile']
        tool_path = os.path.relpath(rel_path)
        tool_id = get_xml_id(tool_path)
        owner = None
        repo = None
        revision = None
    
    if args['remote']:
        owner, repo, tool_id, revision_raw = args['remote'].split(',')
        revision = revision_raw.rsplit(':', 1)[-1] # incase numeric:revision
        assert(owner)
        assert(repo)
        assert(revision)
        rel_path = fetch_wrapper(owner, repo, revision, tool_id)
        tool_path = os.path.relpath(rel_path)

def update_wrapper(wrapper: Wrapper) -> None:
    global tool_path
    global tool_id
    global owner
    global repo
    global revision

    tool_id = wrapper.tool_id
    owner = wrapper.owner
    repo = wrapper.repo
    revision = wrapper.revision
    rel_path = fetch_wrapper(wrapper.owner, wrapper.repo, wrapper.revision, wrapper.tool_id)
    tool_path = os.path.relpath(rel_path)

