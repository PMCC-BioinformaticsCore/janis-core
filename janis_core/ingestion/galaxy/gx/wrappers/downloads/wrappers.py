
from janis_core.ingestion.galaxy.logs import logging
import requests
import tarfile
import os
from typing import Optional

from janis_core.ingestion.galaxy.utils import galaxy as utils
from janis_core.ingestion.galaxy.gx.wrappers.downloads.cache import DownloadCache


CACHE: DownloadCache = DownloadCache()


def fetch_wrapper(owner: str, repo: str, revision: str, tool_id: str) -> str:
    """gets the wrapper locally or from toolshed then returns path to xml"""
    path: Optional[str] = None
    if not path:
        path = _fetch_cache(repo, revision, tool_id)
    if not path:
        path = _fetch_builtin(tool_id)
    if not path:
        path = _fetch_toolshed(owner, repo, revision, tool_id)
    if not path:
        raise RuntimeError(f'could not find wrapper for {tool_id}:{revision}')
    else:
        return path

def get_builtin_tool_path(tool_id: str) -> Optional[str]:
    """returns path to xml file with id='tool_id'"""
    tool_directories = _get_builtin_tool_directories()
    for directory in tool_directories:
        xmlfile = utils.get_xml_by_id(directory, tool_id)
        if xmlfile:
            return f'{directory}/{xmlfile}'
    return None

def _get_builtin_tool_directories() -> list[str]:
    out: list[str] = []
    out += _get_builtin_tools_directories()
    out += _get_datatype_converter_directories()
    return out

def _get_builtin_tools_directories() -> list[str]:
    import galaxy.tools
    tools_folder = str(galaxy.tools.__file__).rsplit('/', 1)[0]
    bundled_folders = os.listdir(f'{tools_folder}/bundled')
    bundled_folders = [f for f in bundled_folders if not f.startswith('__')]
    bundled_folders = [f'{tools_folder}/bundled/{f}' for f in bundled_folders]
    bundled_folders = [f for f in bundled_folders if os.path.isdir(f)]
    return [tools_folder] + bundled_folders

def _get_datatype_converter_directories() -> list[str]:
    import galaxy.datatypes
    datatypes_folder = str(galaxy.datatypes.__file__).rsplit('/', 1)[0]
    converters_folder = f'{datatypes_folder}/converters'
    return [converters_folder]

def _fetch_builtin(tool_id: str) -> Optional[str]:
    return get_builtin_tool_path(tool_id)

def _fetch_cache(repo: str, revision: str, tool_id: str) -> Optional[str]:
    wrapper = CACHE.get(repo, revision)
    if wrapper:
        xml = utils.get_xml_by_id(wrapper, tool_id)
        if xml:
            return f'{wrapper}/{xml}'
    return None

def _fetch_toolshed(owner: str, repo: str, revision: str, tool_id: str) -> Optional[str]:
    # download and add to cache
    url = _get_url_via_revision(owner, repo, revision)
    # logging.msg_downloading_tool(url)
    tar = _download_wrapper(url)
    CACHE.add(tar)
    # fetch from cache
    return _fetch_cache(repo, revision, tool_id)

def _get_url_via_revision(owner: str, repo: str, revision: str) -> str:
    return f'https://toolshed.g2.bx.psu.edu/repos/{owner}/{repo}/archive/{revision}.tar.gz'

def _download_wrapper(url: str) -> tarfile.TarFile:
    response = requests.get(url, stream=True)
    return tarfile.open(fileobj=response.raw, mode='r:gz')


