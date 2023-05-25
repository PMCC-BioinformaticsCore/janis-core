

from typing import Any
from janis_core.ingestion.galaxy.model.workflow import StepMetadata
from janis_core.ingestion.galaxy.gx.wrappers import Wrapper
from janis_core.ingestion.galaxy.gx.wrappers import WrapperCache
from janis_core.ingestion.galaxy.gx.wrappers.requests.versions import request_single_wrapper


"""
This module parses a galaxy step (.ga format) and produces StepMetadata.
During the process it identifies which toolshed Wrapper will be used.

The Wrapper is an important part of StepMetadata. 
It contains information about the tool_id, tool_build, revision,
date_created, conda requirements etc

Local information about Wrappers is stored in data/wrappers/wrappers.json.
If information about which Wrapper to use is not available locally, 
API requests are made to the galaxy toolshed to determine the most suitable wrapper. 
its information is then known and saved to the json file. 

The WrapperCache class is responsible for presenting our local knowledgebank of
wrappers during this phase.
"""


def parse_step_metadata(gxstep: dict[str, Any]) -> StepMetadata:
    return StepMetadata(
        wrapper=get_wrapper(gxstep),
        step_id=int(gxstep['id']),
        step_name=gxstep['name'],
        tool_state=gxstep['tool_state'], 
        workflow_outputs=gxstep['workflow_outputs'],
        _label=gxstep['label']
    )

def get_wrapper(gxstep: dict[str, Any]) -> Wrapper:
    """
    gets a galaxy Wrapper.
    wrapper isn't downloaded - just details are fetched in like owner, repo etc. 
    wrapper download and parsing happens later. 
    wrapper may be pulled from cache of known wrappers, or may a new one may need
    to be created using toolshed API calls. 
    """
    if 'tool_shed_repository' not in gxstep:
        return get_wrapper_builtin(gxstep)
    else:
        return get_wrapper_toolshed(gxstep)

def get_wrapper_builtin(gxstep: dict[str, Any]) -> Wrapper:
    details: dict[str, Any] = {
        'owner': 'Galaxy',
        'repo': 'None',
        'revision': 'None',
        'tool_id': gxstep['tool_id'],
        'tool_build': gxstep['tool_version'],
        'date_created': '1960-01-01 00:00:00',
        'requirements': [],
    }
    wrapper = Wrapper(details)
    wrapper.inbuilt = True
    return wrapper

def get_wrapper_toolshed(gxstep: dict[str, Any]) -> Wrapper:
    wrappers = get_local(gxstep)
    if not wrappers:
        wrappers = get_toolshed(gxstep)
    return most_recent(wrappers)
        
def get_local(gxstep: dict[str, Any]) -> list[Wrapper]:
    cache = WrapperCache()
    return cache.get(
        owner=gxstep['tool_shed_repository']['owner'],
        repo=gxstep['tool_shed_repository']['name'],
        tool_id=gxstep['tool_id'].rsplit('/', 2)[-2],
        tool_build=gxstep['tool_version']
    )

def get_toolshed(gxstep: dict[str, Any]) -> list[Wrapper]:
    # scrape toolshed for wrappers and update cache
    wrapper = request_single_wrapper(
        tool_shed=gxstep['tool_shed_repository']['tool_shed'],
        owner=gxstep['tool_shed_repository']['owner'],
        repo=gxstep['tool_shed_repository']['name'],
        tool_id=gxstep['tool_id'].rsplit('/', 2)[-2],
        tool_build=gxstep['tool_version']
    )
    cache = WrapperCache()
    cache.add(wrapper)
    
    # confirm can now load wrapper details from cache
    local_wrappers = get_local(gxstep)
    if not local_wrappers:
        raise RuntimeError()
    
    return local_wrappers

def most_recent(wrappers: list[Wrapper]) -> Wrapper:
    """return the most recent wrapper which matches [version]"""
    wrappers.sort(key=lambda x: x.date_created, reverse=True)
    return wrappers[0]

