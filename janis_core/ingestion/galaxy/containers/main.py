



from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition
from janis_core.ingestion.galaxy.gx.gxtool.requirements.model import Requirement

from .mulled import create_mulled_container
from .fetch import fetch_container


DEFAULT_CONTAINER = 'python:3.7.16'

def resolve_dependencies_as_container(xmltool: XMLToolDefinition) -> str:
    # for each of the galaxy requirements, find a useable container from quay.io.
    # if more than 1 requirement, mull the containers together, build, and push to our dockerhub repo.
    # else, we have the container we need. 
    # return the container url. 

    if len(xmltool.metadata.requirements) == 0:
        return DEFAULT_CONTAINER
    else:
        return _resolve_single_requirement(xmltool)
    # elif len(xmltool.metadata.requirements) == 1:
    #     return _resolve_single_requirement(xmltool)
    # else:
    #     return _resolve_multiple_requirements(xmltool)

def _resolve_single_requirement(xmltool: XMLToolDefinition) -> str:
    main_req = xmltool.metadata.get_main_requirement()
    main_req_uri = fetch_container(main_req)
    if main_req_uri is None:
        return DEFAULT_CONTAINER
    return main_req_uri 

def _resolve_multiple_requirements(xmltool: XMLToolDefinition) -> str:
    base_container_uri = None
    other_reqs: list[Requirement] = []
    
    # start with main requirement as base container if we can find one
    main_req = xmltool.metadata.get_main_requirement()
    main_req_uri = fetch_container(main_req)

    if main_req_uri is not None:
        base_container_uri = main_req_uri
        other_reqs = [req for req in xmltool.metadata.requirements if req.name != main_req.name]

    else:
        container_uris = [fetch_container(req) for req in xmltool.metadata.requirements]
        container_uris = [uri for uri in container_uris if uri is not None]
        if len(container_uris) == 0:
            base_container_uri = DEFAULT_CONTAINER
        else:
            base_container_uri = container_uris[0]
    
    mulled_uri = create_mulled_container(xmltool, base_container_uri, other_reqs)
    return mulled_uri


