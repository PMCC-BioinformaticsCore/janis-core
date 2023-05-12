

from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition

from .mulled import gen_mulled_image
from .fetch import fetch_online
from .fetch import DEFAULT_CONTAINER
from .cache import init_cache


def resolve_dependencies_as_container(xmltool: XMLToolDefinition) -> str:
    # for each of the galaxy requirements, find a useable container from quay.io.
    # if more than 1 requirement, mull the containers together, build, and push to our dockerhub repo.
    # else, we have the container we need. 
    # return the container url. 
    if len(xmltool.metadata.requirements) == 0:
        return DEFAULT_CONTAINER.uri
    
    versioned_id = xmltool.metadata.versioned_id
    cache = init_cache()
    container = cache.get(versioned_id)
    
    # already had a container for this tool version
    if container is not None:
        return container.uri
    
    # get container online or gen mulled image
    if len(xmltool.metadata.requirements) == 1:
        container = fetch_online(xmltool.metadata.requirements[0])
    else:
        container = gen_mulled_image(xmltool)

    if container:
        cache.add(versioned_id, container)
    
    container = container if container else DEFAULT_CONTAINER
    return container.uri


    



