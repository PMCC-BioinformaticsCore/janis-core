
from janis_core import settings
from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition

from .mulled import gen_mulled_image
from .mulled import get_base_image_uri
from .fetch import fetch_online
from .fetch import DEFAULT_CONTAINER
from .cache import init_cache

from .Container import Container


def resolve_dependencies_as_container(xmltool: XMLToolDefinition) -> str:
    # for each of the galaxy requirements, find a useable container from quay.io.
    # return the container url. 
    if len(xmltool.metadata.requirements) == 0:
        return DEFAULT_CONTAINER.uri
    
    versioned_id = xmltool.metadata.versioned_id
    cache = init_cache()
    container = cache.get(versioned_id)
    
    # already had a container for this tool version (early exit)
    if container is not None:
        return container.uri
    
    # get container online (1 requirement)
    if len(xmltool.metadata.requirements) == 1:
        container = fetch_online(xmltool.metadata.requirements[0])
    
    # gen mulled image for all requirements (2+ requirements)
    elif settings.ingest.galaxy.GEN_IMAGES:
        container = gen_mulled_image(xmltool)
    
    # get container online for 1 requirement, prioritising what we think is the main requirement (2+ requirements)
    else:
        uri, _ = get_base_image_uri(xmltool)
        container = Container({
            'image_type': 'docker',
            'repo': uri.split('/')[-1].split(':')[0],
            'tag': uri.split('/')[-1].split(':')[1],
            'uri': uri,
            '_timestamp': '',
        })

    if container:
        cache.add(versioned_id, container)
    
    container = container if container else DEFAULT_CONTAINER
    return container.uri


    



