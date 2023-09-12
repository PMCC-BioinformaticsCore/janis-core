

from galaxy.tool_util.deps.mulled.util import v2_image_name
from galaxy.tool_util.deps.mulled.util import build_target

from janis_core import settings
from janis_core.ingestion.galaxy.gxtool.model import XMLTool


DEFAULT_CONTAINER = 'quay.io/biocontainers/python:3.10.1'


def resolve_dependencies_as_container(xmltool: XMLTool) -> str:
    # for each of the galaxy requirements, find a useable container from quay.io.
    # return the container url. 
    if settings.testing.TESTING_USE_DEFAULT_CONTAINER:
        return DEFAULT_CONTAINER

    elif len(xmltool.metadata.requirements) == 0:
        return DEFAULT_CONTAINER
    
    elif len(xmltool.metadata.requirements) == 1:
        req = xmltool.metadata.requirements[0]
        return f'quay.io/biocontainers/{req.name}:{req.version}'
    
    else:
        items = [build_target(req.name, version=req.version) for req in xmltool.metadata.requirements]
        resource = v2_image_name(items)
        return f'quay.io/biocontainers/{resource}'