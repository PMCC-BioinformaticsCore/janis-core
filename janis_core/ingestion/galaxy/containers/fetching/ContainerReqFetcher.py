

from janis_core.ingestion.galaxy.gxtool.model import XMLContainerRequirement

from ..Container import Container
from .Fetcher import Fetcher


class ContainerReqFetcher(Fetcher):
    
    def fetch(self, requirement: XMLContainerRequirement) -> list[Container]:
        info: dict[str, str] = {
            'image_type': requirement.image_type,
            'repo': requirement.name,
            'tag': requirement.version,
            'uri': requirement.uri,
            '_timestamp': 'Tue, 14 June 1994 23:56:01 -0000', # my bday lol
        }
        return [Container(info)]
