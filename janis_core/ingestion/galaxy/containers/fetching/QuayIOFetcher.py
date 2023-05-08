


from typing import Any
from janis_core.ingestion.galaxy.gx.gxtool.requirements.model import Requirement

from ..Container import Container
from ..fetching.Fetcher import Fetcher
from . import utils


class QuayIOFetcher(Fetcher):

    def fetch(self, requirement: Requirement) -> list[Container]:
        interactor = QuayInteractor()
        tool_data = interactor.request_tool_data(requirement.name)
        if tool_data:
            return self.response_to_containers(tool_data, requirement)
        else:
            return []

    def response_to_containers(self, repo: dict[str, Any], requirement: Requirement) -> list[Container]:
        out: list[Container] = []
        for details in repo['tags'].values():
            info: dict[str, str] = {
                'image_type': 'docker',
                'repo': repo['name'],
                'tag': details['name'],
                'url': self.format_url(repo, details),
                '_timestamp': details['last_modified']
            }
            out.append(Container(info))
        return out
    
    def format_url(self, repo: dict[str, Any], tag: dict[str, str]) -> str:
        return f"quay.io/biocontainers/{repo['name']}:{tag['name']}"


class QuayInteractor:

    def request_tool_data(self, name: str) -> dict[str, Any]:
        endpoint = f'https://quay.io/api/v1/repository/biocontainers/{name}'
        return utils.make_api_request(endpoint)


    