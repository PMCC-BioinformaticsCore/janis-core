

from typing import Any, Optional

from janis_core.ingestion.galaxy.gx.gxtool.requirements.model import Requirement
from janis_core.ingestion.galaxy.utils.general import global_align

from .Fetcher import Fetcher
from ..Container import Container
from . import utils


class GA4GHFetcher(Fetcher):
    """uses the GA4GH trs api to get tool containers and versions"""

    def fetch(self, requirement: Requirement) -> list[Container]:
        interactor = GA4GHInteractor()
        tools_list = interactor.request_tool_data(requirement.name)
        # TODO finish refactoring. GA4GH was down. 
        if tools_list:
            tool_data = self.select_tool(tools_list)
            version_data = self.get_version_data(tool_data)
            if version_data:
                return self.response_to_containers()
        return []

    def select_tool(self, tool_data: list[dict[str, Any]]) -> dict[str, Any]:
        raise NotImplementedError()

    def get_version_data(self, tool_data: dict[str, Any]) -> Optional[dict[str, Any]]:
        tool_data = self.most_similar_tool(api_results)
        version = self.most_similar_version(tool_data)
        interactor = GA4GHInteractor()
        return interactor.request_version_data(version['url'])

    def response_to_containers(self, data: dict[str, Any]) -> list[Container]:
        raise NotImplementedError()


    def most_similar_tool(containers: list[Container]) -> dict[str, Any]:
        """returns data for the api tool result most similar to the query tool name"""
        result_similarities: list[SimilarityScore] = []
        query_name = requirement.get_text()
        for tool_data in api_results:
            score = global_align(query_name, tool_data['name'])
            result_similarities.append(SimilarityScore(score, tool_data))
        result_similarities.sort(key = lambda x: x.score, reverse=True)
        return result_similarities[0].obj

class GA4GHInteractor:

    def request_tool_data(self, query_tool: str) -> Optional[list[dict[str, Any]]]:
        """
        make api request which returns dict of search results
        select the correct tool from results & return
        only makes requests which return lists of dicts
        """
        api_uri = self.format_tool_request_url(query_tool)
        return utils.make_api_request(api_uri)

    def request_version_data(self, query_url: str) -> list[dict[str, Any]]:
        return utils.make_api_request(query_url)

    def format_tool_request_url(self, tool_query: str) -> str:
        return f'https://api.biocontainers.pro/ga4gh/trs/v2/tools?name={tool_query}&limit=10&sort_field=id&sort_order=asc'

    


