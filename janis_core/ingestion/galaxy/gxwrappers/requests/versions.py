

from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Optional
from datetime import datetime
from bioblend.toolshed import ToolShedInstance

from janis_core.ingestion.galaxy.gxtool.model import XMLRequirement
from janis_core.ingestion.galaxy.gxtool.model import XMLCondaRequirement
from janis_core.ingestion.galaxy.gxwrappers import Wrapper
from janis_core.ingestion.galaxy.runtime.dates import JANIS_DATE_FMT
from janis_core.ingestion.galaxy.runtime.dates import TOOLSHED_DATE_FMT


def request_single_wrapper(tool_shed: str, owner: str, repo: str, tool_id: str, tool_build: str) -> Wrapper:
    api_interactor = ToolshedAPIInteractor(tool_shed, owner, repo, tool_id, tool_build)
    print(f"making galaxy API request for {api_interactor.tool_id} v{api_interactor.tool_build}")
    wrapper = api_interactor.get_single_wrapper()
    return wrapper


@dataclass
class ToolshedAPIInteractor:
    tool_shed: str
    owner: str
    repo: str
    tool_id: str
    tool_build: str

    # TODO MULTIPLE REQUEST ATTEMPTS WITH TIMEOUT???
    def get_single_wrapper(self) -> Wrapper:
        ts = ToolShedInstance(f"https://{self.tool_shed}/")
        revisions = ts.repositories.get_ordered_installable_revisions(self.repo, self.owner) # type: ignore
        if not revisions:
            raise RuntimeError
        for revision in reversed(revisions):
            revision_data = ts.repositories.get_repository_revision_install_info(self.repo, self.owner, revision) # type: ignore
            if revision_data and len(revision_data[0]) > 0:  # type: ignore
                repo_wrappers = self._galaxy_revision_info_to_wrappers(revision_data) # type: ignore
                tool_wrapper = self._select_version_match_strict(repo_wrappers)
                if tool_wrapper:
                    return tool_wrapper
        raise RuntimeError

    def _select_version_match_strict(self, wrappers: list[Wrapper]) -> Optional[Wrapper]:
        for wrapper in wrappers:
            if wrapper.tool_id == self.tool_id and wrapper.tool_build == self.tool_build:
                return wrapper
        return None
    
    def _select_version_match_lenient(self, wrappers: list[Wrapper]) -> Wrapper:
        for wrapper in wrappers:
            if wrapper.tool_id == self.tool_id and wrapper.tool_build == self.tool_build:
                return wrapper
        return wrappers[0]

    def _galaxy_revision_info_to_wrappers(self, revision_data: list[dict[str, Any]]) -> list[Wrapper]:
        wrappers: list[Wrapper] = []
        repository, metadata, install_info = revision_data
        for tool_info in metadata['valid_tools']:
            details = {
                'owner': repository['owner'],
                'repo': repository['name'],
                'revision': self._parse_revision(install_info),
                'tool_id': tool_info['id'],
                'tool_build': tool_info['version'],
                'date_created': self._parse_date_created(revision_data),
                'requirements': self._parse_requirements(tool_info),
            }
            wrappers.append(Wrapper(details))
        return wrappers

    def _parse_requirements(self, tool_info: dict[str, Any]) -> list[XMLRequirement]:
        return [XMLCondaRequirement(req['name'], req['version']) for req in tool_info['requirements'] if req['type'] == 'package']

    def _parse_revision(self, install_info: dict[str, Any]) -> str:
        return list(install_info.values())[0][2]

    def _parse_date_created(self, data: list[dict[str, Any]]) -> str:
        """cast the galaxy toolshed date format to janis format"""
        repository = data[0]
        new_date = datetime.strptime(repository['create_time'], TOOLSHED_DATE_FMT)
        return new_date.strftime(JANIS_DATE_FMT)



