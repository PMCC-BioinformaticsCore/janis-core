

from dataclasses import dataclass
from typing import Any
from datetime import datetime
from bioblend.toolshed import ToolShedInstance

from janis_core.ingestion.galaxy.gx.gxtool.requirements.model import CondaRequirement, ContainerRequirement, Requirement
from janis_core.ingestion.galaxy.gx.wrappers import Wrapper
from janis_core.ingestion.galaxy.runtime.dates import JANIS_DATE_FMT
from janis_core.ingestion.galaxy.runtime.dates import TOOLSHED_DATE_FMT


def request_installable_revision(gxstep: dict[str, Any]) -> list[Wrapper]:
    api_interactor = interactor_factory(gxstep)
    print(f"making galaxy API revision requests for {api_interactor.tool_id} v{api_interactor.tool_build}")
    return get_installable_revision(api_interactor)


@dataclass
class ToolshedAPIInteractor:
    tool_shed: str
    owner: str
    repo: str
    tool_id: str
    tool_build: str

    # TODO MULTIPLE REQUEST ATTEMPTS WITH TIMEOUT???

    def set_revision(self, revision: str) -> None:
        self.revision = revision

    def get_installable_revision(self) -> list[Wrapper]:
        ts = ToolShedInstance(f"https://{self.tool_shed}/")
        response = ts.repositories.get_repository_revision_install_info(self.repo, self.owner, self.revision) # type: ignore
        if response and len(response[0]) > 0:  # type: ignore
            return galaxy_response_to_wrappers(response) # type: ignore
        return []

    def get_ordered_revisions(self) -> list[str]:
        ts = ToolShedInstance(f"https://{self.tool_shed}/")
        return ts.repositories.get_ordered_installable_revisions(self.repo, self.owner) # type: ignore


def galaxy_response_to_wrappers(data: list[dict[str, Any]]) -> list[Wrapper]:
    wrappers: list[Wrapper] = []
    repository, metadata, _ = data
    for tool in metadata['valid_tools']:
        details = {
            'owner': repository['owner'],
            'repo': repository['name'],
            'revision': parse_revision(data),
            'tool_id': tool['id'],
            'tool_build': tool['version'],
            'date_created': parse_date_created(data),
            'requirements': parse_requirements(tool),
        }
        wrappers.append(Wrapper(details))
    return wrappers

def parse_requirements(tool: dict[str, Any]) -> list[Requirement]:
    out: list[Requirement] = []
    for req in tool['requirements']:
        if req['type'] == 'package':
            out.append(CondaRequirement(req['name'], req['version']))
        else:
            raise RuntimeError('got something other than conda requirement')
    return out

def parse_revision(data: list[dict[str, Any]]) -> str:
    install_info = data[2]
    assert(len(install_info)) == 1
    return list(install_info.values())[0][2]

def parse_date_created(data: list[dict[str, Any]]) -> str:
    """cast the galaxy toolshed date format to janis format"""
    repository = data[0]
    new_date = datetime.strptime(repository['create_time'], TOOLSHED_DATE_FMT)
    return new_date.strftime(JANIS_DATE_FMT)




def interactor_factory(gxstep: dict[str, Any]) -> ToolshedAPIInteractor:
    return ToolshedAPIInteractor(
        tool_shed = gxstep['tool_shed_repository']['tool_shed'],
        owner = gxstep['tool_shed_repository']['owner'],
        repo = gxstep['tool_shed_repository']['name'],
        tool_id = gxstep['tool_id'].rsplit('/', 2)[-2],
        tool_build = gxstep['tool_version']
    )

def get_installable_revision(api_interactor: ToolshedAPIInteractor) -> list[Wrapper]:
    revisions = api_interactor.get_ordered_revisions()
    for revision in reversed(revisions):
        wrappers = get_revision_wrappers(api_interactor, revision)
        if wrappers and version_match(wrappers, api_interactor):
            return wrappers
    raise RuntimeError

def get_revision_wrappers(api_interactor: ToolshedAPIInteractor, revision: str) -> list[Wrapper]:
    api_interactor.set_revision(revision)
    return api_interactor.get_installable_revision()



def version_match(wrappers: list[Wrapper], api_interactor: ToolshedAPIInteractor) -> bool:
    for wrapper in wrappers:
        if wrapper.tool_id == api_interactor.tool_id and wrapper.tool_build == api_interactor.tool_build:
            return True
    return False

# def get_installable_revision_old(gxstep: dict[str, Any]) -> list[Wrapper]:
#     api_interactor = interactor_factory(gxstep)
#     revision = gxstep['tool_shed_repository']['changeset_revision']
#     data = request_revision_exact(api_interactor, revision)
#     if not data:
#         data = request_revision_general(api_interactor)
#     return wrapper_factory(data)








