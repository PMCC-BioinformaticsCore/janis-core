
from datetime import datetime
from typing import Any


from janis_core.ingestion.galaxy.runtime.dates import JANIS_DATE_FMT
from janis_core.ingestion.galaxy.gx.gxtool.requirements import Requirement


class Wrapper:
    def __init__(self, details: dict[str, Any]):
        self.owner: str = details['owner']
        self.repo: str = details['repo']
        self.revision: str = details['revision']
        self.tool_id: str = details['tool_id']
        self.tool_build: str = details['tool_build']
        self.date_created: datetime = datetime.strptime(details['date_created'], JANIS_DATE_FMT)
        self.requirements: list[Requirement] = details['requirements']
        self.inbuilt: bool = False
        
    @property
    def tool_version(self) -> str:
        return self.tool_build.split('+galaxy')[0]
    
    @property
    def url(self) -> str:
        return f'https://toolshed.g2.bx.psu.edu/repos/{self.owner}/{self.repo}/archive/{self.revision}.tar.gz'

    def to_dict(self) -> dict[str, Any]:
        return {
            'owner': self.owner,
            'repo': self.repo,
            'revision': self.revision,
            'tool_id': self.tool_id,
            'tool_version': self.tool_version,
            'tool_build': self.tool_build,
            'date_created': self.date_created.strftime(JANIS_DATE_FMT),
            'requirements': [req.__dict__ for req in self.requirements],
            'url': self.url,
        }
