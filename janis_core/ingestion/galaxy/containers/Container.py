



from datetime import datetime

from janis_core.ingestion.galaxy.runtime.dates import JANIS_DATE_FMT
from janis_core.ingestion.galaxy.runtime.dates import QUAY_DATE_FMT



class Container: 
    def __init__(self, info: dict[str, str]):
        self.image_type: str = info['image_type']
        self.repo: str = info['repo']
        self.tag: str = info['tag']
        self.uri: str = info['uri']
        self._timestamp: str = info['_timestamp']

    @property
    def registry_host(self) -> str:
        return self.uri.split('/', 1)[0]

    @property
    def tool_name(self) -> str:
        return self.repo

    @property
    def tool_version(self) -> str:
        if self.registry_host == 'quay.io':
            return self.tag.split('--', 1)[0]
        raise NotImplementedError()

    @property
    def timestamp(self) -> str:
        """map the timestamp from different formats to unified output"""
        if self.registry_host == 'quay.io': 
            date = datetime.strptime(self._timestamp, QUAY_DATE_FMT)
        else:
            raise NotImplementedError()
        return date.strftime(JANIS_DATE_FMT)


