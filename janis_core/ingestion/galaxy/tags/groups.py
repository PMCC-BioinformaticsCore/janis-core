


from typing import Any, Optional
from .strategies import format_tag


class TagGroup:

    def __init__(self, section: str):
        self.section = section
        self.uuids_basetags: dict[str, str] = {}
        self.basetags_uuids: dict[str, list[str]] = {}
        self.active: bool = False

    def exists(self, uuid: str) -> bool:
        if uuid in self.uuids_basetags:
            return True
        return False
    
    def get(self, uuid: str) -> Optional[str]:
        if uuid in self.uuids_basetags:
            basetag = self.uuids_basetags[uuid]
            return self._generate_tag(basetag, uuid)
    
    def get_base_tag(self, uuid: str) -> str:
        return self.uuids_basetags[uuid]

    def register(self, starting_text: str, entity: Any) -> None:
        basetag = format_tag(starting_text, entity)
        uuid = entity.uuid

        self.uuids_basetags[uuid] = basetag
        if basetag not in self.basetags_uuids:
            self.basetags_uuids[basetag] = []
        self.basetags_uuids[basetag].append(uuid)

    def _generate_tag(self, basetag: str, query_uuid: str) -> Optional[str]:
        stored_uuids = self.basetags_uuids[basetag]
        if len(stored_uuids) <= 1:
            return basetag # only 1 object using this basetag
        for i, uuid in enumerate(stored_uuids):
            if uuid == query_uuid:
                return f'{basetag}{i+1}' # appends '1', '2' etc if basetag is shared by multiple objects
        return None
