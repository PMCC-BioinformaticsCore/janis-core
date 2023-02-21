

from typing import Any
from .groups import TagGroup

groups: dict[str, TagGroup] = {}


def register(entity: Any) -> None:  
    """register a tag in the active TagGroup"""
    group = _get_active()
    group.register(entity.name, entity) 

def get(uuid: str) -> str:
    """get a tag from any TagGroup"""
    global groups
    for group in groups.values():
        tag = group.get(uuid)
        if tag:
            return tag
    raise RuntimeError('No tag registered')

def new_group(section: str, uuid: str):
    """
    create new TagGroup (ie for a new workflow or tool), 
    then switch to that TagGroup to start registering tags
    """
    global groups
    groups[uuid] = TagGroup(section)
    switch_group(uuid)

def switch_group(uuid: str):
    """swap to a specific TagGroup to register new tags in that group"""
    _clear_active()
    _set_active(uuid)
            

def _get_active() -> TagGroup:
    global groups
    for group in groups.values():
        if group.active:
            return group
    raise RuntimeError('no active group')

def _set_active(query_uuid: str):
    global groups
    for uuid, group in groups.items():
        if uuid == query_uuid:
            group.active = True
            return
    raise RuntimeError('group doesnt exist')

def _clear_active():
    global groups
    for group in groups.values():
        group.active = False






