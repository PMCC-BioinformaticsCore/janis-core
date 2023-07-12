

import json
import os
from typing import Optional
from janis_core.ingestion.galaxy.gxwrappers import Wrapper
from janis_core.ingestion.galaxy.fileio import safe_init_file
from janis_core.settings.ingest.galaxy import WRAPPER_CACHE


"""
WrapperCache flat file structure:
{
    tool_id: [
            {
                'owner': 'fastqc',
                'repo': 'fastqc',
                'revision': 'fastqc',
                'tool_id': 'fastqc',
                'tool_version': 'fastqc',
                'tool_build': 'fastqc',
                'date_created': [revision date],
                'url': [toolshed download link],
            },
            ...
        ]
    }
}
"""



class WrapperCache:

    def get(
        self, 
        tool_id: Optional[str]=None,
        tool_build: Optional[str]=None,
        owner: Optional[str]=None,
        repo: Optional[str]=None,
        revision: Optional[str]=None
    ) -> list[Wrapper]:
        """loads cache then returns each Wrapper satisfying the query"""
        # load cache and get relevant wrapper in dict form
        cache = self._load()

        wrappers: list[dict[str, str]] = []
        for wrapper_list in cache.values():
            wrappers += wrapper_list
        
        # this can be more elegant ofc
        if tool_id:
            wrappers = [x for x in wrappers if x['tool_id'] == tool_id]
        if tool_build:
            wrappers = [x for x in wrappers if x['tool_build'] == tool_build]
        if owner:
            wrappers = [x for x in wrappers if x['owner'] == owner]
        if repo:
            wrappers = [x for x in wrappers if x['repo'] == repo]
        if revision:
            wrappers = [x for x in wrappers if x['revision'] == revision]
        
        # cast to Wrapper instances
        return [Wrapper(w) for w in wrappers]
    
    def add(self, wrapper: Wrapper) -> None:
        """adds the Wrapper to our cache and saves to file"""
        # load cache 
        cache = self._load() 

        # alter cache if needed
        self._ensure_key(wrapper.tool_id, cache)
        if not self.exists(wrapper, cache):
            cache[wrapper.tool_id].append(wrapper.to_dict())
        
        # write
        self._save(cache) 

    def exists(self, wrapper: Wrapper, cache: dict[str, list[dict[str, str]]]) -> bool:
        """checks if the wrapper is already in cache"""
        if wrapper.tool_id in cache:
            for entry in cache[wrapper.tool_id]:
                if entry['revision'] == wrapper.revision:
                    return True
        return False

    def _load(self) -> dict[str, list[dict[str, str]]]:
        """loads cache from flat file"""
        path = WRAPPER_CACHE
        if not os.path.exists(path):
            safe_init_file(path, contents='{}')
        with open(path, 'r') as fp:
            return json.load(fp)
    
    def _save(self, cache: dict[str, list[dict[str, str]]]) -> None:
        """saves our cache to file"""
        path = WRAPPER_CACHE
        if not os.path.exists(path):
            safe_init_file(path, contents='{}')
        with open(path, 'w') as fp:
            json.dump(cache, fp)

    def _ensure_key(self, key: str, cache: dict[str, list[dict[str, str]]]) -> None:    
        """ensures { key: [] } in cache"""
        if key not in cache:
            cache[key] = []
            
    
