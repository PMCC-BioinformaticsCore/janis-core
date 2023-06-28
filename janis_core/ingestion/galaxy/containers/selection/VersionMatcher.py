


from dataclasses import dataclass
from typing import Any, Optional

from janis_core.ingestion.galaxy import expressions


@dataclass
class Version:
    text: str

    def get_numeric(self) -> str:
        matches = expressions.get_matches(self.text, expressions.patterns.VERSION)
        if matches:
            match = matches[0]
            return match[0] # type: ignore
        raise RuntimeError()

    def get_numeric_levels(self, fill_size: int=0) -> list[int]:
        numeric_version = self.get_numeric()
        split_version = numeric_version.split('.')
        split_version = [int(x) for x in split_version]
        while len(split_version) < fill_size:
            split_version.append(0)
        return split_version
        

class VersionMatcher:

    def get_version_exact(self, tool_data: dict[str, Any], target_version: str) -> Optional[dict[str, str]]:
        for ver in tool_data['versions']:
            if ver['meta_version'] == target_version:
                return ver
        return None

    def get_version_trimmed(self, tool_data: dict[str, Any], target_version: str) -> Optional[dict[str, str]]:
        target = Version(target_version)

        for ver in tool_data['versions']:
            query = Version(ver['meta_version'])
            if query.get_numeric():
                if query.get_numeric() == target.get_numeric():
                    return ver
        return None

    def get_version_trimmed_inexact(self, tool_data: dict[str, Any], target_version: str) -> Optional[dict[str, str]]:
        """
        selects the version by trimming to numeric, then matching most similar, allowing wobble 
        just does a version comparision sort for each digit, right to left. 
        sorts are stable, so ordering comes out correct
        """
        target = Version(target_version)
        queries = [Version(ver['meta_version']) for ver in tool_data['versions']]
        array_size = max([len(ver.get_numeric_levels()) for ver in [target] + queries])

        for i in range(array_size -1, -1, -1):
            queries.sort(key=lambda x: abs(target.get_numeric_levels()[i] - x.get_numeric_levels(fill_size=array_size)[i]))
        
        selected = queries[0].text
        #return the right version dict
        for query in tool_data['versions']:
            if query['meta_version'] == selected:
                return query

    def get_most_recent(self, tool_data: dict[str, Any]) -> dict[str, str]:
        return tool_data['versions'][-1]
        