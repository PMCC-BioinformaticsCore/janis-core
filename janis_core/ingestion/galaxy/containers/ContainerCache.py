

import json
import filelock
from typing import Any

from .Container import Container


class ContainerCache:
    def __init__(self, cache_path: str):
        self.cache_path = cache_path

    def get(self, tool: str) -> list[Container]:
        """returns all containers with this tool id"""
        cache = self._load_cache()
        out: list[Container] = []
        if tool in cache:
            for details_dict in cache[tool]:
                out.append(Container(details_dict))
        return out

    def add(self, container: Container, override: bool=False):
        tool = container.tool_name
        if self.exists(container):
            return
        else:
            cache = self._load_cache()
            if tool not in cache:
                cache[tool] = []
            cache[tool].append(container.__dict__)
            self._write_cache(cache)

    def exists(self, container: Container) -> bool:
        cache = self._load_cache()
        if container.tool_name in cache:
            for details in cache[container.tool_name]:
                if details['url'] == container.url:
                    return True
        return False

    def _load_cache(self) -> dict[str, list[dict[str, str]]]:
        try:
            lockpath = f"{self.cache_path.rsplit('.', 1)[0]}.lock"
            lock = filelock.FileLock(lockpath)
            with lock.acquire(timeout=10):
                with open(self.cache_path, 'r') as fp:
                    return json.load(fp)
        except FileNotFoundError:
            return {}

    def _write_cache(self, cache: dict[str, Any]) -> None:
        filepath = self.cache_path
        lockpath = f"{filepath.rsplit('.', 1)[0]}.lock"
        lock = filelock.FileLock(lockpath)
        with lock.acquire(timeout=10):
            with open(self.cache_path, 'w') as fp:
                json.dump(cache, fp)

