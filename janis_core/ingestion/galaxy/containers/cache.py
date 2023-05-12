

from __future__ import annotations
import json
import filelock
import tempfile
from typing import Any, Optional

from .. import settings
from .Container import Container

DISABLE_CACHE = False


def init_cache() -> ContainerCache:
    if DISABLE_CACHE:
        temp = tempfile.TemporaryFile()
        cache_path = f'{tempfile.gettempdir()}/{temp.name}'
    else:
        cache_path = settings.paths.CONTAINER_CACHE
    return ContainerCache(cache_path)


class ContainerCache:
    def __init__(self, cache_path: str):
        self.cache_path = cache_path

    def get(self, versioned_tool_id: str) -> Optional[Container]:
        cache = self._load()
        if versioned_tool_id in cache:
            return Container(cache[versioned_tool_id])
        return None

    def add(self, versioned_tool_id:str, container: Container):
        cache = self._load()
        cache[versioned_tool_id] = container.__dict__
        self._write(cache)

    def _load(self) -> dict[str, dict[str, str]]:
        try:
            lockpath = f"{self.cache_path.rsplit('.', 1)[0]}.lock"
            lock = filelock.FileLock(lockpath)
            with lock.acquire(timeout=10):
                with open(self.cache_path, 'r') as fp:
                    return json.load(fp)
        except FileNotFoundError:
            return {}

    def _write(self, cache: dict[str, Any]) -> None:
        filepath = self.cache_path
        lockpath = f"{filepath.rsplit('.', 1)[0]}.lock"
        lock = filelock.FileLock(lockpath)
        with lock.acquire(timeout=10):
            with open(self.cache_path, 'w') as fp:
                json.dump(cache, fp)

