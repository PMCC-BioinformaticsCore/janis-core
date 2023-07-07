

from __future__ import annotations
import os
import json
import filelock
import tempfile
from typing import Any, Optional

from janis_core import settings
from janis_core.ingestion.galaxy.fileio import safe_init_file

from .Container import Container
from janis_core.ingestion.galaxy.fileio import safe_init_file


def init_cache() -> ContainerCache:
    if settings.ingest.galaxy.DISABLE_CONTAINER_CACHE:
        temp = tempfile.TemporaryFile()
        cache_path = os.path.join(tempfile.gettempdir(), str(temp.name))
        os.remove(cache_path)
        with open(cache_path, 'w') as fp:
            fp.write('{}')
    else:
        cache_path = settings.ingest.galaxy.CONTAINER_CACHE
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
        if not os.path.exists(self.cache_path):
            safe_init_file(self.cache_path, contents='{}')
        try:
            lockpath = f"{self.cache_path.rsplit('.', 1)[0]}.lock"
            lock = filelock.FileLock(lockpath)
            with lock.acquire(timeout=10):
                with open(self.cache_path, 'r') as fp:
                    return json.load(fp)
        except FileNotFoundError:
            return {}

    def _write(self, cache: dict[str, Any]) -> None:
        if not os.path.exists(self.cache_path):
            safe_init_file(self.cache_path, contents='{}')
        lockpath = f"{self.cache_path.rsplit('.', 1)[0]}.lock"
        lock = filelock.FileLock(lockpath)
        with lock.acquire(timeout=10):
            with open(self.cache_path, 'w') as fp:
                json.dump(cache, fp)

