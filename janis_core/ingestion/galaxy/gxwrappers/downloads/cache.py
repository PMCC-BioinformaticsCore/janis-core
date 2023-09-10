

import os
import tarfile
from typing import Optional

from janis_core import settings
from janis_core.ingestion.galaxy.fileio import safe_init_folder


class DownloadCache: 
    """
    keeps track of the location of downloaded wrapper folders.
    DownloadCache.get() will return the local path to a tool xml if already downloaded
    DownloadCache.add() saves a tar as a download and notes its path. 
    """

    @property
    def path(self) -> str:
        if settings.testing.TESTMODE:
            return settings.ingest.galaxy.TESTING_WRAPPERS_DIR
        else:
            return settings.ingest.galaxy.DEFAULT_WRAPPERS_DIR
    
    def get(self, query_repo: str, query_revision: str) -> Optional[str]:
        """returns the local file path for the tool xml if already downloaded or None"""
        for folder in self._load():
            repo, revision = folder.split('-', 1)
            if repo == query_repo and revision == query_revision:
                return f'{self.path}{os.sep}{folder}'
        return None

    def add(self, tar: tarfile.TarFile) -> None:
        self._save(tar)

    def _save(self, tar: tarfile.TarFile) -> None:
        safe_init_folder(self.path)
        tar.extractall(path=self.path)

    def _load(self) -> set[str]:
        safe_init_folder(self.path)
        folders = os.listdir(self.path)
        folders = [f for f in folders if os.path.isdir(f'{self.path}{os.sep}{f}')]
        return set(folders)

