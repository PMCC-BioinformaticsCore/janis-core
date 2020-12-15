from enum import Enum
import requests
import os
import hashlib


class FileSchemeType(Enum):
    http = "http"
    local = "local"


class FileScheme:
    def __init__(self, source: str):
        """
        :param source: full path to the file source
        :type source: string

        """
        self.source = source

        if source.startswith("http://") or source.startswith("https://"):
            self.type = FileSchemeType.http
        else:
            self.type = FileSchemeType.local

    def download(self, local_full_path: str):
        resp = requests.get(
            self.source,
        )

        with open(local_full_path, "wb") as f:
            f.write(resp.content)

    def is_local(self):
        return self.type == FileSchemeType.local

    def basename(self):
        return os.path.basename(self.source)

    def hash_filename(self):
        hash_md5 = hashlib.md5(str.encode(self.source))
        hashed_filename = hash_md5.hexdigest()

        return hashed_filename
