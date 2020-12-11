from enum import Enum
import requests
import os
import hashlib


class FileSourceType(Enum):
    http = "http"
    local = "local"


class FileSource:
    def __init__(self, file_source: str):
        """
        :param file_source: full path to the file source
        :type file_source: string

        """
        self.file_source = file_source

        if file_source.startswith("http://") or file_source.startswith("https://"):
            self.type = FileSourceType.http
        else:
            self.type = FileSourceType.local

    def download(self, local_full_path: str):
        print(f"downloading file {self.file_source}")
        resp = requests.get(
            self.file_source,
        )

        with open(local_full_path, "wb") as f:
            f.write(resp.content)

    def is_local(self):
        return self.type == FileSourceType.local

    def basename(self):
        return os.path.basename(self.file_source)

    def hash_filename(self):
        hash_md5 = hashlib.md5(str.encode(self.file_source))
        hashed_filename = hash_md5.hexdigest()

        return hashed_filename
