import requests
import os
import hashlib

from typing import Optional
from enum import Enum
from datetime import datetime


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
        elif (
            source.startswith(".")
            or source.startswith("/")
            or source.startswith("file://")
        ):
            self.type = FileSchemeType.local
        else:
            raise Exception(f"This file scheme is not currently supported {source}")

    def download(self, local_full_path: str) -> str:
        """
        Download remote file (currently only handling http/https)

        :param local_full_path: Full path to save the file to
        :type local_full_path: string
        :return: file path of the downloaded file
        :rtype: string
        """

        if self.type == FileSchemeType.http:
            resp = requests.get(
                self.source,
            )

            if resp.status_code != requests.codes.ok:
                raise Exception(
                    f"Failed to download remote file {self.source}: {resp.text}"
                )

            with open(local_full_path, "wb") as f:
                f.write(resp.content)

        elif self.type == FileSchemeType.local:
            pass

        if os.path.isfile(local_full_path):
            return local_full_path
        else:
            raise Exception(f"Failed to download file to {local_full_path}")

    def is_local(self) -> bool:
        """
        Check if a file is a local file

        :return: file is a local file or not
        :rtype: bool
        """
        return self.type == FileSchemeType.local

    def basename(self) -> str:
        """
        Find basename of a local/remote file

        :return: basename of a file
        :rtype: string
        """
        return os.path.basename(self.source)

    def hash_filename(self) -> str:
        """

        :return:
        :rtype:
        """
        # If we cannot find lat modified date, we will always download again
        components = self.source + (self.last_modified() or str(datetime.now()))
        hash_md5 = hashlib.md5(str.encode(components))
        hashed_filename = hash_md5.hexdigest()

        return hashed_filename

    def last_modified(self) -> Optional[str]:
        """
        Get the last modified date of this file

        :return: last modified datetime
        :rtype: string
        """
        if self.type == FileSchemeType.http:
            try:
                response = requests.head(self.source)
                last_modified = response.headers.get("Last-Modified")
                if last_modified:
                    return last_modified
                else:
                    return None
            except:
                return None

        return None
