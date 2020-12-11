import os
from swiftclient import Connection
from swiftclient.service import SwiftService
from swiftclient import client


class ObjectStoreFile:
    def __init__(self, container: str, object_key: str, md5: str):
        self.container = container
        self.object_key = object_key
        self.md5 = md5

    def download(self, conn: Connection, dir: str):
        conn.get_object(container=self.container, obj=self.object_key)

    @staticmethod
    def connect():
        _authurl = "https://keystone.rc.nectar.org.au:5000/v3/"
        _auth_version = "3"
        _user = os.environ.get("OS_USERNAME")
        _key = os.environ.get("OS_PASSWORD")
        _os_options = {
            "user_domain_name": os.environ.get("OS_USER_DOMAIN_NAME"),
            "project_domain_name": os.environ.get("OS_PROJECT_DOMAIN_ID"),
            "project_name": os.environ.get("OS_PROJECT_NAME"),
        }

        conn = Connection(
            authurl=_authurl,
            user=_user,
            key=_key,
            os_options=_os_options,
            auth_version=_auth_version,
        )

        return conn


if __name__ == "__main__":
    conn = ObjectStoreFile.connect()
    o = ObjectStoreFile("janis-test-data", "bioinformatics")

    o.download(conn)
