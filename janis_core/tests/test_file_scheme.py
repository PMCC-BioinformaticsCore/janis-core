import os
from typing import Optional, Dict
from unittest import TestCase, mock
from janis_core.utils.file_scheme import FileSchemeType, FileScheme

valid_url = "https://abc.com/some_dir/myfile.txt"
valid_url_last_modified = "Thu, 03 Dec 2020 02:07:45 GMT"
valid_url_content = "test data content"


class MockResponse:
    def __init__(
        self,
        json_data: Dict,
        status_code: int,
        headers: Optional[Dict] = None,
        content: Optional[str] = None,
    ):
        self.json_data = json_data
        self.status_code = status_code
        self.headers = headers
        self.content = content.encode()

    def json(self):
        return self.json_data


# This method will be used by the mock to replace requests.get
def mocked_requests_head_and_get(*args, **kwargs):

    if args[0] == valid_url:
        return MockResponse(
            {}, 200, {"Last-Modified": valid_url_last_modified}, valid_url_content
        )

    return MockResponse(None, 404)


class TestFileScheme(TestCase):
    def setUp(self):
        self.test_dir = os.path.join(os.getcwd(), "test_remote_files")
        os.makedirs(self.test_dir, exist_ok=True)

    def test_file_type(self):
        url = "http://abc.com/some_dir/myfile.txt"
        fs = FileScheme(url)
        assert fs.type == FileSchemeType.http
        assert fs.is_local() == False

        url = "https://abc.com/some_dir/myfile.txt"
        fs = FileScheme(url)
        assert fs.type == FileSchemeType.http
        assert fs.is_local() == False

        url = "/some_dir/myfile.txt"
        fs = FileScheme(url)
        assert fs.type == FileSchemeType.local
        assert fs.is_local() == True

        url = "./some_dir/myfile.txt"
        fs = FileScheme(url)
        assert fs.type == FileSchemeType.local
        assert fs.is_local() == True

        url = "file://some_dir/myfile.txt"
        fs = FileScheme(url)
        assert fs.type == FileSchemeType.local
        assert fs.is_local() == True

    @mock.patch("requests.head", side_effect=mocked_requests_head_and_get)
    def test_last_modified(self, mock_head):
        url = valid_url
        fs = FileScheme(url)

        assert fs.last_modified() == valid_url_last_modified

    @mock.patch("requests.head", side_effect=mocked_requests_head_and_get)
    def test_hash_filename(self, mock_head):
        url = valid_url
        fs = FileScheme(url)

        assert fs.hash_filename() == "3a30c4d465eff8d4ff2123a9d2d0169c"

    @mock.patch("requests.get", side_effect=mocked_requests_head_and_get)
    def test_download(self, mock_get):
        url = valid_url
        fs = FileScheme(url)

        local_full_path = os.path.join(self.test_dir, "downloaded.txt")
        fs.download(local_full_path)
        assert os.path.exists(local_full_path)

        with open(local_full_path, "rb") as f:
            content = f.read()

        assert content == valid_url_content.encode()

    def test_basename(self):
        url = valid_url
        fs = FileScheme(url)

        assert fs.basename() == "myfile.txt"
