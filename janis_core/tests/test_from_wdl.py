import unittest

from janis_core.ingestion.fromwdl import WdlParser


class TestFromWdl(unittest.TestCase):
    parser = WdlParser()

    def test_missing_disk(self):
        result = self.parser.parse_disk_requirement(None)
        self.assertEqual(2.14748, result)
