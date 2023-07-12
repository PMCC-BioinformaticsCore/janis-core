
import os 
import unittest
from janis_core.ingestion.wdl import WdlParser

WDL_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/wdl')


class TestFromWdl(unittest.TestCase):
    parser = WdlParser()

    @unittest.skip("no wdl ingest yet")
    def test_ingest_tool(self) -> None:
        raise NotImplementedError

    @unittest.skip("no wdl ingest yet")
    def test_ingest_workflow(self) -> None:
        raise NotImplementedError