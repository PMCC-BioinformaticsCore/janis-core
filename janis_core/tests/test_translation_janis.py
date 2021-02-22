import unittest

run_tests = False


@unittest.skipUnless(run_tests, reason="These tests are fully fleshed out yet")
class TestJanisMisc(unittest.TestCase):
    def test_str_tool(self):
        from janis_bioinformatics.tools.common import BwaAligner

        BwaAligner().translate("janis")

    def test_str_python_tool(self):
        from janis_bioinformatics.tools.pmac.generatevardictheaderlines import (
            GenerateVardictHeaderLines,
        )

        GenerateVardictHeaderLines().translate("janis")

    def test_command_tool(self):
        from janis_unix.tools import Cat

        Cat().translate("janis")

    def test_str_big_workflow(self):
        from janis_pipelines.wgs_germline.wgsgermline import WGSGermlineMultiCallers

        WGSGermlineMultiCallers().translate("janis")
