import unittest

from janis_core.redefinitions.tools import Cat
from janis_core.redefinitions.tools import GenerateVardictHeaderLines
from janis_core.redefinitions.workflows import BwaAligner
from janis_core.redefinitions.workflows import WGSGermlineMultiCallers


class TestJanisMisc(unittest.TestCase):
    def test_str_tool(self):
        BwaAligner().translate("janis")

    def test_str_python_tool(self):
        GenerateVardictHeaderLines().translate("janis")

    def test_command_tool(self):
        Cat().translate("janis")

    def test_str_big_workflow(self):
        WGSGermlineMultiCallers().translate("janis")
