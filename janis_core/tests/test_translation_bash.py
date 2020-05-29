import unittest

from janis_core.tests.testtools import SingleTestTool
from janis_core.translations.bash import BashTranslator


class TestCwlTypesConversion(unittest.TestCase):
    bash = BashTranslator().translate_tool_internal(SingleTestTool())
    print(bash)
