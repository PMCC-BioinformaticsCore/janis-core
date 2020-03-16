import unittest
from janis_core.tool.documentation import InputQualityType
from janis_core.tests.testtools import TestInputQualityTool


class TestCommandToolInputGeneration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        cls.inmap = {
            "user": InputQualityType.user,
            "static": InputQualityType.static,
            "configuration": InputQualityType.configuration,
            "none": None,
        }

    def test_regular(self):
        expected_keys = set(self.inmap.keys())
        inputs = TestInputQualityTool().generate_inputs_override()
        self.assertSetEqual(expected_keys, set(inputs.keys()))

    def test_get_user_inputs(self):
        qualtype = InputQualityType.user
        expected_keys = set(i for i in self.inmap.keys() if self.inmap[i] == qualtype)
        inputs = TestInputQualityTool().generate_inputs_override(
            quality_type=[qualtype]
        )
        self.assertEqual(len(inputs), 1)
        self.assertSetEqual(expected_keys, set(inputs.keys()))

    def test_get_static_inputs(self):
        qualtype = InputQualityType.static
        expected_keys = set(i for i in self.inmap.keys() if self.inmap[i] == qualtype)
        inputs = TestInputQualityTool().generate_inputs_override(
            quality_type=[qualtype]
        )
        self.assertEqual(len(inputs), 1)
        self.assertSetEqual(expected_keys, set(inputs.keys()))

    def test_get_static_and_config_inputs(self):
        qualtypes = [InputQualityType.static, InputQualityType.configuration]
        expected_keys = set(i for i in self.inmap.keys() if self.inmap[i] in qualtypes)
        inputs = TestInputQualityTool().generate_inputs_override(quality_type=qualtypes)
        self.assertEqual(len(inputs), 2)
        self.assertSetEqual(expected_keys, set(inputs.keys()))

    def test_ignore(self):
        tool = TestInputQualityTool()
        ignore_keys = {"none"}

        expected_keys = set(i.id() for i in tool.inputs() if i.id() not in ignore_keys)
        inputs = tool.generate_inputs_override(values_to_ignore=ignore_keys)
        self.assertSetEqual(expected_keys, set(inputs.keys()))
