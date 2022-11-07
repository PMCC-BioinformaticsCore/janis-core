from unittest import TestCase

from janis_core.workflow.workflow import WorkflowBuilder
from janis_core.types.common_data_types import Array, String

from janis_core.tests.testtools import BasicTestTool, ArrayStepTool, VersionTestTool


class TestContainers(TestCase):
    def test_command_tool(self):
        t = BasicTestTool()
        d = t.containers()
        self.assertEqual(t.versioned_id(), next(iter(d.keys())))
        self.assertEqual("ubuntu:latest", t.container())

    def test_workflow(self):
        w = WorkflowBuilder("test_container_workflow")
        w.input("inp", String)
        w.input("aInp", Array(String))

        w.step("stp1", BasicTestTool(testtool=w.inp))
        w.step("stp1_v2", VersionTestTool(testtool=w.inp))
        w.step("stp2", ArrayStepTool(inps=w.aInp))

        cons = w.containers()
        self.assertSetEqual(
            {"ArrayStepTool", "TestTranslationtool", "TestTranslationtool_v0_0_2"},
            set(cons.keys()),
        )

        self.assertIsNone(cons["ArrayStepTool"])
        self.assertEqual("ubuntu:latest", cons["TestTranslationtool"])
        self.assertEqual("ubuntu:latest", cons["TestTranslationtool_v0_0_2"])
