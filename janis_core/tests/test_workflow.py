from unittest import TestCase

from janis_core import File, Array, Logger, String, WorkflowBuilder
from janis_core.graph.steptaginput import StepTagInput, first_value, Edge
from janis_core.tests.testtools import SingleTestTool, ArrayTestTool


class TestWorkflow(TestCase):
    def setUp(self):
        Logger.mute()

    def tearDown(self):
        Logger.unmute()

    def test_name(self):
        wn = "test_name"
        w = WorkflowBuilder(wn)
        self.assertEqual(w.id(), wn)

    def test_rename(self):
        wn1 = "test_rename"
        wn2 = "test_rename2"
        w = WorkflowBuilder(wn1)
        w._identifier = wn2
        self.assertEqual(w.id(), wn2)

    def test_add_input(self):
        w = WorkflowBuilder("test_add_input")
        inp = w.input("inputLabel", str)
        self.assertEqual(len(w.input_nodes), 1)
        self.assertEqual(inp, next(iter(w.input_nodes.values())))
        self.assertIsNotNone(w.nodes[w.inputLabel.id()])

    def test_add_step(self):
        w = WorkflowBuilder("test_add_input")
        step = w.step("catStep", SingleTestTool(), ignore_missing=True)
        self.assertEqual(len(w.step_nodes), 1)
        self.assertEqual(step, next(iter(w.step_nodes.values())))
        self.assertIsNotNone(w.nodes[step.id()])

    def test_add_output(self):
        w = WorkflowBuilder("test_add_input")
        w.step("stp", SingleTestTool(), ignore_missing=True)
        w.output("outputStep", str, source=w.stp)
        self.assertEqual(len(w.output_nodes), 1)
        self.assertEqual(w.outputStep, next(iter(w.output_nodes.values())))
        self.assertIsNotNone(w.nodes["stp"])

    def WorkflowBuilder(self):
        w = WorkflowBuilder("test_add_node")
        inp = w.input("inp", str)
        stp = w.step("stp", SingleTestTool(), ignore_missing=True)

        self.assertEqual(len(w.input_nodes), 1)
        self.assertEqual(len(w.step_nodes), 1)
        self.assertEqual(len(w.output_nodes), 0)
        self.assertEqual(w.nodes["inp"].id(), inp.id())
        self.assertEqual(w.nodes["stp"].id(), stp.id())

    def test_add_qualified_edge(self):
        w = WorkflowBuilder("test_add_edge")
        inp = w.input("inp", str)
        stp = w.step("stp", SingleTestTool(inputs=w.inp))

        e = first_value(stp.sources["inputs"].source_map)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, first_value(stp.inputs()).id())

    def test_add_edge_later(self):
        w = WorkflowBuilder("test_add_edge")
        inp = w.input("inp", str)
        stp = w.step("stp", SingleTestTool(), ignore_missing=True)

        stp["inputs"] = inp

        e = first_value(stp.sources["inputs"].source_map)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, first_value(stp.inputs()).id())

    # def test_pipe(self):
    #     w = Workflow("test_add_edge")
    #     inp = Input("tarred", File())
    #     stp = Step("stp", SingleTestTool())  # Only has one input, with no output
    #     out = Output("outp", Array(File()))
    #
    #     w.add_pipe(inp, stp, out)
    #
    #     # the nodes are usually internal
    #     inp_node = w._nodes[inp.id()]
    #     stp_node = w._nodes[stp.id()]
    #     out_node = w._nodes[out.id()]
    #
    #     self.assertEqual(len(inp_node.connection_map), 0)
    #     self.assertEqual(len(stp_node.connection_map), 1)
    #     self.assertEqual(len(out_node.connection_map), 1)
    #
    #     s1: StepInput = first_value(stp_node.connection_map)
    #     s2: StepInput = first_value(out_node.connection_map)
    #
    #     e1 = first_value(s1.source_map)
    #     e2 = first_value(s2.source_map)
    #
    #     self.assertEqual(e1.start.id(), inp.id())
    #     self.assertEqual(e1.finish.id(), stp.id())
    #     self.assertEqual(e2.start.id(), stp.id())
    #     self.assertEqual(e2.finish.id(), out.id())
    #
    # def test_qualified_pipe(self):
    #     w = Workflow("test_add_edge")
    #     inp = Input("tarred", File())
    #     stp = Step("stp", SingleTestTool())  # Only has one input, with no output
    #     out = Output("outp", Array(File()))
    #
    #     w.add_pipe(inp, stp.inputs, out)
    #
    #     # the nodes are usually internal
    #     inp_node = w._nodes[inp.id()]
    #     stp_node = w._nodes[stp.id()]
    #     out_node = w._nodes[out.id()]
    #
    #     self.assertEqual(len(inp_node.connection_map), 0)
    #     self.assertEqual(len(stp_node.connection_map), 1)
    #     self.assertEqual(len(out_node.connection_map), 1)
    #
    #     s1: StepInput = first_value(stp_node.connection_map)
    #     s2: StepInput = first_value(out_node.connection_map)
    #
    #     e1: Edge = first_value(s1.source_map)
    #     e2: Edge = first_value(s2.source_map)
    #
    #     self.assertEqual(e1.start.id(), inp.id())
    #     self.assertEqual(e1.finish.id(), stp.id())
    #     self.assertEqual(e2.start.id(), stp.id())
    #     self.assertEqual(e2.finish.id(), out.id())

    def test_subworkflow(self):
        w = WorkflowBuilder("test_subworkflow")

        sub_w = WorkflowBuilder("subworkflow")
        sub_w.input("sub_inp", str)
        sub_w.step("sub_stp", SingleTestTool(inputs=sub_w.sub_inp))
        sub_w.output("sub_out", source=sub_w.sub_stp.out)

        w.input("inp", str)
        w.step("stp_workflow", sub_w(sub_inp=w.inp))
        w.output("out", source=w.stp_workflow.sub_out)

        # would be good to come up with some tests
        # w.translate("wdl")
        self.assertTrue(True)

    def test_add_scatter(self):
        w = WorkflowBuilder("scatterededge")
        w.input("inp", Array(str))
        stp = w.step("stp", SingleTestTool(inputs=w.inp), scatter="inputs")

        e = first_value(w.stp.sources["inputs"].source_map)

        self.assertTrue(e.compatible_types)
        self.assertListEqual(["inputs"], stp.scatter.fields)

    def test_add_non_scatter_fail(self):
        w = WorkflowBuilder("scatterededge")
        w.input("inp", Array(str))
        stp = w.step("stp", SingleTestTool(inputs=w.inp))

        e = first_value(w.stp.sources["inputs"].source_map)

        self.assertFalse(e.compatible_types)

    def test_add_scatter_incompatible(self):
        w = WorkflowBuilder("scatterededge")
        w.input("inp", Array(int))
        stp = w.step("stp", SingleTestTool(inputs=w.inp), scatter="inputs")

        e = first_value(w.stp.sources["inputs"].source_map)

        self.assertTrue(e.scatter)
        self.assertFalse(e.compatible_types)

    def test_add_scatter_nested_arrays(self):
        w = WorkflowBuilder("scatterededge")
        w.input("inp", Array(Array(str)))
        stp = w.step("stp", ArrayTestTool(inputs=w.inp), scatter="inputs")

        e = first_value(w.stp.sources["inputs"].source_map)

        self.assertTrue(e.compatible_types)
        self.assertListEqual(["inputs"], stp.scatter.fields)

    def test_add_scatter_nested_arrays_incompatible(self):
        w = WorkflowBuilder("scatterededge")
        w.input("inp", Array(Array(int)))
        stp = w.step("stp", ArrayTestTool(inputs=w.inp), scatter="inputs")

        e = first_value(w.stp.sources["inputs"].source_map)

        self.assertFalse(e.compatible_types)
        self.assertListEqual(["inputs"], stp.scatter.fields)

    def test_add_non_scatter(self):
        w = WorkflowBuilder("scatterededge")
        inp = w.input("inp", str)
        stp = w.step("stp", SingleTestTool(inputs=inp))

        e = first_value(w.stp.sources["inputs"].source_map)
        self.assertIsNone(stp.scatter)

    def test_add_non_scatter2(self):
        w = WorkflowBuilder("scatterededge")
        w.input("inp", Array(String()))
        w.step("stp", ArrayTestTool(inputs=w.inp))

        e = first_value(w.stp.sources["inputs"].source_map)
        self.assertFalse(e.scatter)

    def test_invalid_scatter_field(self):
        w = WorkflowBuilder("scatterededge")
        w.input("inp", Array(String()))
        self.assertRaises(
            Exception,
            w.step,
            identifier="stp",
            tool=ArrayTestTool(inputs=w.inp),
            scatter="randomfield",
        )

    def test_invalid_scatter_field_list(self):
        w = WorkflowBuilder("scatterededge")
        w.input("inp", Array(String()))
        self.assertRaises(
            Exception,
            w.step,
            identifier="stp",
            tool=ArrayTestTool(inputs=w.inp),
            scatter=["inputs", "randomfield"],
        )

    def test_merge(self):
        w = WorkflowBuilder("scatterededge")

        w.input("inp1", Array(String()))
        w.step("scatteredStp1", SingleTestTool(inputs=w.inp1), scatter="inputs")
        stp = w.step("mergeStp2", ArrayTestTool(inputs=w.scatteredStp1))

        e1 = first_value(w.scatteredStp1.sources["inputs"].source_map)
        e2 = first_value(w.mergeStp2.sources["inputs"].source_map)

        self.assertTrue(e1.scatter)
        self.assertFalse(e2.scatter)
        self.assertTrue(e2.compatible_types)

    def test_add_rescatter_scattered(self):
        w = WorkflowBuilder("scatterededge")

        w.input("inp1", Array(String()))
        stp1 = w.step("stp1", SingleTestTool(inputs=w.inp1), scatter="inputs")
        stp2 = w.step("stp2", SingleTestTool(inputs=stp1), scatter="inputs")

        e1 = first_value(stp1.sources["inputs"].source_map)
        e2 = first_value(stp2.sources["inputs"].source_map)

        self.assertTrue(e1.scatter)
        self.assertTrue(e2.scatter)

    def test_add_single_to_array_edge(self):
        w = WorkflowBuilder("test_add_single_to_array_edge")
        w.input("inp1", String())
        w.step("stp1", ArrayTestTool(inputs=w.inp1))

        e = first_value(w.stp1.sources["inputs"].source_map)
        self.assertTrue(w.has_multiple_inputs)
        self.assertTrue(e.compatible_types)
