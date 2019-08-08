from unittest import TestCase

from janis_core import (
    File,
    Array,
    Logger,
    CommandTool,
    ToolInput,
    String,
    Input,
    Output,
    Step,
    Workflow,
    ToolOutput,
)
from janis_core.graph.stepinput import StepInput, first_value, Edge


class SingleTestTool(CommandTool):
    @staticmethod
    def tool():
        return "TestStepTool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [ToolInput("inputs", String())]

    def friendly_name(self):
        return None

    def outputs(self):
        return [ToolOutput("out", String())]

    @staticmethod
    def container():
        return None

    @staticmethod
    def version():
        return None


class ArrayTestTool(CommandTool):
    @staticmethod
    def tool():
        return "ArrayStepTool"

    def friendly_name(self):
        return None

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [ToolInput("inputs", Array(String()))]

    def outputs(self):
        return [ToolOutput("outs", Array(String()))]

    @staticmethod
    def container():
        return None

    @staticmethod
    def version():
        return None


class TestWorkflow(TestCase):
    def setUp(self):
        Logger.mute()

    def tearDown(self):
        Logger.unmute()

    def test_name(self):
        wn = "test_name"
        w = Workflow(wn)
        self.assertEqual(w.identifier, wn)

    def test_rename(self):
        wn1 = "test_rename"
        wn2 = "test_rename2"
        w = Workflow(wn1)
        w.identifier = wn2
        self.assertEqual(w.identifier, wn2)

    def test_add_input(self):
        w = Workflow("test_add_input")
        inp = Input("inputLabel", String())
        w._add_item(inp)
        self.assertEqual(len(w._inputs), 1)
        self.assertEqual(w._inputs[0].input, inp)
        self.assertIsNotNone(w._nodes[inp.id()])

    def test_add_step(self):
        w = Workflow("test_add_input")
        step = Step("catStep", SingleTestTool())
        w._add_item(step)
        self.assertEqual(len(w._steps), 1)
        self.assertEqual(w._steps[0].step, step)
        self.assertIsNotNone(w._nodes[step.id()])

    def test_add_output(self):
        w = Workflow("test_add_input")
        outp = Output("outputStep", String())
        w._add_item(outp)
        self.assertEqual(len(w._outputs), 1)
        self.assertEqual(w._outputs[0].output, outp)
        self.assertIsNotNone(w._nodes[outp.id()])

    def test_add_node(self):
        w = Workflow("test_add_node")
        inp = Input("inp", String())
        stp = Step("stp", SingleTestTool())
        w.add_items([inp, stp])
        self.assertEqual(len(w._inputs), 1)
        self.assertEqual(len(w._steps), 1)
        self.assertEqual(len(w._outputs), 0)
        self.assertEqual(w._nodes[inp.id()].id(), inp.id())
        self.assertEqual(w._nodes[stp.id()].id(), stp.id())

    def test_add_qualified_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", SingleTestTool())  # Only has one input, with no output
        e = w.add_edge(inp, stp.inputs)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_add_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", SingleTestTool())  # Only has one input, with no output
        w.add_items([inp, stp])
        e = w.add_edge(inp, stp)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_anonymous_add_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", SingleTestTool())  # Only has one input, with no output
        # w.add_items([inp, stp])
        e = w.add_edge(inp, stp)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_anonymous_add_qualified_edge(self):
        w = Workflow("test_add_edge")
        inp = Input("inp", String())
        stp = Step("stp", SingleTestTool())  # Only has one input, with no output
        e = w.add_edge(inp, stp.inputs)

        self.assertEqual(e.start.id(), inp.id())
        self.assertEqual(e.finish.id(), stp.id())
        self.assertIsNone(e.stag)
        self.assertEqual(e.ftag, next(iter(stp.requires())))

    def test_pipe(self):
        w = Workflow("test_add_edge")
        inp = Input("tarred", File())
        stp = Step("stp", SingleTestTool())  # Only has one input, with no output
        out = Output("outp", Array(File()))

        w.add_pipe(inp, stp, out)

        # the nodes are usually internal
        inp_node = w._nodes[inp.id()]
        stp_node = w._nodes[stp.id()]
        out_node = w._nodes[out.id()]

        self.assertEqual(len(inp_node.connection_map), 0)
        self.assertEqual(len(stp_node.connection_map), 1)
        self.assertEqual(len(out_node.connection_map), 1)

        s1: StepInput = first_value(stp_node.connection_map)
        s2: StepInput = first_value(out_node.connection_map)

        e1 = first_value(s1.source_map)
        e2 = first_value(s2.source_map)

        self.assertEqual(e1.start.id(), inp.id())
        self.assertEqual(e1.finish.id(), stp.id())
        self.assertEqual(e2.start.id(), stp.id())
        self.assertEqual(e2.finish.id(), out.id())

    def test_qualified_pipe(self):
        w = Workflow("test_add_edge")
        inp = Input("tarred", File())
        stp = Step("stp", SingleTestTool())  # Only has one input, with no output
        out = Output("outp", Array(File()))

        w.add_pipe(inp, stp.inputs, out)

        # the nodes are usually internal
        inp_node = w._nodes[inp.id()]
        stp_node = w._nodes[stp.id()]
        out_node = w._nodes[out.id()]

        self.assertEqual(len(inp_node.connection_map), 0)
        self.assertEqual(len(stp_node.connection_map), 1)
        self.assertEqual(len(out_node.connection_map), 1)

        s1: StepInput = first_value(stp_node.connection_map)
        s2: StepInput = first_value(out_node.connection_map)

        e1: Edge = first_value(s1.source_map)
        e2: Edge = first_value(s2.source_map)

        self.assertEqual(e1.start.id(), inp.id())
        self.assertEqual(e1.finish.id(), stp.id())
        self.assertEqual(e2.start.id(), stp.id())
        self.assertEqual(e2.finish.id(), out.id())

    def test_subworkflow(self):

        w = Workflow("test_subworkflow")

        sub_w = Workflow("subworkflow")
        sub_inp = Input("sub_inp", File())
        sub_stp = Step("sub_stp", SingleTestTool())
        sub_out = Output("sub_out", Array(File()))
        sub_w.add_pipe(sub_inp, sub_stp, sub_out)

        inp = Input("inp", File())
        stp = Step("stp_workflow", sub_w)
        out = Output("out", Array(File()))
        w.add_items([inp, stp, out])
        w.add_pipe(inp, stp, out)

        # w.dump_cwl(to_disk=True)

        self.assertTrue(True)

    def test_add_scatter(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", Array(String()))
        step = Step("stp", SingleTestTool())

        e = w.add_edge(inp1, step)
        self.assertTrue(e.scatter)

    def test_add_scatter_nested_arrays(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", Array(Array(String())))
        step = Step("stp", ArrayTestTool())

        e = w.add_edge(inp1, step)
        self.assertTrue(e.scatter)

    def test_add_non_scatter(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", String())
        step = Step("stp", SingleTestTool())

        e = w.add_edge(inp1, step)
        self.assertFalse(e.scatter)

    def test_add_non_scatter2(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", Array(String()))
        step = Step("stp", ArrayTestTool())

        e = w.add_edge(inp1, step)
        self.assertFalse(e.scatter)

    def test_merge(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", Array(String()))
        step1 = Step("scatteredStp1", SingleTestTool())
        step2 = Step("mergeStp2", ArrayTestTool())

        e1 = w.add_edge(inp1, step1)
        e2 = w.add_edge(step1.out, step2)

        self.assertTrue(e1.scatter)
        self.assertFalse(e2.scatter)

    def test_add_rescatter_scattered(self):
        w = Workflow("scatterededge")

        inp1 = Input("inp1", Array(String()))
        step1 = Step("stp1", SingleTestTool())
        step2 = Step("stp2", SingleTestTool())
        e1 = w.add_edge(inp1, step1)
        e2 = w.add_edge(step1.out, step2)

        self.assertTrue(e1.scatter)
        self.assertTrue(e2.scatter)

    def test_add_single_to_array_edge(self):
        w = Workflow("test_add_single_to_array_edge")
        inp1 = Input("inp1", String())
        step1 = Step("stp1", ArrayTestTool())

        e = w.add_edge(inp1, step1.inputs)
        self.assertTrue(w.has_multiple_inputs)
