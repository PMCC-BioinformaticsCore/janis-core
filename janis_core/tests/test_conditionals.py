import unittest

from janis_core.translations import WdlTranslator
from janis_core.utils import is_module_available
from janis_core.workflow.workflow import WorkflowBuilder


@unittest.skipUnless(is_module_available("janis_unix"), "janis_unix is not available")
class TestConditionals(unittest.TestCase):
    def test_1(self):
        from janis_unix import Echo, Cat

        w = WorkflowBuilder("conditionalTest")

        w.input("inp", int, value=1)
        w.input("name", str, value="Michael")

        w.step("echo", Echo(inp=w.name), when=w.inp > 1)
        w.step("cat", Cat(file=w.echo.out), when=w.echo.out.equals("Hello, Michael"))

        w.output("out", source=w.echo.out)

        w.translate(
            "wdl"  # to_disk=True, export_path="~/Desktop/tmp/{name}", validate=True
        )

    def test_switch(self):
        from janis_unix import Echo, Cat

        w = WorkflowBuilder("switchTest")

        w.input("inp", int, value=2)
        w.input("inp1", str, value="Hello")
        w.input("inp2", str, value="Hi there")

        w.conditional(
            "echoswitch",
            [
                (w.inp > 1, Echo(inp=w.inp1)),
                (w.inp.equals(1), Echo(inp="tasy case")),
                Echo(inp=w.inp2),
            ],
        )

        w.output("out", source=w.echoswitch)

        _, wdl_tools = WdlTranslator.translate_workflow(w)
        expected = """\
version development

import "echo_v1_0_0.wdl" as E

workflow echoswitch {
  input {
    Int cond_inp
    String switch_case_1_inp
    String? switch_case_2_inp = "tasy case"
    String switch_case_3_inp
  }
  if ((cond_inp > 1)) {
     call E.echo as switch_case_1 {
      input:
        inp=switch_case_1_inp
    }
  }
  if (((cond_inp == 1) && !((cond_inp > 1)))) {
     call E.echo as switch_case_2 {
      input:
        inp=select_first([switch_case_2_inp, "tasy case"])
    }
  }
  if (!(((cond_inp > 1) || (cond_inp == 1)))) {
     call E.echo as switch_case_3 {
      input:
        inp=switch_case_3_inp
    }
  }
  output {
    File out = select_first([switch_case_1.out, switch_case_2.out, switch_case_3.out])
  }
}"""

        echoswitch = wdl_tools["echoswitch"].get_string()
        self.assertEqual(expected, echoswitch)

        # print(wdl_wf)
