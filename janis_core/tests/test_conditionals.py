import unittest

from janis_core.translations import WdlTranslator
from janis_core.workflow.workflow import WorkflowBuilder

from janis_core.redefinitions.tools import Echo, Cat
from janis_core.translations.common import to_builders


class TestConditionals(unittest.TestCase):
    def test_1(self):
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
        w = to_builders(w)

        translator = WdlTranslator()
        translator.translate_workflow_internal(w)
        assert translator.main is not None
        wdltool = translator.main[1]
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
        inp=switch_case_1_inp   # String    
    }
  }

  if (((cond_inp == 1) && !((cond_inp > 1)))) {
     call E.echo as switch_case_2 {
      input:
        inp=select_first([switch_case_2_inp, "tasy case"])   # String    
    }
  }

  if (!(((cond_inp > 1) || (cond_inp == 1)))) {
     call E.echo as switch_case_3 {
      input:
        inp=switch_case_3_inp   # String    
    }
  }

  output {
    File out = select_first([switch_case_1.out, switch_case_2.out, switch_case_3.out])
  }

}"""

        echoswitch = wdltool.get_string()
        print(echoswitch)
        self.assertEqual(expected, echoswitch)

