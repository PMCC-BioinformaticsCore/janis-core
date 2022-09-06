import unittest
import os 

from janis_core import (
    InputSelector,
    BasenameOperator,
    FileSizeOperator,
    ReadContents,
)

from janis_core.ingestion.fromcwl import CWlParser
from janis_core.ingestion.fromwdl import WdlParser
from janis_core.ingestion.main import ingest_galaxy


class TestFromWdl(unittest.TestCase):
    parser = WdlParser()

    def test_ingest_tool(self) -> None:
        raise NotImplementedError

    def test_ingest_workflow(self) -> None:
        raise NotImplementedError



class TestFromGalaxy(unittest.TestCase):
    TOOL_PATH = os.path.abspath('./janis_core/tests/data/galaxy/abricate/abricate.xml')
    WORKFLOW_PATH = os.path.abspath('./janis_core/tests/data/galaxy/assembly.ga')

    def test_ingest_tool(self) -> None:
        jtool = ingest_galaxy(self.TOOL_PATH)
        self.assertEquals(len(jtool.inputs()), 5)
        self.assertEquals(len(jtool.outputs()), 1)
        self.assertEquals(jtool.base_command(), ['abricate'])

    def test_ingest_workflow(self) -> None:
        jworkflow = ingest_galaxy(self.WORKFLOW_PATH)
        self.assertEquals(len(jworkflow.step_nodes), 6)
        self.assertEquals(len(jworkflow.output_nodes), 19)
        self.assertIn('inForwardReads', jworkflow.input_nodes)
        self.assertIn('inReverseReads', jworkflow.input_nodes)
        self.assertIn('inLongReads', jworkflow.input_nodes)



class TestFromCwlExpressions(unittest.TestCase):
    parser = CWlParser(cwl_version="v1.0")

    def test_number(self):
        expr = "${return 16 }"
        result = self.parser.parse_basic_expression(expr)
        self.assertEqual(16, result)

    def test_number_with_semicolon(self):
        expr = "${return 16;}"
        result = self.parser.parse_basic_expression(expr)
        self.assertEqual(16, result)

    def test_number_with_spaces(self):
        expr = "${ return 80000 }"
        result = self.parser.parse_basic_expression(expr)
        self.assertEqual(80000, result)

    def test_string(self):
        expr = '${ return "over 80000" }'
        result = self.parser.parse_basic_expression(expr)
        self.assertEqual("over 80000", result)

    def test_input_selector(self):
        expr = "$(inputs.my_input)"
        result = self.parser.parse_basic_expression(expr)
        self.assertIsInstance(result, InputSelector)
        self.assertEqual("my_input", result.input_to_select)

    def test_input_selector_with_basename(self):
        expr = "$(inputs.my_input.basename)"
        result = self.parser.parse_basic_expression(expr)
        self.assertIsInstance(result, BasenameOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)

    def test_input_selector_with_filesize(self):
        expr = "$(inputs.my_input.size)"
        result = self.parser.parse_basic_expression(expr)
        self.assertIsInstance(result, FileSizeOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)

    def test_input_selector_with_contents(self):
        expr = "$(inputs.my_input.contents)"
        result = self.parser.parse_basic_expression(expr)
        self.assertIsInstance(result, ReadContents)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)
