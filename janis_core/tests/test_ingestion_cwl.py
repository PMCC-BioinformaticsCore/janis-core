
import unittest
import os 

from janis_core import (
    InputSelector,
    BasenameOperator,
    FileSizeOperator,
    ReadContents,
)

from janis_core import (
    InputSelector,
    BasenameOperator,
    FileSizeOperator,
    ReadContents,
) 

from janis_core.ingestion.cwl import CWlParser
from janis_core.ingestion.cwl.loading import load_cwl_document
from janis_core.ingestion.cwl.preprocessing import handle_inline_cltool_identifiers


class TestPreprocessing(unittest.TestCase):

    def test_inline_cltool_ids(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/super_enhancer_wf.cwl'
        workflow = load_cwl_document(filepath)
        workflow = handle_inline_cltool_identifiers(workflow)
        for step in workflow.steps:
            tool_id = step.run.id
            self.assertTrue(tool_id.endswith('_tool.cwl'))





class TestDocumentLoading(unittest.TestCase):

    def test_workflow1(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/super_enhancer_wf.cwl'
        workflow = load_cwl_document(filepath)
        self.assertIsNotNone(workflow)
    
    def test_workflow2(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/structuralvariants/workflow.cwl'
        workflow = load_cwl_document(filepath)
        self.assertIsNotNone(workflow)

    def test_workflow_subworkflow(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/structuralvariants/workflow.cwl'
        dirpath = os.path.dirname(filepath)  # abs path to main file
        workflow = load_cwl_document(filepath)
        step = workflow.steps[0]
        relpath = step.run # relative path to tool
        subwf = load_cwl_document(relpath, dirpath)  
        self.assertIsNotNone(subwf)
        
    def test_workflow_tool(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/structuralvariants/workflow.cwl'
        dirpath = os.path.dirname(filepath)  # abs path to main file
        workflow = load_cwl_document(filepath)
        step = workflow.steps[9]
        relpath = step.run # relative path to tool
        tool = load_cwl_document(relpath, dirpath)  
        self.assertIsNotNone(tool)


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
