
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

from janis_core.ingestion.cwl.loading import load_cwl_document
from janis_core.ingestion.cwl.loading import load_cwl_version
from janis_core.ingestion.cwl.preprocessing import handle_inline_cltool_identifiers

from janis_core.ingestion.cwl.identifiers import get_cwl_reference
from janis_core.ingestion.cwl.parsing import parse_basic_expression


class TestIdentifiers(unittest.TestCase):

    def test_filename_id(self):
        ident = 'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl'
        ref = get_cwl_reference(ident)
        self.assertEqual(ref.filename, 'super_enhancer_wf')
        self.assertEqual(ref.internal_path, None)
        self.assertEqual(ref.entity, None)
    
    def test_entity_id(self):
        ident = 'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl#annotation_file'
        ref = get_cwl_reference(ident)
        self.assertEqual(ref.filename, 'super_enhancer_wf')
        self.assertEqual(ref.internal_path, None)
        self.assertEqual(ref.entity, 'annotation_file')
    
    def test_internal_path_id(self):
        ident = 'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl#assign_genes/result_file'
        ref = get_cwl_reference(ident)
        self.assertEqual(ref.filename, 'super_enhancer_wf')
        self.assertEqual(ref.internal_path, 'assign_genes')
        self.assertEqual(ref.entity, 'result_file')
    
    def test_inline_step_id(self):
        ident = 'add_island_names_tool.cwl'
        ref = get_cwl_reference(ident)
        self.assertEqual(ref.filename, 'add_island_names_tool')
        self.assertEqual(ref.internal_path, None)
        self.assertEqual(ref.entity, None)
    

            
            
class TestPreprocessing(unittest.TestCase):

    def test_inline_cltool_ids(self):
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
        workflow = handle_inline_cltool_identifiers(workflow)
        for step in workflow.steps:
            tool_id = step.run.id
            self.assertTrue(tool_id.endswith('_tool.cwl'))





class TestDocumentLoading(unittest.TestCase):

    def test_load_cwl_document(self):
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
        self.assertIsNotNone(workflow)
    
    def test_load_cwl_version(self):
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/structuralvariants/workflow.cwl'
        version = load_cwl_version(doc)
        self.assertIsNotNone(version)
        


class TestFromCwlExpressions(unittest.TestCase):

    def test_number(self):
        expr = "${return 16 }"
        result = parse_basic_expression(expr)
        self.assertEqual(16, result)

    def test_number_with_semicolon(self):
        expr = "${return 16;}"
        result = parse_basic_expression(expr)
        self.assertEqual(16, result)

    def test_number_with_spaces(self):
        expr = "${ return 80000 }"
        result = parse_basic_expression(expr)
        self.assertEqual(80000, result)

    def test_string(self):
        expr = '${ return "over 80000" }'
        result = parse_basic_expression(expr)
        self.assertEqual("over 80000", result)

    def test_input_selector(self):
        expr = "$(inputs.my_input)"
        result = parse_basic_expression(expr)
        self.assertIsInstance(result, InputSelector)
        self.assertEqual("my_input", result.input_to_select)

    def test_input_selector_with_basename(self):
        expr = "$(inputs.my_input.basename)"
        result = parse_basic_expression(expr)
        self.assertIsInstance(result, BasenameOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)

    def test_input_selector_with_filesize(self):
        expr = "$(inputs.my_input.size)"
        result = parse_basic_expression(expr)
        self.assertIsInstance(result, FileSizeOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)

    def test_input_selector_with_contents(self):
        expr = "$(inputs.my_input.contents)"
        result = parse_basic_expression(expr)
        self.assertIsInstance(result, ReadContents)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)
