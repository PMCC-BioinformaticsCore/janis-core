
import unittest
from typing import Any, Tuple

from janis_core import (
    InputSelector,
    BasenameOperator,
    FileSizeOperator,
    ReadContents,
    WildcardSelector,
    StringFormatter
)

from janis_core.ingestion.cwl.loading import load_cwl_document
from janis_core.ingestion.cwl.loading import load_cwl_version
from janis_core.ingestion.cwl.loading import load_cwlgen_from_version
from janis_core.ingestion.cwl.preprocessing import handle_inline_cltool_identifiers

from janis_core.ingestion.cwl.identifiers import get_cwl_reference
from janis_core.ingestion.cwl.expressions import parse_basic_expression

from janis_core.ingestion.cwl.types import ingest_cwl_type
from janis_core.ingestion import ingest
from janis_core.types import File
from janis_core import settings

from janis_core.ingestion.cwl.parsing.tool import CLTArgumentParser
from janis_core.ingestion.cwl.parsing.tool import CLTInputParser
from janis_core.ingestion.cwl.parsing.tool import CLTOutputParser
from janis_core.ingestion.cwl.parsing.tool import CLTToolParser
from janis_core.ingestion.cwl.parsing.workflow import WorkflowInputParser
from janis_core.ingestion.cwl.parsing.workflow import WorkflowOutputParser
from janis_core.ingestion.cwl.parsing.workflow import WorkflowStepInputsParser
from janis_core.ingestion.cwl.parsing.workflow import WorkflowStepScatterParser
from janis_core.ingestion.cwl.parsing.common import RequirementsParser

from janis_core.types import GenericFileWithSecondaries

from janis_core.messages import get_messages


class TestDatatypeErrorHandling(unittest.TestCase):
    
    def test_enum_type(self):
        from cwl_utils.parser import cwl_v1_0 as cwlutils
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        datatype = 'file:///home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/gatk_haplotype_tool.cwl#annotation_type'
        dtype, error_messages = ingest_cwl_type(datatype, cwlutils, secondary_files=None)
        self.assertIsInstance(dtype, File)
        self.assertGreater(len(error_messages), 0)

    def test_enum_type_alt(self):
        filepath = 'file:///home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/gatk_haplotype_tool.cwl'
        tool = ingest(filepath, 'cwl')



def _load_cwl_tool(filepath: str) -> Tuple[Any, Any]:
    cwl_version = load_cwl_version(filepath)
    cwl_utils = load_cwlgen_from_version(cwl_version)
    clt = load_cwl_document(filepath)
    return clt, cwl_utils


class TestJavascriptExpressionErrorHandling(unittest.TestCase):

    def test_resource_requirements(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/resources.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        
        parser = RequirementsParser(cwl_utils)
        requirements = parser.parse(clt, janis_uuid='temp')

        expected_time = '<js>[inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0]</js>'
        self.assertEqual(requirements['time'], expected_time)
        
        expected_cpus = '<js>[inputs.runtime_cpu, (2 * inputs.outputFiles), 1].filter(function (inner) { return inner != null })[0]</js>'
        self.assertEqual(requirements['cpus'], expected_cpus)
        
        expected_mem = '<js>Math.round((953.674 * [inputs.runtime_memory, ((inputs.inputFile.size / 1048576) > 1024) ? 4 : 2, 4].filter(function (inner) { return inner != null })[0]))</js>'
        self.assertEqual(requirements['memory'], expected_mem)

    def test_initial_workdir_requirement_listing(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)

        parser = CLTToolParser(cwl_utils)
        jtool = parser.parse(clt)
        self.assertIn('unnamed_1', jtool.files_to_create)
        self.assertEqual(jtool.files_to_create['unnamed_1'], '<js>${    return [{"class": "Directory",            "basename": "subdir",            "listing": [ inputs.example ]            }]}</js>')
        
        error_msgs = get_messages(jtool.uuid)
        expected_msg = "untranslated javascript expression in file / directory localisation value"
        self.assertIn(expected_msg, error_msgs)

    def test_initial_workdir_requirement_dirent_entryname(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)

        parser = CLTToolParser(cwl_utils)
        jtool = parser.parse(clt)
        self.assertIn('<js>1 + 2</js>', jtool.files_to_create)

        error_msgs = get_messages(jtool.uuid)
        expected_msg = "untranslated javascript expression in file / directory localisation name"
        self.assertIn(expected_msg, error_msgs)
    
    def test_initial_workdir_requirement_dirent_entry(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)

        parser = CLTToolParser(cwl_utils)
        jtool = parser.parse(clt)
        self.assertEqual(jtool.files_to_create['<js>1 + 2</js>'], 'workdir')

        error_msgs = get_messages(jtool.uuid)
        expected_msg = "untranslated javascript expression in file / directory localisation value"
        self.assertIn(expected_msg, error_msgs)

    def test_envvar_requirement_envdef_envvalue(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)

        parser = CLTToolParser(cwl_utils)
        jtool = parser.parse(clt)
        self.assertEqual(jtool.env_vars['test1'], '<js>9 + 10</js>')

    def test_clt_stdin(self):
        raise NotImplementedError
    
    def test_clt_stdout(self):
        raise NotImplementedError
    
    def test_clt_stderr(self):
        raise NotImplementedError
    
    def test_clt_argument_valuefrom(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        parser = CLTArgumentParser(cwl_utils)

        arg = parser.parse(clt.arguments[0])
        expected_value = '<js>runtime.outdir</js>'
        self.assertEqual(arg.value, expected_value)
        
        arg = parser.parse(clt.arguments[1])
        self.assertIsInstance(arg.value, BasenameOperator)
        self.assertEquals(arg.value.args[0].input_to_select, 'inFile')
        
        arg = parser.parse(clt.arguments[2])
        expected_value = '--noextract'
        self.assertEqual(arg.value, expected_value)
        
        arg = parser.parse(clt.arguments[3])
        expected_value = '<js>1+1</js>'
        self.assertEqual(arg.value, expected_value)
        
        arg = parser.parse(clt.arguments[4])
        expected_value = '<js>"/foo/bar/baz".split(\'/\').slice(-1)[0]</js>'
        self.assertEqual(arg.value, expected_value)
        
        arg = parser.parse(clt.arguments[5])
        expected_value = '<js>${  var r = [];  for (var i = 10; i >= 1; i--) {    r.push(i);  }  return r;}</js>'
        self.assertEqual(arg.value, expected_value)
    
    @unittest.skip("unnecessary")
    def test_clt_input_format(self):
        # don't need to worry about File format in cwl - type is ingested as File
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        parser = CLTInputParser(cwl_utils)
        tinput = parser.parse(clt.inputs[3])
    
    def test_clt_input_secondaryfiles(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        parser = CLTInputParser(cwl_utils)
        tinput = parser.parse(clt.inputs[1])
        self.assertIsInstance(tinput.input_type, GenericFileWithSecondaries)
        error_msgs = get_messages(tinput.uuid)
        expected_msg = "could not parse secondaries format from javascript expression: <js>self.basename + self.nameext.replace('m','i')</js>"
        self.assertIn(expected_msg, error_msgs)

    def test_clt_input_inputbinding_valuefrom(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        
        parser = CLTInputParser(cwl_utils)
        tinput = parser.parse(clt.inputs[5])
        expected_value = '<js>[inputs.runtime_cpu, 16, 1].filter(function (inner) { return inner != null })[0]</js>'
        self.assertEqual(tinput.value, expected_value)
    
    @unittest.skip("unnecessary")
    def test_clt_output_format(self):
        # don't need to worry about File format in cwl - type is ingested as File
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        parser = CLTOutputParser(cwl_utils)
        tout = parser.parse(clt.outputs[1])

    def test_clt_output_secondaryfiles(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        parser = CLTOutputParser(cwl_utils)
        tout = parser.parse(clt.outputs[1])
        self.assertIsInstance(tout.output_type, GenericFileWithSecondaries)
        error_msgs = get_messages(tout.uuid)
        expected_msg = "could not parse secondaries format from javascript expression: <js>inputs.bam.basename + inputs.bam.nameext.replace('m','i')</js>"
        self.assertIn(expected_msg, error_msgs)
    
    def test_clt_output_outputbinding_glob(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        parser = CLTOutputParser(cwl_utils)
        tout = parser.parse(clt.outputs[0])
        self.assertIsInstance(tout.selector, WildcardSelector)
        self.assertEqual(tout.selector.wildcard, 'output.txt')

        tout = parser.parse(clt.outputs[2])
        self.assertIsInstance(tout.selector.wildcard, StringFormatter)
        self.assertEqual(tout.selector.wildcard._format, '{JANIS_CWL_TOKEN_1}.trimmed.fastq')
        self.assertEqual(tout.selector.wildcard.kwargs['JANIS_CWL_TOKEN_1'], "<js>inputs.fastq.nameroot.replace(/\\b.fastq\\b/g, '')</js>")
    
    def test_clt_output_outputbinding_outputEval(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_tool(filepath)
        parser = CLTOutputParser(cwl_utils)
        tout = parser.parse(clt.outputs[3])
        self.assertIsInstance(tout.selector, ReadContents)
        self.assertEqual(tout.selector.args[0], '<js>self[0]</js>')
    
    @unittest.skip("unnecessary")
    def test_wf_input_format(self) -> None:
        # don't need to worry about File format in cwl - type is ingested as File
        raise NotImplementedError

    def test_wf_input_secondaryfiles(self):
        raise NotImplementedError
    
    @unittest.skip("unnecessary")
    def test_wf_output_format(self) -> None:
        # don't need to worry about File format in cwl - type is ingested as File
        raise NotImplementedError

    def test_wf_output_secondaryfiles(self):
        raise NotImplementedError
  
    def test_step_input_valuefrom(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/super_enhancer_wf.cwl'
        wf = ingest(filepath, 'cwl')
        raise NotImplementedError







class TestIdentifiers(unittest.TestCase):

    def test_filename_id(self):
        ident = 'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/workflows/super_enhancer_wf.cwl'
        ref = get_cwl_reference(ident)
        self.assertEqual(ref.filename, 'super_enhancer_wf')
        self.assertEqual(ref.internal_path, None)
        self.assertEqual(ref.entity, None)
    
    def test_entity_id(self):
        ident = 'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/workflows/super_enhancer_wf.cwl#annotation_file'
        ref = get_cwl_reference(ident)
        self.assertEqual(ref.filename, 'super_enhancer_wf')
        self.assertEqual(ref.internal_path, None)
        self.assertEqual(ref.entity, 'annotation_file')
    
    def test_internal_path_id(self):
        ident = 'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/workflows/super_enhancer_wf.cwl#assign_genes/result_file'
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
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
        workflow = handle_inline_cltool_identifiers(workflow)
        for step in workflow.steps:
            tool_id = step.run.id
            self.assertTrue(tool_id.endswith('_tool.cwl'))



class TestDocumentLoading(unittest.TestCase):

    def test_load_cwl_document(self):
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
        self.assertIsNotNone(workflow)
    
    def test_load_cwl_version(self):
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/structuralvariants/workflow.cwl'
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
