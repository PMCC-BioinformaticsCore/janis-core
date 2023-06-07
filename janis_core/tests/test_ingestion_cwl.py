
import unittest
from typing import Any, Tuple

from janis_core import (
    InputSelector,
    BasenameOperator,
    FileSizeOperator,
    ReadContents,
    WildcardSelector,
    StringFormatter,
    InputNodeSelector,
    NamerootOperator
)

from janis_core.ingestion.cwl.loading import load_cwl_document
from janis_core.ingestion.cwl.loading import load_cwl_version
from janis_core.ingestion.cwl.loading import load_cwl_utils_from_version

from janis_core.ingestion.cwl.identifiers import get_cwl_reference
from janis_core.ingestion.cwl.expressions import parse_basic_expression

from janis_core.ingestion.cwl.parsing.tool import CLTArgumentParser
from janis_core.ingestion.cwl.parsing.tool import CLTInputParser
from janis_core.ingestion.cwl.parsing.tool import CLTOutputParser
from janis_core.ingestion.cwl.parsing.tool import CLTParser
from janis_core.ingestion.cwl.parsing.tool import CLTRequirementsParser
from janis_core.ingestion.cwl.types import ingest_cwl_type
from janis_core.ingestion.cwl import parse as parse_cwl
from janis_core.ingestion.cwl import CWlParser

from janis_core.types import File
from janis_core.types import GenericFileWithSecondaries
from janis_core.messages import get_messages
from janis_core import settings




def _load_cwl_doc(filepath: str) -> Tuple[Any, Any]:
    cwl_version = load_cwl_version(filepath)
    cwl_utils = load_cwl_utils_from_version(cwl_version)
    loaded_doc = load_cwl_document(filepath)
    return loaded_doc, cwl_utils



class TestRequirementsParsing(unittest.TestCase):

    def test_directories_to_create1(self):
        pass

    def test_directories_to_create2(self):
        pass
        
    def test_directories_to_create3(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/prepare_fasta_db.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, uuid='test_directories_to_create3')
        reqs = parser.parse()
        msgs = get_messages('test_directories_to_create3')
        print()
 
    def test_files_to_create1(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/bundle_secondaryfiles.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt)
        reqs = parser.do_parse()
        self.assertDictEqual(reqs['files_to_create'], {})

    def test_files_to_create2(self):
        # expression tool
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/return_directory.cwl'
        etool, cwl_utils = _load_cwl_doc(filepath)
        
        parser = CWlParser(filepath)
        cltool = parser.parse_etool_to_cltool(etool)
        
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=cltool, entity=cltool, is_expression_tool=True)
        reqs = parser.do_parse()
        self.assertIn('return_directory.js', reqs['files_to_create'])
        self.assertIsInstance(reqs['files_to_create']['return_directory.js'], str)
    
    def test_files_to_create8(self):
        # expression tool
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressiontools/check_value.cwl'
        etool, cwl_utils = _load_cwl_doc(filepath)
        
        parser = CWlParser(filepath)
        cltool = parser.parse_etool_to_cltool(etool)
        
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=cltool, entity=cltool, is_expression_tool=True)
        reqs = parser.do_parse()

        # check js script file created, and is a string
        self.assertIn('check_value.js', reqs['files_to_create'])
        self.assertIsInstance(reqs['files_to_create']['check_value.js'], str)

        # check js script file modified correctly
        script = reqs['files_to_create']['check_value.js']
        self.assertNotIn('"use strict";\nvar inputs=$(inputs);\nvar runtime=$(runtime);', script)
        self.assertIn('var inputs = JSON.parse( process.argv[2] );\n', script)
 
    def test_files_to_create3(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/check_threshold.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt)
        reqs = parser.do_parse()
        self.assertIn('check_threshold.py', reqs['files_to_create'])
        self.assertIsInstance(reqs['files_to_create']['check_threshold.py'], str)
        self.assertIn('#!/usr/bin/env python3', reqs['files_to_create']['check_threshold.py'])
 
    def test_files_to_create4(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/amplicon_architect.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, uuid='test_files_to_create4')
        reqs = parser.parse()
        msgs = get_messages('test_files_to_create4')
        self.assertIn('likely untranslated cwl / js in script: setup_vars.sh', msgs)
    
    def test_files_to_create5(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/control-freec-11-6-sbg.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, uuid='test_files_to_create5')
        reqs = parser.parse()
        self.assertIn('config.txt', reqs['files_to_create'])
        msgs = get_messages('test_files_to_create5')
        self.assertIn('config.txt: js code to dynamically create runtime file. please address', msgs)
        print()
    
    def test_files_to_create6(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/antismash_v4.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, uuid='test_files_to_create6')
        reqs = parser.parse()
        msgs = get_messages('test_files_to_create6')
        print()
    
    def test_files_to_create7(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/add_crossmapped2resmapped.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, uuid='test_files_to_create7')
        reqs = parser.parse()
        msgs = get_messages('test_files_to_create7')
        print()
    
    def test_env_vars1(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/requirements/amplicon_architect.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt)
        reqs = parser.do_parse()
        self.assertIn('AA_SRC', reqs['env_vars'])
        self.assertEqual(reqs['env_vars']['AA_SRC'], '/home/programs/AmpliconArchitect-master/src')




class TestFallbacksErrorHandling(unittest.TestCase):
    @unittest.skip('not implemented')
    def test_workflow_input_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_workflow_step_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_workflow_step_scatter_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_workflow_step_input_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_workflow_output_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')    
    def test_clt_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_clt_argument_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')    
    def test_clt_input_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')    
    def test_clt_output_parser(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')    
    def test_requirements_parser(self) -> None:
        raise NotImplementedError





class TestDatatypeErrorHandling(unittest.TestCase):

    def test_multiple_datatypes(self):
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        from cwl_utils.parser import cwl_v1_0 as cwl_utils
        cwl_type = ['null', 'File', 'string']
        dtype, error_messages = ingest_cwl_type(cwl_type, cwl_utils, secondary_files=None)
        self.assertIsInstance(dtype, File)
        self.assertIn("entity supports multiple datatypes: ['File', 'String']. selected File as fallback. this may affect pipeline execution", error_messages)
    
    def test_unparseable_generic(self):
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        from cwl_utils.parser import cwl_v1_0 as cwl_utils
        datatype = 'kitten'
        dtype, error_messages = ingest_cwl_type(datatype, cwl_utils, secondary_files=None)
        self.assertIsInstance(dtype, File)
        self.assertIn('Unsupported datatype: kitten. Treated as a file.', error_messages)
    
    def test_enum_type(self):
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        from cwl_utils.parser import cwl_v1_0 as cwl_utils
        datatype = 'file:///home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/gatk_haplotype_tool.cwl#annotation_type'
        dtype, error_messages = ingest_cwl_type(datatype, cwl_utils, secondary_files=None)
        self.assertIsInstance(dtype, File)
        self.assertGreater(len(error_messages), 0)

    def test_unparseable_secondary_type(self):
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        cinp = clt.inputs[1]
        dtype, error_messages = ingest_cwl_type(cinp.type, cwl_utils, secondary_files=cinp.secondaryFiles)
        self.assertIsInstance(dtype, GenericFileWithSecondaries)
        self.assertGreater(len(error_messages), 0)




class TestJavascriptExpressionErrorHandling(unittest.TestCase):

    def test_resource_requirements(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/resources.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        
        parser = CLTRequirementsParser(cwl_utils, clt=clt, entity=clt, uuid='placeholder')
        requirements = parser.parse()

        expected_time = '<js>[inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0]</js>'
        self.assertIsInstance(requirements['time'], StringFormatter)
        self.assertEqual(requirements['time']._format, "{JANIS_CWL_TOKEN_1}")
        self.assertEqual(requirements['time'].kwargs["JANIS_CWL_TOKEN_1"], expected_time)
        
        expected_cpus = '<js>[inputs.runtime_cpu, (2 * inputs.outputFiles), 1].filter(function (inner) { return inner != null })[0]</js>'
        self.assertIsInstance(requirements['cpus'], StringFormatter)
        self.assertEqual(requirements['cpus']._format, "{JANIS_CWL_TOKEN_1}")
        self.assertEqual(requirements['cpus'].kwargs["JANIS_CWL_TOKEN_1"], expected_cpus)
        
        expected_mem = '<js>Math.round((953.674 * [inputs.runtime_memory, ((inputs.inputFile.size / 1048576) > 1024) ? 4 : 2, 4].filter(function (inner) { return inner != null })[0]))</js>'
        self.assertEqual(requirements["memory"], expected_mem)

    def test_initial_workdir_requirement_listing(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)

        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        msgs = get_messages(jtool.uuid)
        self.assertIsNone(jtool.files_to_create())
        
        expected_msg = 'directory required for runtime generated using cwl / js: <js>${    return [{"class": "Directory",            "basename": "subdir",            "listing": [ inputs.example ]            }]}</js>. please address'
        self.assertIn(expected_msg, msgs)
        
    def test_initial_workdir_requirement_dirent_entryname(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)

        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        msgs = get_messages(jtool.uuid)
        self.assertIsNone(jtool.files_to_create())

        expected_msg = 'untranslated javascript expression in environment variable value'
        self.assertIn(expected_msg, msgs)
    
    def test_initial_workdir_requirement_dirent_entry(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)

        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        msgs = get_messages(jtool.uuid)
        self.assertIsNone(jtool.files_to_create())

        expected_msg = 'untranslated javascript expression in environment variable value'
        self.assertIn(expected_msg, msgs)

    def test_envvar_requirement_envdef_envvalue(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)

        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        self.assertEqual(jtool.env_vars()['test1'], '<js>9 + 10</js>')

    def test_clt_stdout(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/streams.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        
        arg = jtool.arguments()[1]
        self.assertEqual(arg.prefix, '>')
        self.assertEqual(arg.value, 'stdout.txt')

        msgs = get_messages(jtool.uuid)
        expected_msg = 'untranslated javascript expression in stdout filename. used stdout.txt as fallback'
        self.assertIn(expected_msg, msgs)
    
    def test_clt_stderr(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/streams.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        
        arg = jtool.arguments()[0]
        self.assertEqual(arg.prefix, '2>')
        self.assertEqual(arg.value, 'stderr.txt')

        msgs = get_messages(jtool.uuid)
        expected_msg = 'untranslated javascript expression in stdout filename. used stdout.txt as fallback'
        self.assertIn(expected_msg, msgs)

    def test_clt_stdin(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/streams.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        
        arg = jtool.arguments()[2]
        self.assertEqual(arg.prefix, '<')
        self.assertIsInstance(arg.value, InputSelector)
        self.assertEqual(arg.value.input_to_select, 'inFile')

    def test_clt_argument_valuefrom(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)

        parser = CLTArgumentParser(cwl_utils, clt=clt, entity=clt.arguments[0])
        arg = parser.parse()
        expected_value = '.'
        self.assertEqual(arg.value, expected_value)
        
        parser = CLTArgumentParser(cwl_utils, clt=clt, entity=clt.arguments[1])
        arg = parser.parse()
        self.assertIsInstance(arg.value, BasenameOperator)
        self.assertEquals(arg.value.args[0].input_to_select, 'inFile')
        
        parser = CLTArgumentParser(cwl_utils, clt=clt, entity=clt.arguments[2])
        arg = parser.parse()
        expected_value = '--noextract'
        self.assertEqual(arg.value, expected_value)
        
        parser = CLTArgumentParser(cwl_utils, clt=clt, entity=clt.arguments[3])
        arg = parser.parse()
        expected_value = '<js>1+1</js>'
        self.assertEqual(arg.value, expected_value)
        
        parser = CLTArgumentParser(cwl_utils, clt=clt, entity=clt.arguments[4])
        arg = parser.parse()
        expected_value = '<js>"/foo/bar/baz".split(\'/\').slice(-1)[0]</js>'
        self.assertEqual(arg.value, expected_value)
        
        parser = CLTArgumentParser(cwl_utils, clt=clt, entity=clt.arguments[5])
        arg = parser.parse()
        expected_value = '<js>${  var r = [];  for (var i = 10; i >= 1; i--) {    r.push(i);  }  return r;}</js>'
        self.assertEqual(arg.value, expected_value)
    
    @unittest.skip("TODO leave log_info message")
    def test_clt_input_format(self):
        # don't need to worry about File format in cwl - type is ingested as File
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTInputParser(cwl_utils, clt=clt, entity=clt.inputs[3])
        tinput = parser.parse()
    
    def test_clt_input_secondaryfiles(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTInputParser(cwl_utils, clt=clt, entity=clt.inputs[1])
        tinput = parser.parse()
        self.assertIsInstance(tinput.input_type, GenericFileWithSecondaries)
        error_msgs = get_messages(tinput.uuid)
        expected_msg = "could not parse secondaries format from javascript expression: {<js>self.basename + self.nameext.replace('m','i')</js>}"
        self.assertIn(expected_msg, error_msgs)

    def test_clt_input_inputbinding_valuefrom(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        
        parser = CLTInputParser(cwl_utils, clt=clt, entity=clt.inputs[5])
        tinput = parser.parse()
        expected_value = '<js>[inputs.runtime_cpu, 16, 1].filter(function (inner) { return inner != null })[0]</js>'
        self.assertEqual(tinput.value, expected_value)
    
    @unittest.skip("TODO leave log_info message")
    def test_clt_output_format(self):
        # don't need to worry about File format in cwl - type is ingested as File
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[1])
        tout = parser.parse()

    def test_clt_output_secondaryfiles(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[1])
        tout = parser.parse()
        self.assertIsInstance(tout.output_type, GenericFileWithSecondaries)
        error_msgs = get_messages(tout.uuid)
        expected_msg = "could not parse secondaries format from javascript expression: <js>inputs.bam.basename + inputs.bam.nameext.replace('m','i')</js>"
        self.assertIn(expected_msg, error_msgs)
    
    def test_clt_output_outputbinding_glob(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[0])
        tout = parser.parse()
        self.assertIsInstance(tout.selector, WildcardSelector)
        self.assertEqual(tout.selector.wildcard, 'output.txt')

        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[2])
        tout = parser.parse()
        self.assertIsInstance(tout.selector.wildcard, StringFormatter)
        self.assertEqual(tout.selector.wildcard._format, '{JANIS_CWL_TOKEN_1}.trimmed.fastq')
        self.assertEqual(tout.selector.wildcard.kwargs['JANIS_CWL_TOKEN_1'], "<js>inputs.fastq.nameroot.replace(/\\b.fastq\\b/g, '')</js>")
    
    def test_clt_output_outputbinding_outputEval(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[3])
        tout = parser.parse()
        self.assertIsInstance(tout.selector, ReadContents)
        self.assertEqual(tout.selector.args[0], '<js>self[0]</js>')
    
    @unittest.skip("TODO leave log_info message")
    def test_wf_input_format(self) -> None:
        # don't need to worry about File format in cwl - type is ingested as File
        raise NotImplementedError

    def test_wf_input_secondaryfiles(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/expressions.cwl'
        wf = parse_cwl(filepath)
        winp = wf.input_nodes['bambai_pair2']
        self.assertIsInstance(winp.datatype, GenericFileWithSecondaries)
        self.assertEqual(len(winp.datatype.secondaries), 0)
        error_msgs = get_messages(winp.uuid)
        expected_msg = "could not parse secondaries format from javascript expression: <js>self.basename + self.nameext.replace('m','i')</js>"
        self.assertIn(expected_msg, error_msgs)
    
    def test_step_input_valuefrom(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/expressions.cwl'
        wf = parse_cwl(filepath)
        
        # example 1
        jstep = wf.step_nodes['echo2']
        inp_name = 'text'
        
        self.assertIn(inp_name, jstep.tool.connections)
        self.assertEqual(jstep.tool.connections[inp_name], "<js>$(inputs.prefix_str + '_' + inputs.suffix_str)</js>")
        self.assertIn(inp_name, jstep.sources)
        self.assertEqual(jstep.sources[inp_name].source_map[0].source.input_node.default, "<js>$(inputs.prefix_str + '_' + inputs.suffix_str)</js>")
        
        msgs = get_messages(jstep.uuid)
        expected_msg = "'text' input value contains untranslated javascript expression: <js>$(inputs.prefix_str + '_' + inputs.suffix_str)</js>"
        self.assertIn(expected_msg, msgs)

        # example 2
        jstep = wf.step_nodes['rename_png1']
        inp_name = 'target_filename'
        
        self.assertIn(inp_name, jstep.tool.connections)
        self.assertIsInstance(jstep.tool.connections[inp_name], InputNodeSelector)
        self.assertEqual(jstep.tool.connections[inp_name].input_node.id(), 'bambai_pair1')
        
        self.assertIn(inp_name, jstep.sources)
        self.assertIsNone(jstep.sources[inp_name].source_map[0].source.input_node.default)
        
        msgs = get_messages(jstep.uuid)
        self.assertIn('\'target_filename\' input value contains untranslated javascript expression: <js>$(inputs.bambai_pair1.location.split(\'/\').slice(-1)[0].split(\'.\').slice(0,-1).join(\'.\')+"_default_s_enhcr.png")</js>', msgs)

    @unittest.skip("TODO leave log_info message")
    def test_wf_output_format(self) -> None:
        # TODO implement this - log info message about the format which was discarded
        # don't need to worry about File format in cwl - type is ingested as File
        raise NotImplementedError

    def test_wf_output_secondaryfiles(self):
        filepath = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/expressions.cwl'
        wf = parse_cwl(filepath)
        winp = wf.output_nodes['out3']
        self.assertIsInstance(winp.datatype, GenericFileWithSecondaries)
        self.assertEqual(len(winp.datatype.secondaries), 0)
        error_msgs = get_messages(winp.uuid)
        expected_msg = "could not parse secondaries format from javascript expression: <js>self.basename + self.nameext.replace('m','i')</js>"
        self.assertIn(expected_msg, error_msgs)
  







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

    def test_convert_cwl_types_to_python_clt(self):
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/gatk_haplotype_tool.cwl'
        clt = load_cwl_document(doc)
        # no idea how to test :/
    
    def test_convert_cwl_types_to_python_workflow(self):
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
        # no idea how to test :/
    
    def test_inline_cltool_ids(self):
        doc = '/home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/workflows/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
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
        result, success = parse_basic_expression(expr)
        self.assertEqual(16, result)

    def test_number_with_semicolon(self):
        expr = "${return 16;}"
        result, success = parse_basic_expression(expr)
        self.assertEqual(16, result)

    def test_number_with_spaces(self):
        expr = "${ return 80000 }"
        result, success = parse_basic_expression(expr)
        self.assertEqual(80000, result)

    def test_string(self):
        expr = '${ return "over 80000" }'
        result, success = parse_basic_expression(expr)
        self.assertEqual("over 80000", result)

    def test_input_selector(self):
        expr = "$(inputs.my_input)"
        result, success = parse_basic_expression(expr)
        self.assertIsInstance(result, InputSelector)
        self.assertEqual("my_input", result.input_to_select)

    def test_input_selector_with_basename(self):
        expr = "$(inputs.my_input.basename)"
        result, success = parse_basic_expression(expr)
        self.assertIsInstance(result, BasenameOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)

    def test_input_selector_with_filesize(self):
        expr = "$(inputs.my_input.size)"
        result, success = parse_basic_expression(expr)
        self.assertIsInstance(result, FileSizeOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)

    def test_input_selector_with_contents(self):
        expr = "$(inputs.my_input.contents)"
        result, success = parse_basic_expression(expr)
        self.assertIsInstance(result, ReadContents)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual("my_input", result.args[0].input_to_select)
    
    def test_runtime_outdir(self):
        expr = "$(runtime.outdir)"
        result, success = parse_basic_expression(expr)
        self.assertTrue(success)
        self.assertEquals(result, '.')
    
    def test_composite_1(self):
        expr = "$(runtime.outdir)/$(inputs.bam.basename)"
        result, success = parse_basic_expression(expr)
        self.assertTrue(success)
        self.assertIsInstance(result, BasenameOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertEqual(result.args[0].input_to_select, 'bam')
    
    def test_composite_2(self):
        expr = "$(inputs.bam.basename).bai"
        result, success = parse_basic_expression(expr)
        self.assertTrue(success)
        self.assertIsInstance(result, StringFormatter)
        self.assertEqual(result._format, "{JANIS_CWL_TOKEN_1}.bai")
        self.assertIsInstance(result.kwargs['JANIS_CWL_TOKEN_1'], BasenameOperator) 
        self.assertIsInstance(result.kwargs['JANIS_CWL_TOKEN_1'].args[0], InputSelector) 
        self.assertEqual(result.kwargs['JANIS_CWL_TOKEN_1'].args[0].input_to_select, "bam") 
    
    def test_composite_3(self):
        expr = "$(runtime.outdir)/$(inputs.bam.nameroot).bai"
        result, success = parse_basic_expression(expr)
        self.assertTrue(success)
        self.assertIsInstance(result, StringFormatter)
        self.assertEqual(result._format, "{JANIS_CWL_TOKEN_1}.bai")
        self.assertIsInstance(result.kwargs['JANIS_CWL_TOKEN_1'], NamerootOperator) 
        self.assertIsInstance(result.kwargs['JANIS_CWL_TOKEN_1'].args[0], InputSelector) 
        self.assertEqual(result.kwargs['JANIS_CWL_TOKEN_1'].args[0].input_to_select, "bam") 
