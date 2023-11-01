
import unittest
import os
from typing import Any, Tuple

from janis_core import (
    InputSelector,
    CommandToolBuilder,
    ToolInput,
)

from janis_core.operators.stringformatter import StringFormatter
from janis_core.operators.standard import (
    ReadContents,
    BasenameOperator,
    SplitOperator,
    FileSizeOperator,
    NamerootOperator,
    NameextOperator,
    DirnameOperator,
    JoinOperator,
    LengthOperator,
    FlattenOperator,
    FirstOperator,
    SliceOperator,
    ReplaceOperator,
)
from janis_core.operators.operator import (
    IndexOperator,
)
from janis_core.operators.selectors import (
    WildcardSelector,
    InputNodeSelector,
    CpuSelector,
    MemorySelector,  
)
from janis_core.operators.logical import (
    IsDefined,
    If,
    AssertNotNull,
    GroupOperator,
    NotOperator,
    AndOperator,
    OrOperator,
    EqualityOperator,
    InequalityOperator,
    GtOperator,
    GteOperator,
    LtOperator,
    LteOperator,
    AddOperator,
    SubtractOperator,
    MultiplyOperator,
    DivideOperator,
    FloorOperator,
    CeilOperator,
    RoundOperator,
)

from janis_core.ingestion.cwl.loading import load_cwl_document
from janis_core.ingestion.cwl.loading import load_cwl_version
from janis_core.ingestion.cwl.loading import load_cwl_utils_from_version

from janis_core.ingestion.cwl.identifiers import get_cwl_reference
from janis_core.ingestion.cwl.expressions import parse_expression

from janis_core.ingestion.cwl.parsing.tool import CLTArgumentParser
from janis_core.ingestion.cwl.parsing.tool import CLTInputParser
from janis_core.ingestion.cwl.parsing.tool import CLTOutputParser
from janis_core.ingestion.cwl.parsing.tool import CLTParser
from janis_core.ingestion.cwl.parsing.tool import CLTRequirementsParser
from janis_core.ingestion.cwl.types import ingest_cwl_type
from janis_core.ingestion.cwl import parse as parse_cwl

from janis_core.types import File
from janis_core.types import GenericFileWithSecondaries
from janis_core.messages import configure_logging
from janis_core.messages import load_loglines
from janis_core.messages import ErrorCategory
from janis_core import settings


def _do_setup() -> None:
    configure_logging()
    settings.ingest.SAFE_MODE = False


CWL_TESTDATA_DIR = os.path.join(os.getcwd(), 'janis_core/tests/data/cwl')

def _load_cwl_doc(filepath: str) -> Tuple[Any, Any]:
    cwl_version = load_cwl_version(filepath)
    cwl_utils = load_cwl_utils_from_version(cwl_version)
    loaded_doc = load_cwl_document(filepath)
    return loaded_doc, cwl_utils


class TestBasicFunctionality(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()

    def test_tool_bwa_index(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/BWA-Index.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils=cwl_utils, clt=clt, entity=clt)
        tool = parser.parse()
        msgs = load_loglines()
        self.assertEqual(tool.id(), 'BWA_Index')
        self.assertEqual(len(tool._inputs), 3)
        self.assertEqual(len(tool._arguments), 2)

    def test_tool_bwa_mem(self):
        pass



class TestRequirementsParsing(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()

    def test_directories_to_create1(self):
        pass

    def test_directories_to_create2(self):
        pass
        
    def test_directories_to_create3(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/prepare_fasta_db.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_directories_to_create3')
        reqs = parser.parse()
        msgs = load_loglines(tool_uuid='test_directories_to_create3')
        print()

    def test_files_to_create0(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/bundle_secondaryfiles.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create0')
        reqs = parser.parse()
        msgs = load_loglines(tool_uuid='test_files_to_create0')
        self.assertIn('inputs.primary_file', reqs['files_to_create'])
        self.assertIsInstance(reqs['files_to_create']['inputs.primary_file'], InputSelector)
        self.assertIn('inputs.secondary_files', reqs['files_to_create'])
        self.assertIsInstance(reqs['files_to_create']['inputs.secondary_files'], InputSelector)
    
    def test_files_to_create1(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/bismark_prepare_genome.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create1')
        reqs = parser.parse()
        msgs = load_loglines(tool_uuid='test_files_to_create1')
    
    def test_files_to_create2(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/BWA-Mem2-index.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create2')
        reqs = parser.parse()
        msgs = load_loglines(tool_uuid='test_files_to_create2')
    
    def test_files_to_create3(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/bzip2_compress.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create3')
        reqs = parser.parse()
        msgs = load_loglines(tool_uuid='test_files_to_create3')
    
    def test_files_to_create4(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/component_segmentation.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create4')
        reqs = parser.parse()
        msgs = load_loglines(tool_uuid='test_files_to_create4')
    
    def test_files_to_create5(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/cellranger-aggr.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create5')
        reqs = parser.parse()
        msgs = load_loglines(tool_uuid='test_files_to_create5')
 
    # def test_files_to_create9(self):
    #     filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/bundle_secondaryfiles.cwl'
    #     clt, cwl_utils = _load_cwl_doc(filepath)
    #     parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt)
    #     reqs = parser.do_parse()
    #     self.assertDictEqual(reqs['files_to_create'], {})

    # def test_files_to_create10(self):
    #     # expression tool
    #     filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/return_directory.cwl'
    #     etool, cwl_utils = _load_cwl_doc(filepath)
        
    #     parser = CWlParser(filepath)
    #     cltool = parser.parse_etool_to_cltool(etool)
        
    #     parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=cltool, entity=cltool, is_expression_tool=True)
    #     reqs = parser.do_parse()
    #     self.assertIn('return_directory.js', reqs['files_to_create'])
    #     self.assertIsInstance(reqs['files_to_create']['return_directory.js'], str)
 
    # def test_files_to_create11(self):
    #     filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/check_threshold.cwl'
    #     clt, cwl_utils = _load_cwl_doc(filepath)
    #     parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt)
    #     reqs = parser.do_parse()
    #     self.assertIn('check_threshold.py', reqs['files_to_create'])
    #     self.assertIsInstance(reqs['files_to_create']['check_threshold.py'], str)
    #     self.assertIn('#!/usr/bin/env python3', reqs['files_to_create']['check_threshold.py'])
 
    # def test_files_to_create13(self):
    #     filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/amplicon_architect.cwl'
    #     clt, cwl_utils = _load_cwl_doc(filepath)
    #     parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create13')
    #     reqs = parser.parse()
    #     msgs = get_messages('test_files_to_create13')
    #     self.assertIn('likely untranslated cwl / js in script: setup_vars.sh', msgs)
    
    # def test_files_to_create14(self):
    #     filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/control-freec-11-6-sbg.cwl'
    #     clt, cwl_utils = _load_cwl_doc(filepath)
    #     parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create14')
    #     reqs = parser.parse()
    #     self.assertIn('config.txt', reqs['files_to_create'])
    #     msgs = get_messages('test_files_to_create14')
    #     self.assertIn('config.txt: js code to dynamically create runtime file. please address', msgs)
    #     print()
    
    # def test_files_to_create6(self):
    #     filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/antismash_v4.cwl'
    #     clt, cwl_utils = _load_cwl_doc(filepath)
    #     parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create6')
    #     reqs = parser.parse()
    #     msgs = get_messages('test_files_to_create6')
    #     print()
    
    # def test_files_to_create7(self):
    #     filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/add_crossmapped2resmapped.cwl'
    #     clt, cwl_utils = _load_cwl_doc(filepath)
    #     parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test_files_to_create7')
    #     reqs = parser.parse()
    #     msgs = get_messages('test_files_to_create7')
    #     print()
        
    # def test_files_to_create8(self):
    #     # expression tool
    #     filepath = f'{CWL_TESTDATA_DIR}/tools/expressiontools/check_value.cwl'
    #     etool, cwl_utils = _load_cwl_doc(filepath)
        
    #     parser = CWlParser(filepath)
    #     cltool = parser.parse_etool_to_cltool(etool)
        
    #     parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=cltool, entity=cltool, is_expression_tool=True)
    #     reqs = parser.do_parse()

    #     # check js script file created, and is a string
    #     self.assertIn('check_value.js', reqs['files_to_create'])
    #     self.assertIsInstance(reqs['files_to_create']['check_value.js'], str)

    #     # check js script file modified correctly
    #     script = reqs['files_to_create']['check_value.js']
    #     self.assertNotIn('"use strict";\nvar inputs=$(inputs);\nvar runtime=$(runtime);', script)
    #     self.assertIn('var inputs = JSON.parse( process.argv[2] );\n', script)

    def test_env_vars1(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/requirements/amplicon_architect.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTRequirementsParser(cwl_utils=cwl_utils, clt=clt, entity=clt, tool_uuid='test')
        reqs = parser.do_parse()
        self.assertIn('AA_SRC', reqs['env_vars'])
        self.assertEqual(reqs['env_vars']['AA_SRC'], '/home/programs/AmpliconArchitect-master/src')


class TestErrorHandlingFallbacks(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()

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



class TestDatatypes(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()
    
    def test_primitive_types(self) -> None:
        filepath = f'{CWL_TESTDATA_DIR}/tools/datatypes/basic_types.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        self.assertIsInstance(jtool, CommandToolBuilder)
        self.assertEqual(len(jtool.inputs_map()), 6)
        actual_types = {x.id(): x.intype.__repr__() for x in jtool.inputs_map().values()}
        expected_types = {
            'inFile': 'File',
            'inInt': 'Integer',
            'inStr': 'String',
            'inFileOpt': 'Optional<File>',
            'inIntOpt': 'Optional<Integer>',
            'inStrOpt': 'Optional<String>',
        }
        self.assertDictEqual(actual_types, expected_types)
    
    def test_edam_types(self) -> None:
        filepath = f'{CWL_TESTDATA_DIR}/tools/BWA-Index.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils=cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        self.assertIsInstance(jtool, CommandToolBuilder)
        
        # inputs
        actual_inp_types = {x.id(): x.input_type.__repr__() for x in jtool._inputs}
        expected_inp_types = {
            'algo_type': 'Optional<String>',
            'index_name': 'Optional<String>',
            'sequences': 'Fasta',
        }
        self.assertDictEqual(actual_inp_types, expected_inp_types)
        
        # outputs
        actual_out_types = {x.id(): x.output_type.__repr__() for x in jtool._outputs}
        expected_out_types = {
            'index': 'GenericFileWithSecondaries [^.amb, ^.ann, ^.pac, ^.sa]',
        }
        self.assertDictEqual(actual_out_types, expected_out_types)

    
    def test_array_types(self) -> None:
        filepath = f'{CWL_TESTDATA_DIR}/tools/datatypes/array_types.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        self.assertIsInstance(jtool, CommandToolBuilder)
        self.assertEqual(len(jtool.inputs_map()), 6)
        actual_types = {x.id(): x.intype.__repr__() for x in jtool.inputs_map().values()}
        expected_types = {
            'inFileArr': 'Array<File>',
            'inIntArr': 'Array<Integer>',
            'inStrArr': 'Array<String>',
            'inFileArrOpt': 'Optional<Array<File>>',
            'inIntArrOpt': 'Optional<Array<Integer>>',
            'inStrArrOpt': 'Optional<Array<String>>',
        }
        self.assertDictEqual(actual_types, expected_types)
    
    def test_secondary_types(self) -> None:
        filepath = f'{CWL_TESTDATA_DIR}/tools/datatypes/secondary_types.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        self.assertIsInstance(jtool, CommandToolBuilder)
        self.assertEqual(len(jtool.inputs_map()), 8)
        actual_types = {x.id(): x.intype.__repr__() for x in jtool.inputs_map().values()}
        expected_types = {
            'inSecondary1': 'IndexedBam',
            'inSecondary2': 'GenericFileWithSecondaries []',
            'inSecondaryArr1': 'Array<IndexedBam>',
            'inSecondaryArr2': 'Array<GenericFileWithSecondaries []>',
            'inSecondaryOpt1': 'Optional<IndexedBam>',
            'inSecondaryOpt2': 'Optional<GenericFileWithSecondaries> []',
            'inSecondaryArrOpt1': 'Optional<Array<IndexedBam>>',
            'inSecondaryArrOpt2': 'Optional<Array<GenericFileWithSecondaries []>>',
        }
        self.assertDictEqual(actual_types, expected_types)



class TestErrorHandlingDatatypes(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()

    def test_multiple_datatypes(self):
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        from cwl_utils.parser import cwl_v1_0 as cwl_utils
        cwl_type = ['null', 'File', 'string']
        dtype = ingest_cwl_type(cwl_type, cwl_utils, None, tool_uuid='test', secondaries=None)
        self.assertIsInstance(dtype, File)
        lines = load_loglines(tool_uuid='test')
        msgs = [x.message for x in lines]
        self.assertIn("entity supports multiple datatypes: ['File', 'String']. selected File as fallback.", msgs)
    
    def test_unparseable_generic(self):
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        from cwl_utils.parser import cwl_v1_0 as cwl_utils
        cwl_type = 'kitten'
        dtype = ingest_cwl_type(cwl_type, cwl_utils, None, tool_uuid='test', secondaries=None)
        self.assertIsInstance(dtype, File)
        lines = load_loglines(tool_uuid='test')
        msgs = [x.message for x in lines]
        self.assertIn('unsupported datatype: kitten. treated as generic File.', msgs)
    
    def test_enum_type(self):
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        from cwl_utils.parser import cwl_v1_0 as cwl_utils
        cwl_type = f'file://{CWL_TESTDATA_DIR}/tools/gatk_haplotype_tool.cwl#annotation_type'
        dtype = ingest_cwl_type(cwl_type, cwl_utils, None, tool_uuid='test', secondaries=None)
        self.assertIsInstance(dtype, File)
        lines = load_loglines(tool_uuid='test')
        msgs = [x.message for x in lines]
        self.assertIn('unsupported datatype: file:///home/grace/work/pp/translation/janis-core/janis_core/tests/data/cwl/tools/gatk_haplotype_tool.cwl#annotation_type. treated as generic File.', msgs)

    def test_unparseable_secondary_type(self):
        settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        cinp = clt.inputs[1]
        dtype = ingest_cwl_type(cinp.type, cwl_utils, None, tool_uuid='test', secondaries=cinp.secondaryFiles)
        self.assertIsInstance(dtype, GenericFileWithSecondaries)
        lines = load_loglines(tool_uuid='test')
        msgs = [x.message for x in lines]
        self.assertIn("could not parse secondaries format from javascript expression. treated as generic File with secondaries.", msgs)




class TestErrorHandlingExpressions(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()

    def test_resource_requirements(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/resources.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        
        parser = CLTRequirementsParser(cwl_utils, clt=clt, entity=clt, tool_uuid='test')
        requirements = parser.parse()

        lines = load_loglines(tool_uuid='test', category=ErrorCategory.SCRIPT)
        msgs = [x.message for x in lines]
        
        self.assertEqual(requirements["memory"], '__TOKEN1__')
        self.assertEqual(requirements["cpus"], '__TOKEN2__')
        self.assertEqual(requirements["time"], '__TOKEN3__')

        self.assertIn('__TOKEN1__ = "$(Math.round((953.674 * [inputs.runtime_memory, ((inputs.inputFile.size / 1048576) > 1024) ? 4 : 2, 4].filter(function (inner) { return inner != null })[0])))"', msgs)
        self.assertIn('__TOKEN2__ = "$([inputs.runtime_cpu, (2 * inputs.outputFiles), 1].filter(function (inner) { return inner != null })[0])"', msgs)
        self.assertIn('__TOKEN3__ = "$([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])"', msgs)

    def test_initial_workdir_requirement_listing(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)

        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        lines = load_loglines(jtool.uuid)
        msgs = [x.message for x in lines]
        self.assertIsNone(jtool.files_to_create())
        
        expected_msg = 'directory required for runtime generated using cwl / js: <js>${    return [{"class": "Directory",            "basename": "subdir",            "listing": [ inputs.example ]            }]}</js>. please address'
        self.assertIn(expected_msg, msgs)
        
    def test_initial_workdir_requirement_dirent_entryname(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)

        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        lines = load_loglines(jtool.uuid)
        msgs = [x.message for x in lines]
        self.assertIsNone(jtool.files_to_create())

        expected_msg = 'untranslated javascript expression in environment variable value'
        self.assertIn(expected_msg, msgs)
    
    def test_initial_workdir_requirement_dirent_entry(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)

        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        lines = load_loglines(jtool.uuid)
        msgs = [x.message for x in lines]
        self.assertIsNone(jtool.files_to_create())

        expected_msg = 'untranslated javascript expression in environment variable value'
        self.assertIn(expected_msg, msgs)

    def test_envvar_requirement_envdef_envvalue(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/requirements.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTParser(cwl_utils, clt=clt, entity=clt)
        jtool = parser.parse()
        lines = load_loglines(jtool.uuid)
        msgs = [x.message for x in lines]
        self.assertEqual(jtool.env_vars()['test1'], '<js>9 + 10</js>')

    def test_clt_argument_valuefrom(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        
        parser = CLTArgumentParser(cwl_utils, clt=clt, entity=clt.arguments[5], tool_uuid='test')
        arg = parser.parse()
        lines = load_loglines(tool_uuid='test')
        msgs = [x.message for x in lines]
        self.assertEqual(arg.value, '__TOKEN1__')
        self.assertIn('__TOKEN1__ = "${  var r = [];  for (var i = 10; i >= 1; i--) {    r.push(i);  }  return r;}"', msgs)
    
    @unittest.skip("TODO leave log_info message")
    def test_clt_input_format(self):
        # don't need to worry about File format in cwl - type is ingested as File
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTInputParser(cwl_utils, clt=clt, entity=clt.inputs[3], tool_uuid='test')
        tinput = parser.parse()
    
    def test_clt_input_secondaryfiles(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTInputParser(cwl_utils, clt=clt, entity=clt.inputs[1], tool_uuid='test')
        tinput = parser.parse()
        self.assertIsInstance(tinput.input_type, GenericFileWithSecondaries)
        error_msgs = load_loglines(tinput.uuid)
        expected_msg = "could not parse secondaries format from javascript expressions: $(self.basename + self.nameext.replace('m','i'))"
        self.assertIn(expected_msg, error_msgs)

    def test_clt_input_inputbinding_valuefrom(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/inputs_arguments.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        
        parser = CLTInputParser(cwl_utils, clt=clt, entity=clt.inputs[5], tool_uuid='test')
        tinput = parser.parse()
        expected_value = '<js>[inputs.runtime_cpu, 16, 1].filter(function (inner) { return inner != null })[0]</js>'
        self.assertEqual(tinput.value, expected_value)
    
    @unittest.skip("TODO leave log_info message")
    def test_clt_output_format(self):
        # don't need to worry about File format in cwl - type is ingested as File
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[1], tool_uuid='test')
        tout = parser.parse()

    def test_clt_output_secondaryfiles(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[1], tool_uuid='test')
        tout = parser.parse()
        self.assertIsInstance(tout.output_type, GenericFileWithSecondaries)
        error_msgs = load_loglines(tout.uuid)
        expected_msg = "could not parse secondaries format from javascript expressions: $(inputs.bam.basename + inputs.bam.nameext.replace('m','i'))"
        self.assertIn(expected_msg, error_msgs)
    
    def test_clt_output_outputbinding_glob(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[0], tool_uuid='test')
        tout = parser.parse()
        self.assertIsInstance(tout.selector, WildcardSelector)
        self.assertEqual(tout.selector.wildcard, 'output.txt')

        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[2], tool_uuid='test')
        tout = parser.parse()
        self.assertIsInstance(tout.selector, StringFormatter)
        self.assertEqual(tout.selector._format, '{JANIS_CWL_TOKEN_1}.trimmed.fastq')
        self.assertEqual(tout.selector.kwargs['JANIS_CWL_TOKEN_1'], "<js>inputs.fastq.nameroot.replace(/\\b.fastq\\b/g, '')</js>")
    
    def test_clt_output_outputbinding_outputEval(self):
        filepath = f'{CWL_TESTDATA_DIR}/tools/expressions/outputs.cwl'
        clt, cwl_utils = _load_cwl_doc(filepath)
        parser = CLTOutputParser(cwl_utils, clt=clt, entity=clt.outputs[3], tool_uuid='test')
        tout = parser.parse()
        self.assertIsInstance(tout.selector, ReadContents)
        self.assertEqual(tout.selector.args[0], '<js>self[0]</js>')
    
    @unittest.skip("TODO leave log_info message")
    def test_wf_input_format(self) -> None:
        # don't need to worry about File format in cwl - type is ingested as File
        raise NotImplementedError

    def test_wf_input_secondaryfiles(self):
        filepath = f'{CWL_TESTDATA_DIR}/workflows/expressions.cwl'
        wf = parse_cwl(filepath)
        winp = wf.input_nodes['bambai_pair2']
        self.assertIsInstance(winp.datatype, GenericFileWithSecondaries)
        self.assertEqual(len(winp.datatype.secondaries), 0)
        error_msgs = load_loglines(winp.uuid)
        expected_msg = "could not parse secondaries format from javascript expressions: $(self.basename + self.nameext.replace('m','i'))"
        self.assertIn(expected_msg, error_msgs)
    
    def test_step_input_valuefrom(self):
        filepath = f'{CWL_TESTDATA_DIR}/workflows/expressions.cwl'
        wf = parse_cwl(filepath)
        
        # example 1
        jstep = wf.step_nodes['echo2']
        inp_name = 'text'
        
        self.assertIn(inp_name, jstep.tool.connections)
        self.assertEqual(jstep.tool.connections[inp_name], "<js>$(inputs.prefix_str + '_' + inputs.suffix_str)</js>")
        self.assertIn(inp_name, jstep.sources)
        self.assertEqual(jstep.sources[inp_name].source_map[0].source.input_node.default, "<js>$(inputs.prefix_str + '_' + inputs.suffix_str)</js>")
        
        msgs = load_loglines(jstep.uuid)
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
        
        msgs = load_loglines(jstep.uuid)
        self.assertIn('\'target_filename\' input value contains untranslated javascript expression: <js>$(inputs.bambai_pair1.location.split(\'/\').slice(-1)[0].split(\'.\').slice(0,-1).join(\'.\')+"_default_s_enhcr.png")</js>', msgs)

    @unittest.skip("TODO leave log_info message")
    def test_wf_output_format(self) -> None:
        # TODO implement this - log info message about the format which was discarded
        # don't need to worry about File format in cwl - type is ingested as File
        raise NotImplementedError

    def test_wf_output_secondaryfiles(self):
        filepath = f'{CWL_TESTDATA_DIR}/workflows/expressions.cwl'
        wf = parse_cwl(filepath)
        winp = wf.output_nodes['out3']
        self.assertIsInstance(winp.datatype, GenericFileWithSecondaries)
        self.assertEqual(len(winp.datatype.secondaries), 0)
        error_msgs = load_loglines(winp.uuid)
        expected_msg = "could not parse secondaries format from javascript expressions: $(self.basename + self.nameext.replace('m','i'))"
        self.assertIn(expected_msg, error_msgs)
  





class TestIdentifiers(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()

    def test_filename_id(self):
        ident = f'file://{CWL_TESTDATA_DIR}/workflows/super_enhancer_wf.cwl'
        ref = get_cwl_reference(ident)
        self.assertEqual(ref.filename, 'super_enhancer_wf')
        self.assertEqual(ref.internal_path, None)
        self.assertEqual(ref.entity, None)
    
    def test_entity_id(self):
        ident = f'file://{CWL_TESTDATA_DIR}/workflows/super_enhancer_wf.cwl#annotation_file'
        ref = get_cwl_reference(ident)
        self.assertEqual(ref.filename, 'super_enhancer_wf')
        self.assertEqual(ref.internal_path, None)
        self.assertEqual(ref.entity, 'annotation_file')
    
    def test_internal_path_id(self):
        ident = f'file://{CWL_TESTDATA_DIR}/workflows/super_enhancer_wf.cwl#assign_genes/result_file'
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

    def setUp(self) -> None:
        _do_setup()

    def test_convert_cwl_types_to_python_clt(self):
        doc = f'{CWL_TESTDATA_DIR}/tools/gatk_haplotype_caller.cwl'
        clt = load_cwl_document(doc)
        # no idea how to test :/
    
    def test_convert_cwl_types_to_python_workflow(self):
        doc = f'{CWL_TESTDATA_DIR}/workflows/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
        # no idea how to test :/
    
    def test_inline_cltool_ids(self):
        doc = f'{CWL_TESTDATA_DIR}/workflows/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
        for step in workflow.steps:
            tool_id = step.run.id
            self.assertTrue(tool_id.endswith('_tool.cwl'))
    



class TestDocumentLoading(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()

    def test_load_cwl_document(self):
        doc = f'{CWL_TESTDATA_DIR}/workflows/super_enhancer_wf.cwl'
        workflow = load_cwl_document(doc)
        self.assertIsNotNone(workflow)
    
    def test_load_cwl_version(self):
        doc = f'{CWL_TESTDATA_DIR}/workflows/structuralvariants/workflow.cwl'
        version = load_cwl_version(doc)
        self.assertIsNotNone(version)
        


class TestParseExpression(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup()

    # objects
    def test_input(self):
        expr = "$(inputs.sequences)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, InputSelector)
        self.assertEqual("sequences", result.input_to_select)
        
        expr = "$(inputs.reads_1)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, InputSelector)
        self.assertEqual("reads_1", result.input_to_select)
    
    def test_self(self):
        # should return None, False. 
        # can't do much about CWL obj.self.
        expr = "$(self.location.split('/'))"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, str)
        self.assertEqual(result, '__TOKEN1__')
    
    # runtime 
    def test_runtime_outdir(self):
        expr = "$(runtime.outdir)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, str)
        self.assertEqual(result, '.')
    
    def test_runtime_tmpdir(self):
        expr = "$(runtime.tmpdir)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, str)
        self.assertEqual(result, '.')
    
    def test_runtime_outdir_size(self):
        expr = "$(runtime.outdirSize)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, int)
        self.assertEqual(result, 1024)
    
    def test_runtime_tmpdir_size(self):
        expr = "$(runtime.tmpdirSize)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, int)
        self.assertEqual(result, 1024)
    
    def test_runtime_cores(self):
        expr = "$(runtime.cores)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, CpuSelector)
    
    def test_runtime_ram(self):
        expr = "$(runtime.ram)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, MemorySelector)
    
    # logical - misc
    def test_if(self):
        expr = '$(inputs.index_name !== null ? inputs.index_name : inputs.sequences.nameroot)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, If)
        self.assertIsInstance(result.args[0], InequalityOperator)
        self.assertIsInstance(result.args[0].args[0], InputSelector)
        self.assertIsNone(result.args[0].args[1])
        self.assertIsInstance(result.args[1], InputSelector)
        self.assertIsInstance(result.args[2], NamerootOperator)
        self.assertIsInstance(result.args[2].args[0], InputSelector)
    
    def test_group(self):
        expr = '$((1 + 2))'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, GroupOperator)
        self.assertIsInstance(result.args[0], AddOperator)

    # logical - two value
    def test_and(self):
        expr = '$(inputs.filelist && inputs.filelist_mates)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, AndOperator)
    
    def test_or(self):
        expr = '$(inputs.filelist || inputs.filelist_mates)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, OrOperator)
    
    def test_eq(self):
        expr = '$(inputs.depth == "-bg")'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, EqualityOperator)
        
        expr = '$(inputs.depth === "-bg")'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, EqualityOperator)
    
    def test_ineq(self):
        expr = '$(inputs.depth != "-bg")'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, InequalityOperator)
        
        expr = '$(inputs.depth !== "-bg")'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, InequalityOperator)
    
    def test_gte(self):
        expr = '$(inputs.depth >= 20)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, GteOperator)
    
    def test_gt(self):
        expr = '$(inputs.depth > 20)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, GtOperator)
    
    def test_lte(self):
        expr = '$(inputs.depth <= 20)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, LteOperator)
    
    def test_lt(self):
        expr = '$(inputs.depth < 20)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, LtOperator)
    
    def test_add(self):
        expr = '$(1 + 2)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, AddOperator)
    
    def test_sub(self):
        expr = '$(1 - 2)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, SubtractOperator)
    
    def test_mul(self):
        expr = '$(1 * 2)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, MultiplyOperator)
    
    def test_div(self):
        expr = '$(1 / 2)'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, DivideOperator)
    
    # logical - math
    def test_floor(self):
        expr = '$(Math.floor(2.3))'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, FloorOperator)
    
    def test_ceil(self):
        expr = '$(Math.ceil(2.3))'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, CeilOperator)
    
    def test_round(self):
        expr = '$(Math.round(2.3))'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, RoundOperator)

    # attributes - file
    def test_file_basename(self):
        expr = "$(inputs.reads_1.basename)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, BasenameOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        
    def test_file_nameroot(self):
        expr = "$(inputs.reads_1.nameroot)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, NamerootOperator)
        self.assertIsInstance(result.args[0], InputSelector)
    
    def test_file_nameext(self):
        expr = "$(inputs.reads_1.nameext)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, NameextOperator)
        self.assertIsInstance(result.args[0], InputSelector)
    
    def test_file_dirname(self):
        expr = "$(inputs.reads_1.dirname)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, DirnameOperator)
        self.assertIsInstance(result.args[0], InputSelector)
    
    def test_file_size(self):
        expr = "$(inputs.reads_1.size)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, FileSizeOperator)
        self.assertIsInstance(result.args[0], InputSelector)
    
    def test_file_contents(self):
        expr = "$(inputs.reads_1.contents)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, ReadContents)
        self.assertIsInstance(result.args[0], InputSelector)
    
    # methods - array 
    def test_arr_join(self):
        expr = "$(inputs.reads.join('.'))"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, JoinOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertIsInstance(result.args[1], str)
        self.assertEqual(result.args[1], "'.'")
    
    def test_arr_slice(self):
        expr = "$(inputs.reads.slice(0,-1))"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, SliceOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertIsInstance(result.args[1], int)
        self.assertEqual(result.args[1], 0)
        self.assertIsInstance(result.args[2], int)
        self.assertEqual(result.args[2], -1)
        
        expr = "$(inputs.reads.slice(3))"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, SliceOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertIsInstance(result.args[1], int)
        self.assertEqual(result.args[1], 3)
        self.assertIsNone(result.args[2])
    
    def test_arr_length(self):
        expr = "$(inputs.readgroups.length)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, LengthOperator)
        self.assertIsInstance(result.args[0], InputSelector)
    
    def test_arr_flatten(self):
        expr = "$(inputs.readgroups.flat())"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, FlattenOperator)
        self.assertIsInstance(result.args[0], InputSelector)
    
    def test_arr_index(self):
        expr = "$(inputs.readgroups[1])"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, IndexOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertIsInstance(result.args[1], int)
        self.assertEqual(result.args[1], 1)

    # methods - string
    def test_str_split(self):
        expr = "$(inputs.sequences.split('/'))"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, SplitOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertIsInstance(result.args[1], str)
        self.assertEqual(result.args[1], "'/'")

    def test_str_replace(self):
        expr = "$(inputs.name.replace(/ /g,''))"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, ReplaceOperator)
        self.assertIsInstance(result.args[0], InputSelector)
        self.assertIsInstance(result.args[1], str)
        self.assertEqual(result.args[1], "/ /g")
        self.assertIsInstance(result.args[2], str)
        self.assertEqual(result.args[2], "''")
    
    # complex cases
    def test_complex_stringformatter(self):
        expr = "$(inputs.reads_1.basename).trimmed$(inputs.reads_1.nameext)"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, StringFormatter)
        self.assertEqual(result._format, '{token1}.trimmed{token2}')
        self.assertIsInstance(result.kwargs['token1'], BasenameOperator)
        self.assertIsInstance(result.kwargs['token2'], NameextOperator)

    def test_complex_methodchain(self):
        expr = "$(inputs.index.split('/').slice(-1)[0].split('.').slice(-1)[0])"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, IndexOperator)
        self.assertIsInstance(result.args[0], SliceOperator)
        self.assertEqual(result.args[1], 0)
        result = result.args[0]
        self.assertIsInstance(result.args[0], SplitOperator)
        self.assertEqual(result.args[1], -1)
        self.assertEqual(result.args[2], None)
        result = result.args[0]
        self.assertIsInstance(result.args[0], IndexOperator)
        self.assertEqual(result.args[1], "'.'")
        result = result.args[0]
        self.assertEqual(result.args[1], 0)
        self.assertIsInstance(result.args[0], SliceOperator)
        result = result.args[0]
        self.assertEqual(result.args[1], -1)
        self.assertEqual(result.args[2], None)
        self.assertIsInstance(result.args[0], SplitOperator)
        result = result.args[0]
        self.assertEqual(result.args[1], "'/'")
        self.assertIsInstance(result.args[0], InputSelector)
        result = result.args[0]
    
    def test_complex_concat1(self):
        expr = "$('> ' + 'stats.log')"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, AddOperator)
        self.assertEqual(result.args[0], "'> '")
        self.assertEqual(result.args[1], "'stats.log'")
    
    def test_complex_concat2(self):
        expr = "$('> ' + inputs.index_base_name + '.log')"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, AddOperator)
        self.assertIsInstance(result.args[0], AddOperator)
        self.assertEqual(result.args[1], "'.log'")
        result = result.args[0]
        self.assertIsInstance(result, AddOperator)
        self.assertEqual(result.args[0], "'> '")
        self.assertIsInstance(result.args[1], InputSelector)
    
    def test_complex_math1(self):
        expr = "$(Math.ceil(inputs.target.size/(1024*1024*1024) + 20))"
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, CeilOperator)
        self.assertIsInstance(result.args[0], AddOperator)
        result = result.args[0]
        self.assertIsInstance(result.args[0], DivideOperator)
        self.assertEqual(result.args[1], 20)
        result = result.args[0]
        self.assertIsInstance(result.args[0], FileSizeOperator)
        self.assertIsInstance(result.args[1], GroupOperator)
        self.assertIsInstance(result.args[0].args[0], InputSelector)
        result = result.args[1].args[0]
        self.assertIsInstance(result.args[0], MultiplyOperator)
        self.assertEqual(result.args[1], 1024)
        result = result.args[0]
        self.assertEqual(result.args[0], 1024)
        self.assertEqual(result.args[1], 1024)

    def test_complex_stringmethods(self):
        expr = '$("/foo/bar/baz".split(\'/\').slice(-1)[0])'
        result, success = parse_expression(expr, 'blank')
        print()

    # javascript blocks
    def test_block_number(self):
        expr = "${return 16 }"
        result, success = parse_expression(expr, 'blank')
        self.assertEqual(16, result)

    def test_block_number_with_semicolon(self):
        expr = "${return 16;}"
        result, success = parse_expression(expr, 'blank')
        self.assertEqual(16, result)

    def test_block_number_with_spaces(self):
        expr = "${ return 80000 }"
        result, success = parse_expression(expr, 'blank')
        self.assertEqual(80000, result)

    def test_block_string(self):
        expr = '${ return "over 80000" }'
        result, success = parse_expression(expr, 'blank')
        self.assertEqual('"over 80000"', result)

    def test_block_if(self):
        expr = """${
  if( inputs.output_name == null ){
    return inputs.bedgraph.basename;
  }
  else{
    return inputs.output_name;
  }
}"""
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, If)
        self.assertIsInstance(result.args[0], GroupOperator)
        cond_check = result.args[0].args[0]
        cond_true = result.args[1]
        cond_false = result.args[2]
        self.assertIsInstance(cond_check, EqualityOperator)
        self.assertIsInstance(cond_check.args[0], InputSelector)
        self.assertIsNone(cond_check.args[1])
        self.assertIsInstance(cond_true, BasenameOperator)
        self.assertIsInstance(cond_false, InputSelector)
    
    def test_block_complex(self):
        expr = '${ return(parseInt(runtime.ram/runtime.cores-100).toString() + "M") }'
        result, success = parse_expression(expr, 'blank')
        self.assertIsInstance(result, AddOperator)
        expected = '(str(int(((inputs.runtime_memory / inputs.runtime_cpu) - 100))) + "M")'
        actual = str(result)
        self.assertEqual(actual, expected)
