
from typing import Any, Optional
import unittest
import os 
import json
import xml.etree.ElementTree as et

from janis_core.ingestion.main import ingest_galaxy
from janis_core.ingestion.galaxy.gx.gxtool.text.simplification.aliases import resolve_aliases

from janis_core.ingestion.galaxy import settings
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.tool_state import load_tool_state
from janis_core.ingestion.galaxy.gx.gxtool.text.cheetah.evaluation import sectional_evaluate
from janis_core.ingestion.galaxy.gx.gxtool.load import load_xmltool

from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.tool_step.metadata import parse_step_metadata
from janis_core.ingestion.galaxy import containers
from janis_core.ingestion.galaxy.gx.gxtool.requirements import CondaRequirement

from janis_core.ingestion.galaxy import datatypes
from janis_core.ingestion.galaxy.datatypes.core import file_t, string_t, bool_t

from janis_core import WorkflowBuilder
from janis_core import WorkflowMetadata
from janis_core.workflow.workflow import InputNode
from janis_core.workflow.workflow import OutputNode
from janis_core.workflow.workflow import StepNode
from janis_core.workflow.workflow import ScatterMethod

from janis_core import CommandToolBuilder
from janis_core import CommandTool
from janis_core import ToolInput
from janis_core import ToolOutput
from janis_core import ToolMetadata
from janis_core import (
    InputSelector,
    WildcardSelector
)
from janis_core import (
    Stdout,
    Array,
    Float,
    Boolean,
    File,
    String
)

from janis_core.ingestion.galaxy import datatypes

from janis_core.ingestion.galaxy.janis_mapping.workflow import to_janis_workflow
from janis_core.ingestion.galaxy.janis_mapping.workflow import to_janis_inputs_dict
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_datatype
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_selector
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_metadata
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_tool_input
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_tool_output
from janis_core.ingestion.galaxy.janis_mapping.tool import to_janis_tool

# mock objects
from .mock.mock_components import MOCK_POSITIONAL1
from .mock.mock_components import MOCK_FLAG1
from .mock.mock_components import MOCK_OPTION2
from .mock.mock_components import MOCK_REDIRECT_OUTPUT
from .mock.mock_entities import MOCK_WORKFLOW_INPUT1
from .mock.mock_tool import MOCK_TOOL_ABRICATE
from .mock.mock_workflow import MOCK_WORKFLOW

UNICYCLER_VANILLA_PATH = './janis_core/tests/data/command/manipulation/template/unicycler/unicycler_command.xml'
UNICYCLER_TEMPLATED_PATH = './janis_core/tests/data/command/manipulation/template/unicycler/unicycler_command_templated.xml'
UNICYCLER_INPUTS_PATH = './janis_core/tests/data/command/manipulation/template/unicycler/unicycler_step.json'
 
QUERY1 = CondaRequirement(_name='abricate', _version='1.0.1')
QUERY1_EXPECTED_RESULT = 'quay.io/biocontainers/abricate:1.0.1--ha8f3691_1'

QUERY2 = CondaRequirement(_name='samtools', _version='1.15')
QUERY2_EXPECTED_RESULT = 'quay.io/biocontainers/samtools:1.15--h1170115_1'

QUERY3 = CondaRequirement(_name='cutadapt', _version='3.5')
QUERY3_EXPECTED_RESULT = 'quay.io/biocontainers/cutadapt:3.5--py36h91eb985_1'

from janis_core.ingestion import ingest
from janis_core.translations import translate


def run(filepath: str, srcfmt: str, destfmt: str) -> Optional[str]:
    wf = ingest(filepath, srcfmt)
    return translate(wf, destfmt, allow_empty_container=True, export_path='./translated')


def _load_galaxy_workflow(filepath: str) -> dict[str, Any]:
    with open(filepath, 'r') as fp:
        return json.load(fp)
    
# def _load_tool_command_and_state(step: dict[str, Any], flat: bool=False) -> Tuple[str, dict[str, Any]]:
#     metadata = parse_step_metadata(step)
#     settings.tool.update_via_wrapper(metadata.wrapper)
#     xmltool = load_xmltool(settings.tool.tool_path)
#     tool_state = load_tool_state(step, flat=flat)
#     cmdstr = simplify_xmltool_command(
#         xmltool=xmltool, 
#         inputs_dict=tool_state,
#         # additional_filters=['remove_dataset_attributes', 'remove_dataset_methods']
#     )
#     return cmdstr, tool_state

def _load_tool_state(step: dict[str, Any], additional_filters: list[str]=[]) -> dict[str, Any]:
    metadata = parse_step_metadata(step)
    settings.tool.update_via_wrapper(metadata.wrapper)
    return load_tool_state(step, additional_filters=additional_filters)


class TestLoadToolState(unittest.TestCase):
    """
    tests ability to map or generate janis datatypes / selectors etc
    from internal tool objects.
    """
    def setUp(self) -> None:
        datatypes.populate()
        self.tool = MOCK_TOOL_ABRICATE

    def test_default_filters(self) -> None:
        wf_path = os.path.abspath('./janis_core/tests/data/galaxy/cutadapt_wf.ga')
        gx_workflow = _load_galaxy_workflow(wf_path)
        gx_step = gx_workflow['steps']['2']
        tool_state = _load_tool_state(gx_step)
        self.assertEqual(tool_state['output_selector'], None)
        self.assertEqual(tool_state['adapter_options']['action'], 'trim')
        self.assertEqual(tool_state['library']['r1']['cut'], '0')
    
    def test_null_varname_filter(self) -> None:
        wf_path = os.path.abspath('./janis_core/tests/data/galaxy/cutadapt_wf.ga')
        gx_workflow = _load_galaxy_workflow(wf_path)
        gx_step = gx_workflow['steps']['2']
        tool_state = _load_tool_state(gx_step, additional_filters=['ReplaceNullWithVarname'])
        self.assertEqual(tool_state['output_selector'], '$output_selector')
        self.assertEqual(tool_state['adapter_options']['action'], 'trim')
        self.assertEqual(tool_state['library']['r1']['cut'], '0')

    def test_flat_filter(self) -> None:
        wf_path = os.path.abspath('./janis_core/tests/data/galaxy/cutadapt_wf.ga')
        gx_workflow = _load_galaxy_workflow(wf_path)
        gx_step = gx_workflow['steps']['2']
        tool_state = _load_tool_state(gx_step, additional_filters=['Flatten'])
        self.assertEqual(tool_state['output_selector'], None)
        self.assertEqual(tool_state['adapter_options.action'], 'trim')
        self.assertEqual(tool_state['library.r1.cut'], '0')




class TestJanisGeneralMapping(unittest.TestCase):
    """
    tests ability to map or generate janis datatypes / selectors etc
    from internal tool objects.
    """
    def setUp(self) -> None:
        datatypes.populate()
        self.tool = MOCK_TOOL_ABRICATE

    def test_to_janis_datatype(self) -> None:
        """
        checks that JanisDatatype objects can be mapped to
        janis DataType objects 
        """
        jtype1 = to_janis_datatype(self.tool.inputs[0])  
        jtype2 = to_janis_datatype(self.tool.inputs[1])
        jtype3 = to_janis_datatype(self.tool.inputs[2])
        jtype4 = to_janis_datatype(self.tool.outputs[0])
        # check janis objects are correct
        self.assertIsInstance(jtype1, File)
        self.assertIsInstance(jtype2, Boolean)
        self.assertIsInstance(jtype3, Array)
        self.assertIsInstance(jtype4, Stdout)
        # check attributes are correct
        self.assertEquals(jtype2.optional, True)
        self.assertIsInstance(jtype3.subtype(), Float)
        self.assertIsInstance(jtype4.subtype, File)

    def test_to_janis_selector(self) -> None:
        """
        checks that OutputComponent objects can generate 
        janis Selector objects
        """
        jsel1 = to_janis_selector(self.tool.outputs[0])
        jsel2 = to_janis_selector(self.tool.outputs[1])
        jsel3 = to_janis_selector(self.tool.outputs[2])
        # check janis objects are correct
        self.assertIsInstance(jsel2, WildcardSelector)
        self.assertIsInstance(jsel3, InputSelector)
        self.assertIsNone(jsel1)
        # check attributes are correct
        self.assertEquals(jsel2.wildcard, 'report.txt')
        self.assertEquals(jsel3.input_to_select, 'fileInput')



class TestJanisToolMapping(unittest.TestCase):
    """
    tests ability to map or generate janis tool objects 
    from internal tool objects.
    MOCK_TOOL is used to test mapping functions.
    """
    def setUp(self) -> None:
        datatypes.populate()
        self.tool = MOCK_TOOL_ABRICATE
    
    def test_to_janis_tool(self) -> None:
        """
        checks that Tool objects can be mapped to
        janis CommandToolBuilder objects 
        """
        jtool = to_janis_tool(self.tool)
        self.assertIsInstance(jtool, CommandToolBuilder)

    def test_to_janis_tool_input(self) -> None:
        """
        checks that InputComponent objects can be mapped to
        janis ToolInput objects 
        """
        jinp1 = to_janis_tool_input(self.tool.inputs[0])  
        jinp2 = to_janis_tool_input(self.tool.inputs[1])
        jinp3 = to_janis_tool_input(self.tool.inputs[2])
        jinp4 = to_janis_tool_input(self.tool.inputs[3])
        # check janis objects are correct
        self.assertIsInstance(jinp1, ToolInput)
        self.assertIsInstance(jinp2, ToolInput)
        self.assertIsInstance(jinp3, ToolInput)
        self.assertIsInstance(jinp4, ToolInput)
        # check attributes are correct
        self.assertEquals(jinp1.tag, 'fileInput')
        self.assertEquals(jinp1.prefix, None)
        self.assertEquals(jinp1.separate_value_from_prefix, True)
        self.assertIsInstance(jinp1.input_type, File)
        
        self.assertEquals(jinp2.tag, 'noheader')
        self.assertEquals(jinp2.prefix, '--noheader')
        self.assertIsInstance(jinp2.input_type, Boolean)
        self.assertEquals(jinp2.input_type.optional, True)

        self.assertEquals(jinp3.tag, 'minid')
        self.assertEquals(jinp3.prefix, '--minid=')
        self.assertEquals(jinp3.separate_value_from_prefix, False)
        self.assertIsInstance(jinp3.input_type, Array)
        
        self.assertEquals(jinp4.tag, 'db')
        self.assertEquals(jinp4.prefix, '--db=')
        self.assertEquals(jinp4.separate_value_from_prefix, False)
        self.assertEquals(jinp4.default, 'resfinder')
        self.assertIsInstance(jinp4.input_type, String)

    def test_to_janis_tool_output(self) -> None:
        """
        checks that OutputComponent objects can be mapped to
        janis ToolOutput objects 
        """
        jout1 = to_janis_tool_output(self.tool.outputs[0])
        jout2 = to_janis_tool_output(self.tool.outputs[1])
        jout3 = to_janis_tool_output(self.tool.outputs[2])
        # check janis objects are correct
        self.assertIsInstance(jout1, ToolOutput)
        self.assertIsInstance(jout2, ToolOutput)
        self.assertIsInstance(jout3, ToolOutput)
        # check attributes are correct
        self.assertIsInstance(jout1.output_type, Stdout)
        self.assertEquals(jout1.tag, 'outReport1')
        self.assertIsNone(jout1.selector)
        
        self.assertIsInstance(jout2.output_type, File)
        self.assertEquals(jout2.tag, 'outReport2')
        self.assertIsInstance(jout2.selector, WildcardSelector)
        self.assertEquals(jout2.selector.wildcard, 'report.txt')
        
        self.assertIsInstance(jout3.output_type, File)
        self.assertEquals(jout3.tag, 'outFileInput')
        self.assertIsInstance(jout3.selector, InputSelector)
        self.assertEquals(jout3.selector.input_to_select, 'fileInput')

    def test_to_janis_metadata(self) -> None:
        """
        checks that ToolMetadata objects can be mapped to
        janis ToolXMLMetadata objects 
        """
        # MOCK_TOOL
        jmeta = to_janis_metadata(self.tool.metadata)
        self.assertIsInstance(jmeta, ToolMetadata)
        self.assertIsNotNone(jmeta.contributors)
        self.assertIsNotNone(jmeta.dateCreated)
        self.assertIsNotNone(jmeta.dateUpdated)
        self.assertIsNotNone(jmeta.citation)
        self.assertIsNotNone(jmeta.documentation)
        self.assertIsNotNone(jmeta.short_documentation)
        self.assertIsNotNone(jmeta.version)




class TestJanisWorkflowMapping(unittest.TestCase):
    """
    tests ability to map or generate janis workflow objects 
    from internal workflow objects
    """
    def setUp(self) -> None:
        datatypes.populate()
        self.workflow = MOCK_WORKFLOW
        self.jworkflow = to_janis_workflow(self.workflow)
        self.jinputs = to_janis_inputs_dict(self.workflow)

    def test_to_janis_workflow(self) -> None:
        self.assertIsInstance(self.jworkflow, WorkflowBuilder)
    
    def test_to_janis_inputs_dict(self) -> None:
        # single input, no value
        self.assertEquals(self.jinputs['inFasta'], None)
        # single input, provided value
        self.workflow.inputs[0].value = 'path/to/file.fasta'
        jinputs = to_janis_inputs_dict(self.workflow)
        self.assertEquals(jinputs['inFasta'], 'path/to/file.fasta')

    def test_janis_metadata(self) -> None:
        self.assertIsInstance(self.jworkflow.metadata, WorkflowMetadata)

    def test_janis_inputs(self) -> None:
        target_inputs = {
            'inFasta', 
            'abricate_noheader', 
            'abricate_minid', 
            'abricate_db'
        }
        actual_inputs = set(self.jworkflow.input_nodes.keys())
        self.assertEquals(target_inputs, actual_inputs)

        jinp = self.jworkflow.input_nodes['inFasta']
        self.assertIsInstance(jinp, InputNode)
        self.assertIsInstance(jinp.datatype, File)

    def test_janis_outputs(self) -> None:
        target_outputs = {'abricate_outReport1'}
        actual_outputs = set(self.jworkflow.output_nodes.keys())
        self.assertEquals(target_outputs, actual_outputs)

        jout = self.jworkflow.output_nodes['abricate_outReport1']
        self.assertIsInstance(jout, OutputNode)
        self.assertIsInstance(jout.datatype, Stdout)
        self.assertIsInstance(jout.datatype.subtype, File)
        self.assertEquals(jout.doc.doc, 'report file')

    def test_janis_steps(self) -> None:
        target_steps = {'abricate'}
        actual_steps = set(self.jworkflow.step_nodes.keys())
        self.assertEquals(target_steps, actual_steps)

        # basic object checks
        jstep = self.jworkflow.step_nodes['abricate']
        self.assertIsInstance(jstep, StepNode)
        self.assertIsInstance(jstep.tool, CommandTool)

        # sources (tool input values)
        self.assertIn('fileInput', jstep.sources)
        self.assertIn('noheader', jstep.sources)
        self.assertIn('minid', jstep.sources)
        self.assertIn('db', jstep.sources)
        
        # scatter 
        self.assertEquals(jstep.scatter.fields, ['minid'])
        self.assertEquals(jstep.scatter.method, ScatterMethod.dot)




class TestDatatypeInference(unittest.TestCase):
    """
    tests the datatype which is assigned to an entity
    """
    def setUp(self) -> None:
        datatypes.populate()

    def test_positional(self) -> None:
        self.assertEquals(datatypes.get(MOCK_POSITIONAL1), file_t)
    
    def test_flag(self) -> None:
        self.assertEquals(datatypes.get(MOCK_FLAG1), bool_t)
    
    def test_option(self) -> None:
        self.assertEquals(datatypes.get(MOCK_OPTION2), string_t)
    
    def test_outputs(self) -> None:
        self.assertEquals(datatypes.get(MOCK_REDIRECT_OUTPUT), file_t)
    
    def test_workflow_input(self) -> None:
        self.assertEquals(datatypes.get(MOCK_WORKFLOW_INPUT1), file_t)
    
    # def test_option_typestring(self) -> None:
    #     raise NotImplementedError()




class TestContainerFetching(unittest.TestCase):
    """
    tests biocontainer fetching given a tool id, tool version, and a Requirement()
    """

    def setUp(self) -> None:
        """
        creates a temp container cache for use. the tearDown() method will remove this after tests have run. 
        """
        self.temp_cache_dir = '/tmp/container_cache.json'
        with open(self.temp_cache_dir, 'w') as fp:
            fp.write('{}')
    
    def test_1(self) -> None:
        result = containers.fetch_container(QUERY1)
        self.assertEquals(result, QUERY1_EXPECTED_RESULT)
    
    def test_2(self) -> None:
        result = containers.fetch_container(QUERY2)
        self.assertEquals(result, QUERY2_EXPECTED_RESULT)
    
    def test_3(self) -> None:
        result = containers.fetch_container(QUERY3)
        self.assertEquals(result, QUERY3_EXPECTED_RESULT)
    
    def tearDown(self) -> None:
        if os.path.exists(self.temp_cache_dir):
            os.remove(self.temp_cache_dir)




def read_cmd(path: str) -> str:
    tree = et.parse(path)
    root = tree.getroot()
    assert(root.text)
    return root.text

def read_step_inputs(path: str) -> dict[str, Any]:
    with open(path, 'r') as fp:
        step = json.load(fp)
    step['tool_state'] = load_tool_state(step)
    return step['tool_state']

class TestSectionalCheetah(unittest.TestCase):

    def test_unicycler(self):
        vanilla = read_cmd(UNICYCLER_VANILLA_PATH)
        reference = read_cmd(UNICYCLER_TEMPLATED_PATH)
        inputs = read_step_inputs(UNICYCLER_INPUTS_PATH)
        templated = sectional_evaluate(vanilla, inputs)
        self.assertEquals(reference, templated)




def get_cmd(path: str) -> str:
    tree = et.parse(path)
    root = tree.getroot()
    assert(root.text)
    return root.text

class TestAliases(unittest.TestCase):

    def test_resolve_fastqc(self):
        raw_path = './janis_core/tests/data/command/manipulation/aliases/fastqc/fastqc_command.xml'
        ref_path = './janis_core/tests/data/command/manipulation/aliases/fastqc/fastqc_command_resolved.xml'
        raw_cmd = get_cmd(raw_path)
        ref_cmd = get_cmd(ref_path)
        res_cmd = resolve_aliases(raw_cmd)
        self.assertEquals(ref_cmd, res_cmd)
    
    def test_resolve_unicycler(self):
        raw_path = './janis_core/tests/data/command/manipulation/aliases/unicycler/unicycler_command.xml'
        ref_path = './janis_core/tests/data/command/manipulation/aliases/unicycler/unicycler_command_resolved.xml'
        raw_cmd = get_cmd(raw_path)
        ref_cmd = get_cmd(ref_path)
        res_cmd = resolve_aliases(raw_cmd)
        self.assertEquals(ref_cmd, res_cmd)



class TestFromGalaxy(unittest.TestCase):

    def test_ingest_abricate_tool(self) -> None:
        filepath = os.path.abspath('./janis_core/tests/data/galaxy/abricate/abricate.xml')
        jtool = ingest_galaxy(filepath)
        assert(isinstance(jtool, CommandTool))
        
        self.assertEquals(len(jtool.inputs()), 5)
        self.assertEquals(len(jtool.outputs()), 1)
        self.assertEquals(jtool.base_command(), ['abricate'])

    def test_ingest_cutadapt_wf(self) -> None:
        filepath = os.path.abspath('./janis_core/tests/data/galaxy/cutadapt_wf.ga')
        jworkflow = ingest_galaxy(filepath)
        assert(isinstance(jworkflow, WorkflowBuilder))

        self.assertEquals(len(jworkflow.step_nodes), 6)
        self.assertEquals(len(jworkflow.output_nodes), 19)
        self.assertIn('inForwardReads', jworkflow.input_nodes)
        self.assertIn('inReverseReads', jworkflow.input_nodes)
        self.assertIn('inLongReads', jworkflow.input_nodes)
    
    def test_translate_cutadapt_wf_nextflow(self) -> None:
        srcfmt = 'galaxy'
        destfmt = 'nextflow'
        filepath = os.path.abspath('./janis_core/tests/data/galaxy/cutadapt_wf.ga')
        run(filepath, srcfmt, destfmt)
    
    def test_ingest_unicycler_assembly(self) -> None:
        filepath = os.path.abspath('./janis_core/tests/data/galaxy/unicycler_assembly.ga')
        jworkflow = ingest_galaxy(filepath)
        assert(isinstance(jworkflow, WorkflowBuilder))

        self.assertEquals(len(jworkflow.step_nodes), 6)
        self.assertEquals(len(jworkflow.output_nodes), 19)
        self.assertIn('inForwardReads', jworkflow.input_nodes)
        self.assertIn('inReverseReads', jworkflow.input_nodes)
        self.assertIn('inLongReads', jworkflow.input_nodes)

