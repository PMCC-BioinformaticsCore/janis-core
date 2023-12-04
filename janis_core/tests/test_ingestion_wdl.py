

import os 
import unittest
from janis_core import CommandToolBuilder, WorkflowBuilder, InputSelector, StringFormatter, FilterNullOperator
from janis_core import ScatterDescription, ScatterMethod
from janis_core import IsDefined, NotOperator, Operator, FirstOperator
from janis_core.ingestion.wdl.parsing import parse_task
from janis_core.ingestion.wdl.parsing import parse_container_requirement
from janis_core.ingestion.wdl.parsing import parse_cpus_requirement
from janis_core.ingestion.wdl.parsing import parse_memory_requirement
from janis_core.ingestion.wdl.parsing import parse_disk_requirement
from janis_core.ingestion.wdl import WdlParser
from janis_core.ingestion import ingest
from janis_core import (
    File, 
    String, 
    Int,
    Boolean, 
    Array,
)
import WDL

WDL_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/wdl')
from janis_core import settings 
from janis_core.messages import configure_logging

def _do_setup_unsafe() -> None:
    configure_logging()
    settings.ingest.SAFE_MODE = False
    settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = False
    settings.graph.ALLOW_INCOMPATIBLE_TYPES = False
    settings.graph.ALLOW_INCORRECT_NUMBER_OF_SOURCES = False
    settings.graph.ALLOW_NON_ARRAY_SCATTER_INPUT = False
    settings.graph.ALLOW_UNKNOWN_SCATTER_FIELDS = False
    settings.graph.ALLOW_UNKNOWN_SOURCE = False
    settings.validation.STRICT_IDENTIFIERS = True
    settings.validation.VALIDATE_STRINGFORMATTERS = True
    settings.general.ALLOW_EMPTY_CONTAINER = False

def _do_setup_safe() -> None:
    configure_logging()
    settings.ingest.SAFE_MODE = True
    settings.datatypes.ALLOW_UNPARSEABLE_DATATYPES = True
    settings.graph.ALLOW_INCOMPATIBLE_TYPES = True
    settings.graph.ALLOW_INCORRECT_NUMBER_OF_SOURCES = True
    settings.graph.ALLOW_NON_ARRAY_SCATTER_INPUT = True
    settings.graph.ALLOW_UNKNOWN_SCATTER_FIELDS = True
    settings.graph.ALLOW_UNKNOWN_SOURCE = True
    settings.validation.STRICT_IDENTIFIERS = False
    settings.validation.VALIDATE_STRINGFORMATTERS = False
    settings.general.ALLOW_EMPTY_CONTAINER = True


def _simple_lines(text: str) -> list[str]:
    if not isinstance(text, str):
        text = str(text)
    lines = text.split('\n')
    lines = [ln.strip(' \t') for ln in lines]
    out = []
    for ln in lines:
        if '#' in ln:
            ln = ln.split('#')[0]
        out.append(ln)
    out = [ln for ln in out if ln != '']
    return out

##############
### BASICS ###
##############

class TestBasicFunctionality(unittest.TestCase):
    
    def setUp(self) -> None:
        _do_setup_unsafe()
    
    def test_tool_rename(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/rename_tool.wdl'
        cmdtool = ingest(filepath, 'wdl')
        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(len(cmdtool._inputs), 2)
        self.assertEqual(len(cmdtool._outputs), 1)
    
    def test_tool_io(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/io_tool.wdl'
        cmdtool = ingest(filepath, 'wdl')
        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(len(cmdtool._inputs), 7)
        self.assertEqual(len(cmdtool._outputs), 5)
    
    def test_tool_bwa(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/bwa.wdl'
        cmdtool = ingest(filepath, 'wdl')
        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(len(cmdtool._inputs), 14)
        self.assertEqual(len(cmdtool._outputs), 2)
    
    def test_tool_fastqc(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/fastqc.wdl'
        cmdtool = ingest(filepath, 'wdl')
        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(len(cmdtool._inputs), 21)
        self.assertEqual(len(cmdtool._outputs), 5)
    
    def test_tool_trim_adapters(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/TrimAdapters.wdl'
        cmdtool = ingest(filepath, 'wdl')
        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(len(cmdtool._inputs), 8)
        self.assertEqual(len(cmdtool._outputs), 3)
    
    def test_workflow_io(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/io_wf.wdl'
        wf = ingest(filepath, 'wdl')
        self.assertIsInstance(wf, WorkflowBuilder)
        self.assertEqual(wf.id(), 'main')
        self.assertEqual(len(wf.input_nodes), 7)
        self.assertEqual(len(wf.step_nodes), 1)
        self.assertEqual(len(wf.output_nodes), 5)

    def test_workflow_atac(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/ATAC.wdl'
        wf = ingest(filepath, 'wdl')
        self.assertIsInstance(wf, WorkflowBuilder)
        self.assertEqual(wf.id(), 'ATAC')
        self.assertEqual(len(wf.input_nodes), 17)
        self.assertEqual(len(wf.step_nodes), 13)
        self.assertEqual(len(wf.output_nodes), 4)
    



class TestRequirements(unittest.TestCase):

    def setUp(self) -> None:
        _do_setup_unsafe()
        self.mocktool = CommandToolBuilder(
            tool='testing',
            version='DEV',
            container='',
            base_command=None,
            inputs=[],
            outputs=[]
        )
    
    def test_container1(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/requirements/annotsv_filter.wdl'
        task = WDL.load(filepath).tasks[0]
        actual = parse_container_requirement(task, self.mocktool)
        expected = 'python:3'
        self.assertEqual(actual, expected)
    
    def test_container2(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/requirements/bwa.wdl'
        task = WDL.load(filepath).tasks[0]
        actual = parse_container_requirement(task, self.mocktool)
        expected = 'quay.io/biocontainers/mulled-v2-ad317f19f5881324e963f6a6d464d696a2825ab6:c59b7a73c87a9fe81737d5d628e10a3b5807f453-0'
        self.assertEqual(actual, expected)
    
    def test_cpus1(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/requirements/TrimAdapters.wdl'
        task = WDL.load(filepath).tasks[0]
        actual = parse_cpus_requirement(task, self.mocktool)
        expected = 1
        self.assertEqual(actual, expected)
    
    def test_cpus2(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/requirements/bwa.wdl'
        task = WDL.load(filepath).tasks[0]
        actual = parse_cpus_requirement(task, self.mocktool)
        self.assertIsInstance(actual, InputSelector)
        self.assertEqual(actual.input_to_select, 'threads')
    
    def test_mem1(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/requirements/annotsv_filter.wdl'
        task = WDL.load(filepath).tasks[0]
        actual = parse_memory_requirement(task, self.mocktool)
        expected = 4.0
        self.assertEqual(actual, expected)
    
    def test_mem2(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/requirements/bwa.wdl'
        task = WDL.load(filepath).tasks[0]
        actual = parse_memory_requirement(task, self.mocktool)
        self.assertIsInstance(actual, StringFormatter)
        self.assertEqual(actual._format, '{TOKEN1}G')
        inner = list(actual.kwargs.values())[0]
        self.assertIsInstance(inner, FirstOperator)
        inner1 = inner.args[0][0]
        self.assertIsInstance(inner1, InputSelector)
        self.assertEqual(inner1.input_to_select, 'memoryGb')
        inner2 = inner.args[0][1]
        self.assertIsInstance(inner2, InputSelector)
        self.assertEqual(inner2.input_to_select, 'estimatedMemoryGb')
    
    def test_disk1(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/requirements/alignment_metrics.wdl'
        task = WDL.load(filepath).tasks[0]
        actual = parse_disk_requirement(task, self.mocktool)
        expected = 150
        self.assertEqual(actual, expected)
    
    def test_disk2(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/requirements/annotsv_filter.wdl'
        task = WDL.load(filepath).tasks[0]
        actual = parse_disk_requirement(task, self.mocktool)
        self.assertIsInstance(actual, StringFormatter)
        expected = 'local-disk {inputs.space_needed_gb} SSD'
        self.assertEqual(str(actual), expected)
    


###############################
### DATATYPES / EXPRESSIONS ###
###############################

class TestDatatypes(unittest.TestCase):
    
    def setUp(self) -> None:
        _do_setup_unsafe()


class TestExpressions(unittest.TestCase):
    
    def setUp(self) -> None:
        _do_setup_unsafe()


class TestLexer(unittest.TestCase):
    
    def setUp(self) -> None:
        _do_setup_unsafe()



###############
### COMMAND ###
###############

class TestNativeSimpleCommandParser(unittest.TestCase):
    
    def setUp(self) -> None:
        _do_setup_unsafe()
        settings.ingest.wdl.COMMAND_PARSER = 'native_simple'

    @unittest.skip('TODO implement')
    def test_rename_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/rename_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        cmdtool = parse_task(task)

        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(cmdtool._base_command, ['cp'])
        self.assertEqual(len(cmdtool._inputs), 2)
        self.assertEqual(cmdtool._inputs[0].id(), 'sourceFile')
        self.assertIsInstance(cmdtool._inputs[0].input_type, File)
        self.assertEqual(cmdtool._inputs[0].position, 1)
        self.assertEqual(cmdtool._inputs[1].id(), 'targetFilename')
        self.assertIsInstance(cmdtool._inputs[1].input_type, String)
        self.assertEqual(cmdtool._inputs[1].position, 2)
    
    @unittest.skip('TODO implement')
    def test_io_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/io_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        cmdtool = parse_task(task)
        
        ### basics ###
        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(cmdtool._base_command, ['sh', 'script.sh'])
        self.assertEqual(len(cmdtool._inputs), 7)
        self.assertEqual(len(cmdtool._outputs), 5)
        
        ### inputs ###
        # inInt
        tinp = cmdtool._inputs[0]
        self.assertEqual(tinp.id(), 'inInt')
        self.assertIsInstance(tinp.input_type, Int)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 1)
        self.assertIsNone(tinp.prefix)
        # inStr
        tinp = cmdtool._inputs[1]
        self.assertEqual(tinp.id(), 'inStr')
        self.assertIsInstance(tinp.input_type, String)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 2)
        self.assertIsNone(tinp.prefix)
        # inBool
        tinp = cmdtool._inputs[2]
        self.assertEqual(tinp.id(), 'inBool')
        self.assertIsInstance(tinp.input_type, Boolean)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 3)
        self.assertEqual(tinp.prefix, '--flag1')
        # inFile
        tinp = cmdtool._inputs[3]
        self.assertEqual(tinp.id(), 'inFile')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 4)
        self.assertEqual(tinp.prefix, '--in-file=')
        self.assertFalse(tinp.separate_value_from_prefix)
        # inFileOpt
        tinp = cmdtool._inputs[4]
        self.assertEqual(tinp.id(), 'inFileOpt')
        self.assertIsInstance(tinp.input_type, File)
        self.assertTrue(tinp.input_type.optional)
        self.assertEqual(tinp.position, 5)
        self.assertIsNone(tinp.prefix)
        # inFileArr
        tinp = cmdtool._inputs[5]
        self.assertEqual(tinp.id(), 'inFileArr')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 6)
        self.assertIsNone(tinp.prefix)
        # inSecondary
        tinp = cmdtool._inputs[6]
        self.assertEqual(tinp.id(), 'inSecondary')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 7)
        self.assertIsNone(tinp.prefix)

        ### outputs ###
        tout = cmdtool._outputs[0]

    @unittest.skip('TODO implement')
    def test_fastqc_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/fastqc.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)
    
    @unittest.skip('TODO implement')
    def test_bwa_mem(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/bwa_mem.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)


class TestNativeArgumentsCommandParser(unittest.TestCase):
    
    def setUp(self) -> None:
        _do_setup_unsafe()
        settings.ingest.wdl.COMMAND_PARSER = 'native_arguments'

    @unittest.skip('TODO implement')
    def test_rename_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/rename_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        cmdtool = parse_task(task)

        self.assertIsInstance(cmdtool, CommandToolBuilder)
        print(cmdtool._files_to_create['script.sh'])
        self.assertEqual(cmdtool._base_command, ['sh', 'script.sh'])
        self.assertEqual(len(cmdtool._inputs), 2)
        self.assertEqual(cmdtool._inputs[0].id(), 'sourceFile')
        self.assertIsInstance(cmdtool._inputs[0].input_type, File)
        self.assertEqual(cmdtool._inputs[1].id(), 'targetFilename')
        self.assertIsInstance(cmdtool._inputs[1].input_type, String)
    
    @unittest.skip('TODO implement')
    def test_io_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/io_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        cmdtool = parse_task(task)
        
        ### basics ###
        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(cmdtool._base_command, ['sh', 'script.sh'])
        self.assertEqual(len(cmdtool._inputs), 7)
        self.assertEqual(len(cmdtool._outputs), 5)
        
        ### inputs ###
        # inInt
        tinp = cmdtool._inputs[0]
        self.assertEqual(tinp.id(), 'inInt')
        self.assertIsInstance(tinp.input_type, Int)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 1)
        self.assertIsNone(tinp.prefix)
        # inStr
        tinp = cmdtool._inputs[1]
        self.assertEqual(tinp.id(), 'inStr')
        self.assertIsInstance(tinp.input_type, String)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 2)
        self.assertIsNone(tinp.prefix)
        # inBool
        tinp = cmdtool._inputs[2]
        self.assertEqual(tinp.id(), 'inBool')
        self.assertIsInstance(tinp.input_type, Boolean)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 3)
        self.assertEqual(tinp.prefix, '--flag1')
        # inFile
        tinp = cmdtool._inputs[3]
        self.assertEqual(tinp.id(), 'inFile')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 4)
        self.assertEqual(tinp.prefix, '--in-file=')
        self.assertFalse(tinp.separate_value_from_prefix)
        # inFileOpt
        tinp = cmdtool._inputs[4]
        self.assertEqual(tinp.id(), 'inFileOpt')
        self.assertIsInstance(tinp.input_type, File)
        self.assertTrue(tinp.input_type.optional)
        self.assertEqual(tinp.position, 5)
        self.assertIsNone(tinp.prefix)
        # inFileArr
        tinp = cmdtool._inputs[5]
        self.assertEqual(tinp.id(), 'inFileArr')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 6)
        self.assertIsNone(tinp.prefix)
        # inSecondary
        tinp = cmdtool._inputs[6]
        self.assertEqual(tinp.id(), 'inSecondary')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 7)
        self.assertIsNone(tinp.prefix)

        ### outputs ###
        tout = cmdtool._outputs[0]

    @unittest.skip('TODO implement')
    def test_fastqc_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/fastqc.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)
    
    @unittest.skip('TODO implement')
    def test_bwa_mem(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/bwa_mem.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)




class TestShellCommandParser(unittest.TestCase):
    
    def setUp(self) -> None:
        _do_setup_unsafe()
        settings.ingest.wdl.COMMAND_PARSER = 'shell'
    
    def test_rename_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/rename_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        cmdtool = parse_task(task)

        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(cmdtool._base_command, ['sh', 'script.sh'])
        self.assertEqual(len(cmdtool._inputs), 2)
        self.assertEqual(cmdtool._inputs[0].id(), 'sourceFile')
        self.assertIsInstance(cmdtool._inputs[0].input_type, File)
        self.assertEqual(cmdtool._inputs[1].id(), 'targetFilename')
        self.assertIsInstance(cmdtool._inputs[1].input_type, String)

        actual = _simple_lines(cmdtool._files_to_create['script.sh'])
        expected = [
            'set -e',
            'cp \\',
            '{inputs.sourceFile} \\',
            '{inputs.targetFilename} \\',
        ]
        self.assertListEqual(actual, expected)
    
    def test_io_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/io_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        cmdtool = parse_task(task)
        
        self.assertIsInstance(cmdtool, CommandToolBuilder)
        self.assertEqual(cmdtool._base_command, ['sh', 'script.sh'])
        self.assertEqual(len(cmdtool._inputs), 7)
        self.assertEqual(len(cmdtool._outputs), 5)

        actual = _simple_lines(cmdtool._files_to_create['script.sh'])
        expected = [
            'set -e',
            'echo \\',
            '{inputs.inInt} \\',
            '{inputs.inStr} \\',
            '{inputs.inBool} \\',
            '--in-file={inputs.inFile} \\',
            '{inputs.inFileOpt} \\',
            '{inputs.inFileArr} \\',
            '{inputs.inSecondary} \\',
            '> stdout.txt',
        ]
        self.assertListEqual(actual, expected)

    def test_fastqc_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/fastqc.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)
    
    def test_bwa_mem(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/bwa_mem.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)


class TestPlumbing(unittest.TestCase):
    
    def setUp(self) -> None:
        _do_setup_unsafe()

    def test_step_inputs1(self) -> None:
        settings.ingest.wdl.COMMAND_PARSER = 'shell'
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/subworkflows/create_alignment_from_families_files.wdl'
        wf = ingest(filepath, 'wdl')
        self.assertIsInstance(wf, WorkflowBuilder)
        actual_steps = list(wf.step_nodes.keys())
        expected_steps = [
            'SepareChunksFastqString',
            'RunBwaAlignment',
            'MergeBams',
        ]
        self.assertEqual(actual_steps, expected_steps)
        
        # FIRST STEP 
        step = wf.step_nodes['SepareChunksFastqString']
        self.assertIsNone(step.scatter)
        expected_sources = {
            'families_info': 'inputs.families_info',
            'chunk_size': 'inputs.chunk_size',
        }
        for tinput_id, src in step.sources.items():
            actual_src = str(src.source_map[0].source)
            expected_src = expected_sources[tinput_id]
            self.assertEqual(actual_src, expected_src)
        
        # SECOND STEP 
        step = wf.step_nodes['RunBwaAlignment']
        expected_sources = {
            'sampleName': 'SepareChunksFastqString.chunks[1]',
            'reads': 'SepareChunksFastqString.chunks[0]',
            'libraries': 'SepareChunksFastqString.chunks[2]',
            'references': 'inputs.references',
            'max_cores': 'inputs.max_cores',
            'rm_dupli': 'inputs.rm_dupli',
        }
        for tinput_id, src in step.sources.items():
            actual_src = str(src.source_map[0].source)
            expected_src = expected_sources[tinput_id]
            self.assertEqual(actual_src, expected_src)
        
        # THIRD STEP 
        step = wf.step_nodes['MergeBams']
        self.assertIsNone(step.scatter)
        expected_sources = {
            'bam_files': 'flatten(RunBwaAlignment.bam)',
        }
        for tinput_id, src in step.sources.items():
            actual_src = str(src.source_map[0].source)
            expected_src = expected_sources[tinput_id]
            self.assertEqual(actual_src, expected_src)

    def test_conditional_deps(self) -> None:
        settings.ingest.wdl.COMMAND_PARSER = 'shell'
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/subworkflows/create_alignment_from_read_simulations.wdl'
        wf = ingest(filepath, 'wdl')
        self.assertIsInstance(wf, WorkflowBuilder)

        cond1 = '!(isdefined(inputs.sequencing))'
        cond2 = 'isdefined(inputs.sequencing)'
        cond3 = '((inputs.sequencing == WGS) or (inputs.sequencing == exome))'
        cond4 = '((inputs.sequencing == sdRAD) or (inputs.sequencing == ddRAD))'

        expected_conditions = {
            'GenerateAlternativeGenome': cond1,
            'CreatePedigreeSimulatorInputs': cond1,
            'Vcf2PedigreeSimulator': cond2,
            'RunPedigreeSimulator': None,
            'ConvertPedigreeSimulationToVcf': None,
            'GenerateSampleNames': None,
            'SimuscopProfile': cond3,
            'SimuscopSimulation': cond3,
            'RADinitioSimulation': cond4,
            'SepareChunksFastq': None,
            'RunBwaAlignmentSimu': None,
            'MergeBams': None,
        }

        for stepid, step in wf.step_nodes.items():
            expected = expected_conditions[stepid]
            if expected is not None:
                self.assertIsInstance(step.when, Operator)
                self.assertEqual(str(step.when), expected)
            else:
                self.assertIsNone(step.when)
    
    def test_scatter_deps(self) -> None:
        settings.ingest.wdl.COMMAND_PARSER = 'shell'
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/subworkflows/create_alignment_from_families_files.wdl'
        wf = ingest(filepath, 'wdl')
        self.assertIsInstance(wf, WorkflowBuilder)
        
        step = wf.step_nodes['RunBwaAlignment']
        self.assertIsInstance(step.scatter, ScatterDescription)
        self.assertEqual(step.scatter.method, ScatterMethod.dot)
        self.assertSetEqual(set(step.scatter.fields), set(['sampleName', 'reads', 'libraries']))
    
    def test_nested_deps(self) -> None:
        settings.ingest.wdl.COMMAND_PARSER = 'shell'
        filepath = f'{WDL_TESTDATA_PATH}/Reads2Map/subworkflows/create_alignment_from_read_simulations.wdl'
        wf = ingest(filepath, 'wdl')
        self.assertIsInstance(wf, WorkflowBuilder)

        step = wf.step_nodes['SimuscopSimulation']
        self.assertIsInstance(step.scatter, ScatterDescription)
        self.assertEqual(step.scatter.method, ScatterMethod.dot)
        self.assertSetEqual(set(step.scatter.fields), set(['sampleName']))
        expected = '((inputs.sequencing == WGS) or (inputs.sequencing == exome))'
        self.assertIsInstance(step.when, Operator)
        self.assertEqual(str(step.when), expected)