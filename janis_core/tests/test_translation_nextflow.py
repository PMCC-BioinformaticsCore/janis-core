import unittest
import regex as re

from janis_core.tests.testtools import (
    InputQualityTestTool,
    BasicTestTool,
    FastqcTestTool,
    BwaMemTestTool,
    GridssTestTool,
)

from janis_core.tests.testworkflows import (

    # basics
    BasicIOTestWF,
    WildcardSelectorOutputTestWF,
    InputSelectorTestWF,
    StepInputsTestWF,
    StepInputsWFInputTestWF,
    StepInputsStaticTestWF,
    StepInputsPartialStaticTestWF,
    StepInputsMinimalTestWF,
    StepConnectionsTestWF,
    DirectivesTestWF,

    # arrays
    ArrayIOTestWF,
    ArrayIOExtrasTestWF,
    ArrayStepInputsTestWF,
    ArrayStepConnectionsTestWF,

    # scatter
    BasicScatterTestWF,
    ChainedScatterTestWF,
    ScatterDotTestWF,
    ScatterCrossTestWF,

    # secondaries
    SecondariesTestWF,
    ComprehensiveScatterTestWF,

    # combos
    ScatterSecondariesTestWF,

    # additional features
    StepInputExpressionTestWF,
    ConditionStepTestWF,
    AliasSelectorTestWF,
    ArraysOfSecondaryFilesOutputsTestWF,
    ForEachTestWF,
    IndexOperatorTestWF,
    StringFormatterTestWF,

    # codetools
    InputsPythonToolTestWF,
    OutputsPythonToolTestWF,

    # specific workflows
    AssemblyTestWF,
    SubworkflowTestWF,
    Subworkflow2TestWF,
    Subworkflow3TestWF,
    DataSourceTestWF,
    FilenameTestWF1,
    FilenameTestWF2,
    OutputCollectionTestWF,
    UnwrapTestWF,
    NamingTestWF,
    PlumbingTypeMismatchTestWF,
    EntityTraceTestWF,
    FilePairsTestWF,
    ProcessInputsTestWF,
    OrderingTestWF,
    PlumbingEdgeCaseTestWF,
    OptionalTestWF,
    DuplicateTasksTestWF,
    MinimalTaskInputsTestWF1,
    MinimalTaskInputsTestWF2,
    MinimalTaskInputsTestWF3,
    MinimalTaskInputsTestWF4,
    MinimalTaskInputsTestWF5,
    MinimalTaskInputsTestWF6,
) 

from janis_core import (
    Workflow,
    CommandTool,
    PythonTool,
    InputSelector,
    StringFormatter,
    JoinOperator,
)

from janis_core.translations import translate
from janis_core.translations import NextflowTranslator as translator
from janis_core.translations import nextflow

from janis_core.translations.nextflow.variables import VariableManager

from janis_core import (
    String, 
    Int, 
    Float, 
    File, 
    Boolean,
    UnionType,
    Array,
)

from janis_bioinformatics.data_types.fasta import Fasta, FastaDict, FastaWithIndexes
from janis_bioinformatics.data_types.fastq import Fastq, FastqGzPair
from janis_bioinformatics.data_types.bam import Bam, BamBai
from janis_bioinformatics.data_types.sam import Sam
from janis_bioinformatics.data_types.cram import Cram

from janis_core.translations.nextflow.nfgen_utils import to_groovy
from janis_core import settings


### helper functions

def reset_globals() -> None:
    # nextflow specific
    settings.translate.nextflow.MODE = 'workflow'
    settings.translate.nextflow.MINIMAL_PROCESS = True

    # general
    settings.validation.STRICT_IDENTIFIERS = True 
    settings.translate.ALLOW_EMPTY_CONTAINER = True 
    settings.translate.MERGE_RESOURCES = False
    settings.translate.RENDER_COMMENTS = True 
    settings.translate.SHOULD_VALIDATE = False
    settings.translate.SHOULD_ZIP = False
    settings.translate.TO_DISK = False
    settings.translate.TO_CONSOLE = True 
    settings.translate.TOOL_TO_CONSOLE = False
    settings.translate.WITH_CONTAINER = True 
    settings.translate.WITH_RESOURCE_OVERRIDES = False
    settings.translate.WRITE_INPUTS_FILE = True 
    settings.translate.CONTAINER_OVERRIDE = None
    settings.translate.ADDITIONAL_INPUTS = None
    settings.translate.HINTS = None
    settings.translate.MAX_CORES = None           
    settings.translate.MAX_DURATION = None           
    settings.translate.MAX_MEM = None

    nextflow.task_inputs.clear()

def do_preprocessing_workflow(wf: Workflow) -> None:
    nextflow.preprocessing.populate_task_inputs(wf, wf)

def do_preprocessing_tool(tool: CommandTool | PythonTool) -> None:
    pass

def split_task_call_to_lines(call: str) -> list[str]:
    lines = call.split('\n')                             # split
    lines = lines[1:]                                    # remove task name & opening brace line
    lines = [ln for ln in lines if not ln in ('', ')')]  # remove empty lines, closing brace line
    lines = [ln.strip(' ,') for ln in lines]             # strip indents & commas
    return lines

def split_to_lines(call: str) -> list[str]:
    lines = call.split('\n')                        # split
    lines = [ln.strip(' ') for ln in lines]         # strip indents
    lines = [ln for ln in lines if not ln == '']    # remove empty lines
    return lines



### helper classes

class DataTypeWithSecondary(File):
    @staticmethod
    def name() -> str:
        return "test_secondary"

    @staticmethod
    def secondary_files():
        return [".txt", ".csv"]


class DataTypeNoSecondary(File):
    @staticmethod
    def name() -> str:
        return "test_no_secondary"




### test classes

class TestTaskInputs(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    # no subworkflows
    def test_one_call(self) -> None:
        wf = MinimalTaskInputsTestWF1()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    def test_two_calls(self) -> None:
        wf = MinimalTaskInputsTestWF2()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int2',
            'val in_int3',
            'val in_str4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    def test_three_calls(self) -> None:
        wf = MinimalTaskInputsTestWF3()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int2',
            'val in_int3',
            'val in_str4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    # subworkflows
    def test_one_call_sub(self) -> None:
        # TODO improve? 
        wf = MinimalTaskInputsTestWF4()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        step = step.tool.step_nodes["stp1"]
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)

    def test_two_calls_sub(self) -> None:
        # TODO improve? 
        # main wf process
        wf = MinimalTaskInputsTestWF5()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int1',
            'val in_int2',
            'val in_str1',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)

        # sub wf process
        step = wf.step_nodes["stp2"]
        scope = nextflow.Scope()
        scope.update(step)
        step = step.tool.step_nodes["stp1"]
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int1',
            'val in_int2',
            'val in_str1',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    def test_three_calls_sub(self) -> None:
        # TODO improve? 
        wf = MinimalTaskInputsTestWF6()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int1',
            'val in_int2',
            'val in_int3',
            'val in_str1',
            'val in_str2',
            'val in_str3',
            'val in_str4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
        
        step = wf.step_nodes["stp2"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int1',
            'val in_int2',
            'val in_int3',
            'val in_str1',
            'val in_str2',
            'val in_str3',
            'val in_str4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)

        # sub wf process
        step = wf.step_nodes["stp3"]
        scope = nextflow.Scope()
        scope.update(step)
        step = step.tool.step_nodes["stp1"]
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int1',
            'val in_int2',
            'val in_int3',
            'val in_str1',
            'val in_str2',
            'val in_str3',
            'val in_str4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)




class TestPreprocessingTaskInputs(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_main_wf(self) -> None:
        wf = MinimalTaskInputsTestWF1()
        do_preprocessing_workflow(wf)

        # main wf
        actual_task_inputs = nextflow.task_inputs.task_inputs(wf.id())
        expected_task_inputs: set[str] = set()
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
        actual_param_inputs = nextflow.task_inputs.param_inputs(wf.id())
        expected_param_inputs = {'inFile', 'inStr1', 'inInt1'}
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
        actual_static_inputs = nextflow.task_inputs.static_inputs(wf.id())
        expected_static_inputs: set[str] = set()
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(wf.id())
        expected_ignored_inputs = {
            'inStr2',
            'inInt3',
            'inInt2',
            'stp1_inInt1',
            'inStr3',
        }
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)

    # no subworkflows
    def test_one_call(self) -> None:
        wf = MinimalTaskInputsTestWF1()
        do_preprocessing_workflow(wf)
        tool = wf.step_nodes['stp1'].tool

        actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
        expected_task_inputs = {'inFile'}
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
        actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
        expected_param_inputs = {'inStr1', 'inInt2'}
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
        actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
        expected_static_inputs = {'inInt1'}
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        expected_ignored_inputs = {
            'inStr2',
            'inStr3',
            'inStr4',
            'inInt3',
            'inInt4',
        }
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)
    
    def test_two_calls(self) -> None:
        wf = MinimalTaskInputsTestWF2()
        do_preprocessing_workflow(wf)
        tool = wf.step_nodes['stp1'].tool

        actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
        expected_task_inputs = {'inFile', 'inStr2', 'inStr3', 'inStr4', 'inInt2', 'inInt3'}
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
        actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
        expected_param_inputs = {'inStr1'}
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
        actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
        expected_static_inputs = {'inInt1'}
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        expected_ignored_inputs = {'inInt4'}
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)
    
    def test_three_calls(self) -> None:
        wf = MinimalTaskInputsTestWF3()
        do_preprocessing_workflow(wf)
        tool = wf.step_nodes['stp1'].tool

        actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
        expected_task_inputs = {
            'inFile', 
            'inStr1', 
            'inStr2', 
            'inStr3', 
            'inStr4', 
            'inInt1', 
            'inInt2', 
            'inInt3',
            'inInt4',
        }
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
        actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
        expected_param_inputs = set()
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
        actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
        expected_static_inputs = set()
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        expected_ignored_inputs = set()
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)


    # subworkflows
    def test_one_call_sub(self) -> None:
        wf = MinimalTaskInputsTestWF4()
        do_preprocessing_workflow(wf)
        tool = wf.step_nodes['stp1'].tool

        actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
        expected_task_inputs = {'inFile'}
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
        actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
        expected_param_inputs = {'inStr1', 'inInt2'}
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
        actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
        expected_static_inputs = {'inInt1'}
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        expected_ignored_inputs = {
            'inStr2',
            'inStr3',
            'inInt3',
        }
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)

    def test_two_calls_sub(self) -> None:
        wf = MinimalTaskInputsTestWF5()
        do_preprocessing_workflow(wf)

        # TaskInputsTestTool1
        tool = wf.step_nodes['stp1'].tool
        actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
        expected_task_inputs = {
            'inFile',
            'inStr2',
            'inStr3',
            'inInt1',
            'inInt2',
            'inInt4',
        }
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
        actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
        expected_param_inputs = {'inStr1'}
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
        actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
        expected_static_inputs = set()
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        expected_ignored_inputs = {
            'inStr4',
            'inInt3',
        }
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)
        
        # SubMinimalTaskInputsTestWF
        tool = wf.step_nodes['stp2'].tool

        actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
        expected_task_inputs = {'inFile'}
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
        actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
        expected_param_inputs = {'inStr1', 'inInt2'}
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
        actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
        expected_static_inputs = {'inInt1'}
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        expected_ignored_inputs = {
            'inStr2',
            'inStr3',
            'inInt3',
        }
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)



class TestPreprocessingHelpers(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_get_true_workflow_inputs_1(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)

        # main data source wf
        actual_ids = nextflow.preprocessing.params_channels.get_all_workflow_inputs(wf)
        expected_ids = set([
            'inFile1',
            'inFileOpt1',
            'inStr1',
            'inStrOpt1',
        ])
        self.assertSetEqual(actual_ids, expected_ids)

    def test_get_true_workflow_inputs_2(self) -> None:
        # assembly wf
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        actual_ids = nextflow.preprocessing.params_channels.get_all_workflow_inputs(wf)
        expected_ids = set([
            "inForwardReads",
            "inReverseReads",
            "inLongReads",
            "testInput",
            "fastqc1_adapters",
            "fastqc1_contaminants",
            "fastqc1_limits",
            "fastqc2_adapters",
            "fastqc2_contaminants",
            "fastqc2_limits",
        ])
        self.assertSetEqual(actual_ids, expected_ids)

    def test_get_linkable_params_sub_connected(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        # sub connected
        step = wf.step_nodes['stp4']
        linking_dict = nextflow.preprocessing.params_channels.get_linkable_params(step.tool, step.sources)
        actual_param_names = set([x.name for x in linking_dict.values()])
        expected_ids = set([
            'in_file1',
            'in_file_opt1',
            'in_str1',
            'in_str_opt1',
        ])
        self.assertSetEqual(actual_param_names, expected_ids)
        
    def test_get_linkable_params_sub_static(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        # sub static
        step = wf.step_nodes['stp5']
        linking_dict = nextflow.preprocessing.params_channels.get_linkable_params(step.tool, step.sources)
        actual_param_names = set([x.name for x in linking_dict.values()])
        expected_ids = set([
            'in_file1',
            'in_file_opt1',
        ])
        self.assertSetEqual(actual_param_names, expected_ids)
        
    def test_get_linkable_params_sub_omitted(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        # sub omitted
        step = wf.step_nodes['stp6']
        linking_dict = nextflow.preprocessing.params_channels.get_linkable_params(step.tool, step.sources)
        actual_param_names = set([x.name for x in linking_dict.values()])
        expected_ids = set([
            'in_file1',
        ])
        self.assertSetEqual(actual_param_names, expected_ids)

    def test_get_linkable_params_sub_connected_sub_connected(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        subwf = wf.step_nodes['stp4'].tool
        step = subwf.step_nodes['stp4']
        linking_dict = nextflow.preprocessing.params_channels.get_linkable_params(step.tool, step.sources)
        actual_param_names = set([x.name for x in linking_dict.values()])
        expected_ids = set([
            'in_file1',
            'in_file_opt1',
            'in_str1',
            'in_str_opt1',
        ])
        self.assertSetEqual(actual_param_names, expected_ids)
    
    def test_get_linkable_params_sub_connected_sub_static(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        subwf = wf.step_nodes['stp4'].tool
        step = subwf.step_nodes['stp5']
        linking_dict = nextflow.preprocessing.params_channels.get_linkable_params(step.tool, step.sources)
        actual_param_names = set([x.name for x in linking_dict.values()])
        expected_ids = set([
            'in_file1',
            'in_file_opt1',
        ])
        self.assertSetEqual(actual_param_names, expected_ids)
    
    def test_get_linkable_params_sub_connected_sub_omitted(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        subwf = wf.step_nodes['stp4'].tool
        step = subwf.step_nodes['stp6']
        linking_dict = nextflow.preprocessing.params_channels.get_linkable_params(step.tool, step.sources)
        actual_param_names = set([x.name for x in linking_dict.values()])
        expected_ids = set([
            'in_file1',
        ])
        self.assertSetEqual(actual_param_names, expected_ids)




class TestDataSourceVariables(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()
        self.wf = DataSourceTestWF()
        do_preprocessing_workflow(self.wf)
        self.pdsv = nextflow.data_sources.tvn_register.data_structure
        print(nextflow.data_sources.tvn_register.to_string())

    def test_stp1(self) -> None:
        label = 'main.stp1'
        self.assertEqual(self.pdsv[label]['inFile'].value, 'in_file')
        self.assertEqual(self.pdsv[label]['inFileOpt'].value, 'in_file_opt')
        self.assertEqual(self.pdsv[label]['inStr'].value, 'params.in_str1')
        self.assertEqual(self.pdsv[label]['inStrOpt'].value, 'params.in_str_opt1')
    
    def test_stp2(self) -> None:
        label = 'main.stp2'
        self.assertEqual(self.pdsv[label]['inFile'].value, 'in_file')
        self.assertEqual(self.pdsv[label]['inFileOpt'].value, 'in_file_opt')
        self.assertEqual(self.pdsv[label]['inStr'].value, None)
        self.assertEqual(self.pdsv[label]['inStrOpt'].value, None)
    
    def test_stp3(self) -> None:
        label = 'main.stp3'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': None,
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp1(self) -> None:
        label = 'main.stp4.stp1'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': 'in_file_opt',
            'inStr': 'params.in_str1',
            'inStrOpt': 'params.in_str_opt1',
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp2(self) -> None:
        label = 'main.stp4.stp2'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': 'in_file_opt',
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp3(self) -> None:
        label = 'main.stp4.stp3'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': None,
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp4_stp1(self) -> None:
        label = 'main.stp4.stp4.stp1'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': 'in_file_opt',
            'inStr': 'params.in_str1',
            'inStrOpt': 'params.in_str_opt1',
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp4_stp2(self) -> None:
        label = 'main.stp4.stp4.stp2'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': 'in_file_opt',
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp4_stp3(self) -> None:
        label = 'main.stp4.stp4.stp3'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': None,
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)

    def test_stp4_stp5_stp1(self) -> None:
        label = 'main.stp4.stp5.stp1'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': 'in_file_opt',
            'inStr': 'params.stp4.stp5.in_str3',
            'inStrOpt': 'params.stp4.stp5.in_str_opt3',
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp5_stp2(self) -> None:
        label = 'main.stp4.stp5.stp2'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': 'in_file_opt',
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp5_stp3(self) -> None:
        label = 'main.stp4.stp5.stp3'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': None,
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)

    def test_stp4_stp6_stp1(self) -> None:
        label = 'main.stp4.stp6.stp1'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': 'params.stp4.stp6.in_file_opt3',
            'inStr': 'params.stp4.stp6.in_str3',
            'inStrOpt': 'params.stp4.stp6.in_str_opt3',
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp6_stp2(self) -> None:
        label = 'main.stp4.stp6.stp2'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': 'params.stp4.stp6.in_file_opt3',
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)
    
    def test_stp4_stp6_stp3(self) -> None:
        label = 'main.stp4.stp6.stp3'
        expected_map = {
            'inFile': 'in_file',
            'inFileOpt': None,
            'inStr': None,
            'inStrOpt': None,
        }
        self.assertDictEqual(self.pdsv[label], expected_map)



class TestToGroovyStr(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_string(self) -> None:
        inp = 'Hello'
        expected = "'Hello'"
        actual = to_groovy(inp, String())
        self.assertEquals(expected, actual)
    
    def test_numeric(self) -> None:
        # int
        inp = 5
        expected = '5'
        actual = to_groovy(inp, Int())
        self.assertEquals(expected, actual)
        # float
        inp = 5.2
        expected = '5.2'
        actual = to_groovy(inp, Float())
        self.assertEquals(expected, actual)
    
    def test_bool(self) -> None:
        # true
        inp = True
        expected = 'true'
        actual = to_groovy(inp, Boolean())
        self.assertEquals(expected, actual)
        # false
        inp = False
        expected = 'false'
        actual = to_groovy(inp, Boolean())
        self.assertEquals(expected, actual)
    
    def test_none(self) -> None:
        inp = None
        expected = 'null'
        actual = to_groovy(inp, String())
        self.assertEquals(expected, actual)
    
    def test_empty_array(self) -> None:
        inp = []
        expected = "[]"
        actual = to_groovy(inp, String())
        self.assertEquals(expected, actual)
    
    def test_string_array(self) -> None:
        inp = ['Hello']
        expected = "['Hello']"
        actual = to_groovy(inp, String())
        self.assertEquals(expected, actual)
    
    def test_numeric_array(self) -> None:
        # int
        inp = [1, 2, 3]
        expected = '[1, 2, 3]'
        actual = to_groovy(inp, Int())
        self.assertEquals(expected, actual)
        # float
        inp = [1.1, 2.2, 3.3]
        expected = '[1.1, 2.2, 3.3]'
        actual = to_groovy(inp, Float())
        self.assertEquals(expected, actual)
    
    def test_bool_array(self) -> None:
        inp = [True, False]
        expected = '[true, false]'
        actual = to_groovy(inp)
        self.assertEquals(expected, actual)
    
    def test_none_array(self) -> None:
        inp = [None, None]
        expected = '[null, null]'
        actual = to_groovy(inp)
        self.assertEquals(expected, actual)




class TestSettings(unittest.TestCase):
    def setUp(self) -> None:
        reset_globals()

    def test_get_settings(self):
        self.assertEquals(settings.translate.nextflow.LIB_FILENAME, 'lib.nf')
        self.assertEquals(settings.translate.nextflow.CONFIG_FILENAME, 'nextflow.config')
    
    def test_set_settings(self):
        settings.translate.nextflow.MINIMAL_PROCESS = False
        self.assertEquals(settings.translate.nextflow.MINIMAL_PROCESS, False)




class TestParams(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()
    
    def test_channel_params(self) -> None:
        """
        Every channel requires a param. 
        Channels are created from a subset of workflow inputs.
        """
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        params = nextflow.params.getall()
        param_ids = {p.name for p in params}
        expected_ids = {
            'test_input', 
            'fastqc2_limits', 
            'fastqc1_adapters', 
            'fastqc2_adapters', 
            'fastqc2_contaminants', 
            'in_long_reads', 
            'in_reverse_reads', 
            'in_forward_reads', 
            'fastqc1_contaminants', 
            'fastqc1_limits'
        }
        for inp_id in expected_ids:
            self.assertIn(inp_id, param_ids)

    def test_fed_step_inputs(self) -> None:
        wf = StepInputsWFInputTestWF()
        
        # minimal
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'outdir': "'./outputs'",
            'in_file': 'null',
            'in_file_opt': "'NO_FILE1'",
            'in_str': 'null',
            'in_str_opt': 'null',
            'in_int': 'null',
            'in_bool': 'null',
        }
        self.assertEquals(actual_params, expected_params)
    
    def test_optional_workflow_inputs(self) -> None:
        wf = OptionalTestWF()
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'in_file': "null",
            'in_file_array': "[]",
            'in_file_array_opt': "['NO_FILE5']",
            'in_file_opt': "'NO_FILE6'",
            'in_secondary': "null",
            'in_secondary_array': "[]",
            'in_secondary_array_opt': "[['NO_FILE1', 'NO_FILE2']]",
            'in_secondary_opt': "['NO_FILE3', 'NO_FILE4']",
            'in_file_pair': "[]",
            'in_file_pair_opt': "[]",
            'outdir': "'./outputs'",
        }
        self.assertEquals(actual_params, expected_params)

    def test_static_step_inputs(self) -> None:
        wf = StepInputsStaticTestWF()
        
        # minimal
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'outdir': "'./outputs'",
            'in_file': 'null',
            # 'in_str': 'null',
            # 'in_int': 'null',
            # 'in_bool': 'null',
            # 'stp2_pos_default': '100',
            # 'stp2_pos_optional': "'static'",
            # 'stp2_flag_true': 'false',
            # 'stp2_flag_false': 'true',
            # 'stp2_opt_basic': "'static'",
            # 'stp2_opt_default': '100',
            # 'stp2_opt_optional': "''"
        }
        self.assertEquals(actual_params, expected_params)
    
    def test_static_and_omitted_step_inputs(self) -> None:
        wf = StepInputsPartialStaticTestWF()
        
        # minimal
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'outdir': "'./outputs'",
            'in_file': 'null',
            # 'in_str': 'null',
            # 'in_int': 'null',
            # 'in_bool': 'null',
            # 'stp3_opt_basic': "'static'",
            # 'stp3_opt_default': '100',
            # 'stp3_opt_optional': "''"
        }
        self.assertEquals(actual_params, expected_params)
    
    def test_omitted_step_inputs(self) -> None:
        wf = StepInputsMinimalTestWF()
        
        # minimal
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'outdir': "'./outputs'",
            'in_file': 'null',
            'in_str': 'null',
            # 'in_int': 'null',
            # 'in_bool': 'null',
        }
        self.assertEquals(actual_params, expected_params)

    def test_array_inputs(self) -> None:
        wf = ArrayIOTestWF()
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'outdir': "'./outputs'",
            'in_file_array': '[]',
            'in_str_array': '[]',
            'in_int_array': '[]',
        }
        self.assertEquals(actual_params, expected_params)

    def test_secondaries_inputs(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'in_alignments': 'null',
        }
        self.assertDictContainsSubset(expected_params, actual_params)
    
    def test_secondaries_array_inputs(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'in_alignments_arr': '[]',
        }
        self.assertDictContainsSubset(expected_params, actual_params)
    
    def test_subworkflow(self) -> None:
        wf = Subworkflow2TestWF()
        do_preprocessing_workflow(wf)
        actual_params = nextflow.params.serialize()
        expected_params = {
            'in_file': 'null',
            'new_filename': 'null',
        }
        self.assertDictContainsSubset(expected_params, actual_params)

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError



class TestFiles(unittest.TestCase):

    def setUp(self) -> None:
        self.maxDiff = None
        reset_globals()
        wf = DuplicateTasksTestWF()
        do_preprocessing_workflow(wf)
        self.maintask, self.subtask_dict = translator.translate_workflow_internal(wf)

    def test_duplicate_tool_usage(self) -> None:
        actual_subtasks = len(self.subtask_dict)
        expected_subtasks = 12
        self.assertEqual(actual_subtasks, expected_subtasks)

        # checking values supplied 
        subtask = self.subtask_dict['modules/stp1']
        print(subtask)
        raise NotImplementedError

    def test_duplicate_subworkflow_usage(self) -> None:
        raise NotImplementedError


class TestFileFormats(unittest.TestCase):

    def setUp(self) -> None:
        self.maxDiff = None
        reset_globals()

    def test_main_workflow(self) -> None:
        wf = AssemblyTestWF()
        mainstr, _ = translator.translate_workflow_internal(wf)
        expected_lines = [
            "nextflow.enable.dsl=2",
            "include { FASTQC as FASTQC1 } from './modules/fastqc'",
            "include { FASTQC as FASTQC2 } from './modules/fastqc'",
            "include { FASTQC as FASTQC3 } from './modules/fastqc'",
            "include { CAT_TEST_TOOL } from './modules/cat_test_tool'",
            "include { UNICYCLER } from './modules/unicycler'",
            "ch_in_forward_reads     = Channel.fromPath( params.in_forward_reads )",
            "ch_in_long_reads        = Channel.fromPath( params.in_long_reads )",
            "ch_in_reverse_reads     = Channel.fromPath( params.in_reverse_reads )",
            "ch_test_input           = Channel.fromPath( params.test_input )",
            "ch_fastqc1_adapters     = Channel.fromPath( params.fastqc1_adapters )",
            "ch_fastqc1_contaminants = Channel.fromPath( params.fastqc1_contaminants )",
            "ch_fastqc1_limits       = Channel.fromPath( params.fastqc1_limits )",
            "ch_fastqc2_adapters     = Channel.fromPath( params.fastqc2_adapters )",
            "ch_fastqc2_contaminants = Channel.fromPath( params.fastqc2_contaminants )",
            "ch_fastqc2_limits       = Channel.fromPath( params.fastqc2_limits )",
            "workflow  {",
            "FASTQC1(",
            "ch_in_forward_reads,",
            "ch_fastqc1_adapters,",
            "ch_fastqc1_contaminants,",
            "ch_fastqc1_limits",
            ")",
            "FASTQC2(",
            "ch_in_reverse_reads,",
            "ch_fastqc2_adapters,",
            "ch_fastqc2_contaminants,",
            "ch_fastqc2_limits",
            ")",
            "FASTQC3(",
            "ch_test_input,",
            "null,",
            "null,",
            "null",
            ")",
            "CAT_TEST_TOOL(",
            "FASTQC3.out.outTextFile",
            ")",
            "UNICYCLER(",
            "ch_in_forward_reads,",
            "ch_in_reverse_reads,",
            "ch_in_long_reads",
            ")",
            "}",
        ]
        actual_lines = mainstr.split('\n')
        actual_lines = [ln.strip() for ln in actual_lines]
        actual_lines = [ln for ln in actual_lines if ln != '']
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

    def test_process(self) -> None:
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        _, substr_dict = translator.translate_workflow_internal(wf)
        expected_lines = [
            'nextflow.enable.dsl=2',
            'process FASTQC {',
            'debug true',
            'container "quay.io/biocontainers/fastqc:0.11.8--2"',
            'publishDir "${params.outdir}/fastqc"',
            'input:',
            'path input_file, stageAs: \'input_file\'',
            'path adapters, stageAs: \'adapters\'',
            'path contaminants, stageAs: \'contaminants\'',
            'path limits, stageAs: \'limits\'',
            'output:',
            'path "output.html", emit: outHtmlFile',
            'path "output.txt", emit: outTextFile',
            'script:',
            'def adapters = adapters != \'NO_FILE7\' ? "--adapters ${adapters}" : ""',
            'def contaminants = contaminants != \'NO_FILE8\' ? "--contaminants ${contaminants}" : ""',
            'def limits = limits != \'NO_FILE9\' ? "--limits ${limits}" : ""',
            '"""',
            'fastqc \\',
            '${adapters} \\',
            '${contaminants} \\',
            '${limits} \\',
            '--kmers 7 \\',
            '${input_file} \\',
            '"""',
            '}',
        ]
        process_str = substr_dict['modules/fastqc']
        actual_lines = process_str.split('\n')
        actual_lines = [ln.strip() for ln in actual_lines]
        actual_lines = [ln for ln in actual_lines if ln != '']
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_subworkflow(self) -> None:
        wf = SubworkflowTestWF()
        _, substr_dict = translator.translate_workflow_internal(wf)
        expected_lines = [
            "nextflow.enable.dsl=2",
            "include { STRING_TOOL } from '../modules/string_tool'",
            "include { STRING_OPT_TOOL } from '../modules/string_opt_tool'",
            "include { ORANGES_SUBWORKFLOW } from './oranges_subworkflow'",
            "workflow APPLES_SUBWORKFLOW {",
            "take:",
            "main:",
            "STRING_TOOL()",
            "STRING_OPT_TOOL()",
            "ORANGES_SUBWORKFLOW(",
            "STRING_TOOL.out.out",
            ")",
            "emit:",
            "outStringFile = STRING_TOOL.out.out",
            "outIntFile = ORANGES_SUBWORKFLOW.out.out",
            "}",
        ]
        subwf_str = substr_dict['subworkflows/apples_subworkflow']
        actual_lines = subwf_str.split('\n')
        actual_lines = [ln.strip() for ln in actual_lines]
        actual_lines = [ln for ln in actual_lines if ln != '']
        self.assertEqual(actual_lines, expected_lines)
    
    def test_config(self) -> None:
        # TODO expand this to subworkflow with subworkflow, process specific params
        wf = SubworkflowTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        expected_lines = [
            "docker.enabled = true",
            "params {",
            "// OUTPUT DIRECTORY",
            "outdir  = './outputs'",
            "// INPUTS",
            "in_file     = null",
            "in_int      = null",
            "in_str      = null",
            "in_str_opt  = null",
            "}",
        ]
        actual_lines = config.split('\n')
        actual_lines = [ln.strip() for ln in actual_lines]
        actual_lines = [ln for ln in actual_lines if ln != '']
        self.assertEqual(actual_lines, expected_lines)

class TestChannels(unittest.TestCase):
    
    def setUp(self) -> None:
        self.maxDiff = None
        reset_globals()

    def test_infer_wf_inputs(self) -> None:
        # checks the inferred wf inputs (from total wf inputs) are correct
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        channel_ids = {c.name for c in nextflow.channels.getall(scope)}
        expected_ids = {
            'ch_in_forward_reads',
            'ch_in_reverse_reads',
            'ch_in_long_reads',
            'ch_test_input',
            'ch_fastqc1_adapters',
            'ch_fastqc1_contaminants',
            'ch_fastqc1_limits',
            'ch_fastqc2_adapters',
            'ch_fastqc2_contaminants',
            'ch_fastqc2_limits',
        }
        self.assertEqual(channel_ids, expected_ids)
    
    @unittest.skip('deprecated optionality handling')
    def test_optional_file_channels(self) -> None:
        """
        Every Optional(File) type wf input should have a channel. 
        '.ifEmpty(null)' should appear in the channel string definition.
        """
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        all_channels = nextflow.channels.getall(scope)
        expected_channel_names = {
            'ch_fastqc1_adapters',
            'ch_fastqc1_contaminants',
            'ch_fastqc1_limits',
            'ch_fastqc2_adapters',
            'ch_fastqc2_contaminants',
            'ch_fastqc2_limits',
        }
        optional_channels = [ch for ch in all_channels if ch.name in expected_channel_names]

        # check each expected channel exists
        self.assertEqual(len(optional_channels), len(expected_channel_names))
        
        # for each optional channel, check it has correct format
        for channel in optional_channels:
            self.assertIn('.ifEmpty( null )', channel.get_string())

    def test_nonfile_no_channel(self) -> None:
        """
        Non-File-type wf input should not have nextflow.channels.
        """
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        nonfile_wf_input_ids = [
            'unicycler_kmers',
            'unicycler_scores',
            'unicycler_startGeneCov',
            'unicycler_startGeneId',
        ]
        for inp_id in nonfile_wf_input_ids:
            node = wf.input_nodes[inp_id]
            self.assertFalse(nextflow.channels.exists(scope, node.uuid))

    def test_array_inputs(self) -> None:
        wf = ArrayIOTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        
        # check expected channels are created
        channels_ids = {c.name for c in nextflow.channels.getall(scope)}
        expected_ids = {
            'ch_in_file_array',
        }
        self.assertEqual(channels_ids, expected_ids)

        # check channels can be looked up using janis_uuid
        inp = wf.input_nodes['inFileArray']
        inp_ch = nextflow.channels.get(scope, inp.uuid)
        self.assertIsNotNone(inp_ch)

        # check channel definition is correct
        actual_definition = inp_ch.definition
        expected_definition = 'ch_in_file_array = Channel.fromPath( params.in_file_array ).toList()'
        self.assertEqual(actual_definition, expected_definition)
        
    def test_secondaries_inputs(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        
        # check expected channels are created
        channels_ids = {c.name for c in nextflow.channels.getall(scope)}
        expected_ids = {
            'ch_in_alignments',
        }
        for expected_id in expected_ids:
            self.assertIn(expected_id, channels_ids)

        # check channels can be looked up using janis_uuid
        inp = wf.input_nodes['inAlignments']
        inp_ch = nextflow.channels.get(scope, inp.uuid)
        self.assertIsNotNone(inp_ch)

        # check channel definition is correct
        actual_definition = inp_ch.definition
        expected_definition = 'ch_in_alignments = Channel.fromPath( params.in_alignments ).toList()'
        self.assertEqual(actual_definition, expected_definition)
    
    def test_secondaries_array_inputs(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()

        # check expected channels are created
        channels_ids = {c.name for c in nextflow.channels.getall(scope)}
        expected_ids = {
            'ch_in_alignments_arr',
        }
        for expected_id in expected_ids:
            self.assertIn(expected_id, channels_ids)

        # check channels can be looked up using janis_uuid
        inp = wf.input_nodes['inAlignmentsArr']
        inp_ch = nextflow.channels.get(scope, inp.uuid)
        self.assertIsNotNone(inp_ch)

        # check channel definition is correct
        actual_definition = inp_ch.definition
        expected_definition = 'ch_in_alignments_arr = Channel.fromPath( params.in_alignments_arr.flatten() ).collate( 2 )'
        self.assertEqual(actual_definition, expected_definition)

    def test_filename_types(self) -> None:
        wf = FilenameTestWF1()
        do_preprocessing_workflow(wf)
        print(nextflow.data_sources.tvn_register.to_string())
        scope = nextflow.Scope()
        channels_ids = {c.name for c in nextflow.channels.getall(scope)}
        expected_ids = {
            'ch_in_file',
            'ch_in_file_opt',
        }
        self.assertEqual(channels_ids, expected_ids)
    
    


class TestChannelsSubworkflows(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()
        self.wf = DataSourceTestWF()
        do_preprocessing_workflow(self.wf)
        self.scope = nextflow.Scope()

    def test_subworkflows_main(self) -> None:
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt1 = Channel.fromPath( params.in_file_opt1 )',
            'ch_in_file1     = Channel.fromPath( params.in_file1 )',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
    def test_subworkflows_sub_connected(self) -> None:
        step =self.wf.step_nodes['stp4']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file2',
            'ch_in_file_opt2',
            'ch_in_str2',
            'ch_in_str_opt2',
        ])
        self.assertEqual(channel_declarations, expected_declarations)

    def test_subworkflows_sub_static(self) -> None:
        step =self.wf.step_nodes['stp5']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt2',
            'ch_in_file2',
            'ch_in_str2',
            'ch_in_str_opt2',
        ])
        self.assertEqual(channel_declarations, expected_declarations)

    def test_subworkflows_sub_omitted(self) -> None:
        step =self.wf.step_nodes['stp6']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt2',
            'ch_in_file2',
            'ch_in_str2',
            'ch_in_str_opt2',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
    def test_subworkflows_sub_connected_sub_connected(self) -> None:
        step =self.wf.step_nodes['stp4']
        self.scope.update(step)
        step = step.tool.step_nodes['stp4']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
    def test_subworkflows_sub_connected_sub_static(self) -> None:
        step =self.wf.step_nodes['stp4']
        self.scope.update(step)
        step = step.tool.step_nodes['stp5']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
    def test_subworkflows_sub_connected_sub_omitted(self) -> None:
        step =self.wf.step_nodes['stp4']
        self.scope.update(step)
        step = step.tool.step_nodes['stp6']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
 
    def test_subworkflows_sub_static_sub_connected(self) -> None:
        step =self.wf.step_nodes['stp5']
        self.scope.update(step)
        step = step.tool.step_nodes['stp4']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
    def test_subworkflows_sub_static_sub_static(self) -> None:
        step =self.wf.step_nodes['stp5']
        self.scope.update(step)
        step = step.tool.step_nodes['stp5']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
    def test_subworkflows_sub_static_sub_omitted(self) -> None:
        step =self.wf.step_nodes['stp5']
        self.scope.update(step)
        step = step.tool.step_nodes['stp6']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)

    def test_subworkflows_sub_omitted_sub_connected(self) -> None:
        step =self.wf.step_nodes['stp6']
        self.scope.update(step)
        step = step.tool.step_nodes['stp4']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
    def test_subworkflows_sub_omitted_sub_static(self) -> None:
        step =self.wf.step_nodes['stp6']
        self.scope.update(step)
        step = step.tool.step_nodes['stp5']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
    def test_subworkflows_sub_omitted_sub_omitted(self) -> None:
        step =self.wf.step_nodes['stp6']
        self.scope.update(step)
        step = step.tool.step_nodes['stp6']
        self.scope.update(step)
        
        # main wf channels
        channel_declarations = nextflow.channels.getstr(self.scope)
        channel_declarations = channel_declarations.strip('[]').split('\n')
        channel_declarations = {x for x in channel_declarations if x != ''}
        expected_declarations = set([
            'ch_in_file_opt3',
            'ch_in_file3',
            'ch_in_str3',
            'ch_in_str_opt3',
        ])
        self.assertEqual(channel_declarations, expected_declarations)
    
        
    @unittest.skip('not implemented')
    def test_channel_methods(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError


class TestCmdtoolProcess(unittest.TestCase):
    """
    Tests janis CommandTool can be parsed to nextflow process (end-to-end).
    """

    def setUp(self) -> None:
        reset_globals() 
    
    def test_fastqc(self) -> None:
        tool = FastqcTestTool()
        process = translate(tool, 'nextflow')
        print(process)
        print()
    
    def test_bwamem(self) -> None:
        tool = BwaMemTestTool()
        process = translate(tool, 'nextflow')
        print(process)
        print()
    
    def test_gridss(self) -> None:
        tool = GridssTestTool()
        process = translate(tool, 'nextflow')
        print(process)
        print()
 

class TestCmdtoolProcessDirectives(unittest.TestCase):
    """
    Tests identifying tool inputs which should be process inputs.
    Need a process input for each tool input in step sources

    INCLUDES WORKFLOW OUTPUTS as publishDir
    """

    def setUp(self) -> None:
        reset_globals()

    def test_directives_order(self) -> None:
        wf = DirectivesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        process = translator.handle_container(step.tool, process)
        directives = nextflow.ordering.order_nf_directives(process.directives)
        actual_order = [type(x).__name__ for x in directives]
        expected_order = [
            'NFDebugDirective',
            'NFContainerDirective',
            'NFPublishDirDirective',
            'NFCpusDirective',
            'NFDiskDirective',
            'NFMemoryDirective',
            'NFTimeDirective',
        ]
        for actual, expected in zip(actual_order, expected_order):
            self.assertEqual(actual, expected)

    def test_directives(self) -> None:
        wf = DirectivesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        process = translator.handle_container(step.tool, process)
        actual_directives = {d.get_string() for d in process.directives}
        expected_directives = {
            'container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"',
            'publishDir "${params.outdir}/resources_test_tool"',
            'debug true',
            'disk "${params.resources_test_tool.disk}"',
            'memory "${params.resources_test_tool.memory}"',
            'time "${params.resources_test_tool.time}"'
        }
        for direc in expected_directives:
            self.assertIn(direc, actual_directives)
    
    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError
    



class TestCmdtoolProcessInputs(unittest.TestCase):
    """
    Tests identifying tool inputs which should be process inputs.
    Need a process input for each tool input in step sources
    """

    def setUp(self) -> None:
        reset_globals()
        
    def test_stage_as(self) -> None:
        wf = ProcessInputsTestWF()
        do_preprocessing_workflow(wf)

        step = wf.step_nodes["stp2"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path fastq_inp1, stageAs: \'fastq_inp1.fastq\'',
            'path fastq_inp2, stageAs: \'fastq_inp2.fastq\'',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)

    def test_basic(self) -> None:
        # DO NOT need a process input for each static value in step sources.
        # non-files are fed data via params. 
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["unicycler"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            "path option1, stageAs: 'option1'",
            "path option2, stageAs: 'option2'",
            "path option_l, stageAs: 'option_l'",
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    def test_file_pairs_fmt(self) -> None:
        wf = FilePairsTestWF()
        do_preprocessing_workflow(wf)

        # filepair
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path reads'
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
        
        # filepair array
        step = wf.step_nodes["stp3"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        expected_inputs = {
            'path reads'
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)

    def test_arrays_fmt(self) -> None:
        # definition should be the same as singles. 
        # nextflow doesn't differentiate. 
        wf = ArrayStepInputsTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path pos_basic',
            'path pos_basic2',
            'val pos_default',
            'val pos_optional',
            'val opt_basic',
            'val opt_default',
            'val opt_optional',
        }
        self.assertEqual(actual_inputs, expected_inputs)

    def test_secondaries_fmt(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'tuple path(bam), path(bai)'
        }
        self.assertEqual(actual_inputs, expected_inputs)   
    
    def test_secondaries_array_fmt(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp4"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path indexed_bam_flat',
        }
        self.assertEqual(actual_inputs, expected_inputs)  

    def test_filename_types(self) -> None:
        wf = FilenameTestWF1()
        do_preprocessing_workflow(wf)

        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            "path inp1, stageAs: 'inp1'",
            "val inp2",
        }
        self.assertEqual(actual_inputs, expected_inputs)  

        step = wf.step_nodes["stp2"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            "path inp1, stageAs: 'inp1'",
        }
        self.assertEqual(actual_inputs, expected_inputs)  

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError
    




class TestCmdtoolProcessOutputs(unittest.TestCase):
    """
    Need a process output for each tool output.
    """

    def setUp(self) -> None:
        reset_globals()
        self.wf = OutputCollectionTestWF()
        do_preprocessing_workflow(self.wf)

    def test_get_fmttype(self) -> None:
        pass
        # get_fmttype

    def test_stdout(self):
        wf = BasicIOTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'stdout, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)  

    def test_wildcard(self) -> None:
        step = self.wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "myfile.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_wildcard_array(self) -> None:
        wf = WildcardSelectorOutputTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp2"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "*.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_input_selector(self) -> None:
        wf = InputSelectorTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path inp, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_input_selector_param(self) -> None:
        step = self.wf.step_nodes["stp4"]
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "myfile.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
    def test_input_selector_array(self) -> None:
        wf = InputSelectorTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp3']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path inp, emit: out'}
        print(process.get_string())
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_filenames_generated(self) -> None:
        wf = FilenameTestWF2()
        do_preprocessing_workflow(wf)
        
        # inp1 is in step.sources
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {
            'path "${inp1.simpleName + ".csv"}", emit: out3',
            'path "${"generated" + ".csv"}", emit: out4',
            'path "generated.csv", emit: out5',
            'path "generated.merged.csv", emit: out6'
        }
        print(process.get_string())
        self.assertEqual(actual_outputs, expected_outputs)
        
    def test_filenames_referenced(self) -> None:
        # inp1, inp1_1, inp4, inp5, inp6 are in step.sources
        wf = FilenameTestWF1()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp5']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {
            'path "${inp1.simpleName + ".csv"}", emit: out3',
            'path "${inp4 + ".csv"}", emit: out4',
            'path "${inp5 + ".csv"}", emit: out5',
            'path "${inp6 + ".merged.csv"}", emit: out6'
        }
        print(process.get_string())
        self.assertEqual(actual_outputs, expected_outputs)

    def test_file_pair(self) -> None:
        # eg read1.fastq, read2.fastq
        # collection method is list, len(list) == 2.
        step = self.wf.step_nodes['stp6']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "[${inp.simpleName + "-R1.fastq"}, ${inp.simpleName + "-R2.fastq"}]", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_secondaries(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'tuple path("*.bam"), path("*.bam.bai"), emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_secondaries_replaced(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp3']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'tuple path("*.bam"), path("*.bai"), emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_secondaries_edge_basename(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp5']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        print(process.outputs[0].get_string())
        expected_outputs = {
            'tuple path("${bam.name}"), path("*.bai"), emit: out'
        }
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_secondaries_edge_no_secondaries_present_as(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp6']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        print(process.outputs[0].get_string())
        expected_outputs = {
            'tuple path(inp), path("${inp}.tbi"), emit: out'
        }
        self.assertEqual(actual_outputs, expected_outputs)
    
    @unittest.skip('not implemented')
    def test_secondaries_array(self) -> None:
        # TODO make this an unsupported feature in release version
        # highly unlikely workflow would do this
        raise NotImplementedError
        
    def test_complex_expression(self) -> None:
        # two_value operator etc. uses ${} syntax around whole phrase.
        # strings inside are quoted. 
        step = self.wf.step_nodes['stp5']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "${inp.simpleName + ".gz"}", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_edge_markduplicates_metrics(self) -> None:
        step = self.wf.step_nodes['stp8']
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_outputs = {out.get_string() for out in process.outputs}
        print(process.get_string())
        expected_outputs = {
            'path "${["hello", "generated"].find{ it != null } + ".metrics.txt"}", emit: metrics'
        }
        self.assertEqual(actual_outputs, expected_outputs)

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError



class TestCmdtoolProcessScript(unittest.TestCase):
    """
    Tests of the process prescript and script sections.
    """

    def setUp(self) -> None:
        reset_globals()

    def test_multiple_statements(self) -> None:
        pass
    
    def test_directories_to_create(self) -> None:
        pass
    
    def test_files_to_create_cmdtool(self) -> None:
        pass
    
    def test_files_to_create_codetool(self) -> None:
        pass
    
    def test_files_to_create_cmdtool_exprtool(self) -> None:
        pass

    def test_variables_defined(self) -> None:
        wf = EntityTraceTestWF()
        do_preprocessing_workflow(wf)

        # inputs referencing undefined inputs
        step = wf.step_nodes["stp8"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_prescript = process.pre_script
        assert(actual_prescript)
        expected_lines = {
            'def java_options = null',
        }

        for ln in expected_lines:
            self.assertIn(ln, actual_prescript)
        
        # arguments referencing undefined inputs
        step = wf.step_nodes["stp9"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        
        actual_prescript = process.pre_script
        assert(actual_prescript)
        expected_lines = {
            'def compression_level = null',
            'def java_options = null',
        }

        for ln in expected_lines:
            self.assertIn(ln, actual_prescript)
        
        # outputs referencing undefined inputs
        step = wf.step_nodes["stp10"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        
        actual_prescript = process.pre_script
        assert(actual_prescript)
        expected_lines = {
            'def java_options = null',
        }
        for ln in expected_lines:
            self.assertIn(ln, actual_prescript)

    def test_components_prescript(self) -> None:
        wf = StepInputsTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp1']
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_prescript = process.pre_script
        expected_lines = {
            'def pos_default = pos_default ? pos_default : 95',
            'def pos_optional = pos_optional ? pos_optional : ""',
            'def flag_true = flag_true == false ? "" : "--flag-true"',
            'def flag_false = flag_false ? "--flag-false" : ""',
            'def opt_default = opt_default ? opt_default : 5',
            'def opt_optional = opt_optional ? "--opt-optional ${opt_optional}" : ""',
        }
        for ln in expected_lines:
            self.assertIn(ln, actual_prescript)
    
    def test_components_script(self) -> None:
        wf = StepInputsTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_script = process.script
        expected_lines = [
            'echo',
            '${pos_basic}',
            '${pos_default}',
            '${pos_optional}',
            '${flag_true}',
            '${flag_false}',
            '--opt-basic ${opt_basic}',
            '--opt-default ${opt_default}',
            '${opt_optional}',
        ]
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

    def test_components_array_prescript(self) -> None:
        wf = ArrayStepInputsTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_prescript = set(process.pre_script.split('\n'))
        expected_prescript = set([
            'def pos_basic_joined = pos_basic.join(\' \')',
            'def pos_basic2_joined = pos_basic2 != [\'NO_FILE\'] ? pos_basic2.join(\' \') : ""',
            'def pos_default_joined = pos_default ? pos_default.join(\' \') : "1 2 3"',
            'def pos_optional_joined = pos_optional ? pos_optional.join(\' \') : ""',
            'def opt_basic_joined = opt_basic.join(\' \')',
            'def opt_default_items = opt_default ? opt_default.collect{ "--opt-default " + it }.join(\' \') : "--opt-default 1 --opt-default 2 --opt-default 3"',
            'def opt_optional_joined = opt_optional ? "--opt-optional " + opt_optional.join(\',\') : ""',
        ])
        actual_prescript = sorted(actual_prescript)
        expected_prescript = sorted(expected_prescript)
        for ln in actual_prescript:
            self.assertIn(ln, expected_prescript)
    
    def test_components_array_script(self) -> None:
        wf = ArrayStepInputsTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        actual_script = process.script
        expected_lines = [
            'echo',
            '${pos_basic_joined}',
            '${pos_basic2_joined}',
            '${pos_default_joined}',
            '${pos_optional_joined}',
            '--opt-basic=${opt_basic_joined}',
            '${opt_default_items}',
            '${opt_optional_joined}',
        ]
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)
    
    def test_secondaries(self) -> None:
        # name accession should be different?
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        
        # pre-script
        actual_pre_script = process.pre_script
        expected_pre_script = {}
        for ln in expected_pre_script:
            self.assertIn(ln, actual_pre_script)
        
        # script
        actual_script = process.script
        print(actual_script)
        expected_script = {
            'echo',
            '${bam}',
            '--inp ${bam}',
            '--inp-index-0 ${bam}',
            '--inp-index-1 ${bai}',
        }
        for ln in expected_script:
            self.assertIn(ln, actual_script)
    
    def test_secondaries_array(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp4"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        
        # pre-script
        actual_pre_script = process.pre_script
        expected_pre_script = {
            "def inp = get_primary_files(indexed_bam_flat)",
            "def inp_joined = inp.join(' ')"
        }
        for ln in expected_pre_script:
            self.assertIn(ln, actual_pre_script)
        
        # script
        actual_script = process.script
        expected_script = {
            'echo',
            '${inp_joined}',
            '--inp ${inp_joined}',
            '--inp-index-0 ${inp_joined[0]}',
            '--inp-index-1 ${inp_joined[1]}',
        }
        print(process.get_string())
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_filename_generated_tool(self):
        wf = FilenameTestWF1()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp3"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(f'actual: \n{process.get_string()}')
        expected = """\
process FILENAME_GENERATED_TOOL {
    debug true
    publishDir "${params.outdir}/filename_generated_tool"

    input:
    path file_inp, stageAs: 'file_inp.txt'
    path file_inp_optional, stageAs: 'file_inp_optional.txt'

    output:
    val "*", emit: out

    script:
    \"\"\"
    echo \\
    ${params.filename_generated_tool.inp} \\
    ${params.filename_generated_tool.inp_optional} \\
    ${file_inp.simpleName}.transformed.fnp \\
    ${file_inp_optional.simpleName}.optional.txt \\
    \"\"\"

}
"""     
        print(f'expected: \n{expected}')
        self.assertEqual(expected, process.get_string())
   
    def test_filename_types(self) -> None:
        wf = FilenameTestWF1()
        do_preprocessing_workflow(wf)

        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_script = process.script
        expected_lines = {
            'echo',
            '${inp1}',
            '${inp2}',
        }
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

        step = wf.step_nodes["stp2"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_script = process.script
        expected_lines = {
            'echo',
            '${inp1}',
            '${inp1.simpleName}.processed.txt',
        }
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

    def test_optional_types(self) -> None:
        wf = OptionalTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        print(process.get_string())
        actual_prescript = process.pre_script
        expected_lines = {
            "def in_secondary_array_opt = get_primary_files(indexed_bam_flat)", 
            "def in_secondary_array_opt_joined = in_secondary_array_opt != ['NO_FILE'] ? in_secondary_array_opt.join(' ') : \"\"", 
            "def in_secondary_opt = bam != 'NO_FILE' ? bam : \"\"", 
            "def in_file_pair_opt_joined = in_file_pair_opt != ['NO_FILE', 'NO_FILE'] ? in_file_pair_opt.join(' ') : \"\"", 
            "def in_file_array_opt_joined = in_file_array_opt != ['NO_FILE'] ? in_file_array_opt.join(' ') : \"\"", 
            "def in_file_opt = in_file_opt != 'NO_FILE' ? in_file_opt : \"\"", 
        }
        for ln in expected_lines:
            self.assertIn(ln, actual_prescript)

    @unittest.skip('not implemented')
    def test_translate_commandtool(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_translate_pythontool(self) -> None:
        raise NotImplementedError






class TestWorkflows(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()

    def test_basic(self) -> None:
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        process_dict = nextflow.generate.process.generate_processes(wf)
        workflow_dict = nextflow.generate.main.generate_workflows(wf, process_dict)
        nf_workflow = workflow_dict['UnicyclerAssembly']
        print(nf_workflow.get_string())
        actual_lines = split_to_lines(nf_workflow.get_string())
        expected_lines = [
            'workflow {',
            'FASTQC1(',
            'ch_in_forward_reads,',
            'fastqc1_adapters,',
            'fastqc1_contaminants,',
            'fastqc1_limits',
            ')',
            'FASTQC2(',
            'ch_in_reverse_reads,',
            'fastqc2_adapters,',
            'fastqc2_contaminants,',
            'fastqc2_limits',
            ')',
            'FASTQC3(',
            'ch_test_input,',
            'null,',
            'null,',
            'null',
            ')',
            'CAT_TEST_TOOL(',
            'FASTQC3.out.outTextFile',
            ')',
            'UNICYCLER(',
            'ch_in_forward_reads,',
            'ch_in_reverse_reads,',
            'ch_in_long_reads',
            ')',
            '}',
        ]
        self.assertListEqual(actual_lines, expected_lines)





class TestSubWorkflows(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()







class TestPythontoolProcessInputs(unittest.TestCase):
    
    def setUp(self) -> None:
        self.wf = InputsPythonToolTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    def test_input_generation(self) -> None:
        # File, String, Int input types
        step = self.wf.step_nodes["stp0"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path code_file',
            "path inp1, stageAs: 'inp1'"
        }
        self.assertEqual(actual_inputs, expected_inputs)
        
        # Array(String) input type
        step = self.wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path code_file',
        }
        self.assertEqual(actual_inputs, expected_inputs)
        
        # File (secondaries) input type
        step = self.wf.step_nodes["stp2"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            'path code_file',
            'tuple path(bam), path(bai)',
        }
        self.assertEqual(actual_inputs, expected_inputs)



class TestPythontoolProcessOutputs(unittest.TestCase):

    def setUp(self) -> None:
        self.wf = OutputsPythonToolTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()
    
    def test_output_generation(self) -> None:
        # file output
        step = self.wf.step_nodes['stp0']
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'val "${file("${task.workDir}/" + file("${task.workDir}/out_out").text.replace(\'"\', \'\'))}", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
        # String output
        step = self.wf.step_nodes['stp1']
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'val "${file("${task.workDir}/out_out").text}", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
        # Array(String) output
        step = self.wf.step_nodes['stp2']
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {
            'val "${file("${task.workDir}/out_out").text.replace(\'[\', \'\').replace(\']\', \'\')}", emit: out'
        }
        self.assertEqual(actual_outputs, expected_outputs)



class TestPythontoolProcess(unittest.TestCase):

    def setUp(self) -> None:
        self.wf = InputsPythonToolTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    def test_format(self) -> None:
        step = self.wf.step_nodes["stp0"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_lines = process.get_string().split('\n')
        actual_lines = [x.strip() for x in actual_lines]
        actual_lines = [x for x in actual_lines if x != '\n' and x != '']
        print(actual_lines)
        expected_lines = [
            'process STP0 {',
            'debug true',
            'publishDir "${params.outdir}/stp0"',
            'input:',
            'path code_file',
            'path inp1, stageAs: \'inp1\'',
            'output:',
            'val "${file("${task.workDir}/" + file("${task.workDir}/out_out").text.replace(\'"\', \'\'))}", emit: out',
            'exec:',
            'script:',
            '"""',
            '#!/usr/bin/env python',
            'from ${code_file.simpleName} import code_block',
            'import os',
            'import json',
            'result = code_block(inp1="${inp1}", inp2="${params.in_str}", inp3=${params.in_int})',
            'work_dir = os.getcwd()',
            'for key in result:',
            '    with open(os.path.join(work_dir, f"out_{key}"), "w") as fp:',
            '        fp.write(json.dumps(result[key]))',
            '"""',
            '}',
        ]
        expected_lines = [x.strip() for x in expected_lines]
        self.assertEqual(actual_lines, expected_lines)

class TestPythontoolProcessScript(unittest.TestCase):
    
    def setUp(self) -> None:
        self.wf = InputsPythonToolTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    def test_input_references(self) -> None:
        # File, String, Int input types
        step = self.wf.step_nodes["stp0"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_script = process.script
        expected_lines = {
            'result = code_block(inp1="${inp1}", inp2="${params.in_str}", inp3=${params.in_int})',
        }
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

        # Array(String) input type
        step = self.wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_script = process.script
        expected_lines = {
            'result = code_block(inp="${params.in_str_arr}".split(" "))',
        }
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

        # File (secondaries) input type
        step = self.wf.step_nodes["stp2"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_codetool(step.tool, step.sources, scope)
        actual_script = process.script
        expected_lines = {
            'result = code_block(inp="${bam}")',
        }
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)




class TestEntityTracing(unittest.TestCase):
    """
    This test group checks we can handle array / single datatype mismatches
    in a workflow. This occurs in wgsgermline. 
    """
    def setUp(self) -> None:
        self.maxDiff = None
        self.wf = EntityTraceTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    def test_trace_entity_counts4(self) -> None:
        step_id = 'stp4'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_counts = nextflow.trace.trace_entity_counts(src)
        print(actual_counts)
        expected_counts = {
            'StepTagInput': 1, 
            'Edge': 1, 
            'FirstOperator': 1, 
            'list': 1, 
            'InputNodeSelector': 1, 
            'InputNode': 1, 
            'StepOutputSelector': 1, 
            'StepNode': 1, 
            'str': 1
        }
        self.assertEqual(actual_counts, expected_counts)
    
    def test_trace_entity_counts5(self) -> None:
        step_id = 'stp5'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_counts = nextflow.trace.trace_entity_counts(src)
        print(actual_counts)
        expected_counts = {
            'StepTagInput': 1, 
            'Edge': 1, 
            'IndexOperator': 1, 
            'InputNodeSelector': 1, 
            'InputNode': 1, 
            'int': 1
        }
        self.assertEqual(actual_counts, expected_counts)

    def test_trace_source_datatype1(self) -> None:
        step_id = 'stp1'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = nextflow.trace.trace_source_datatype(src)
        expected_dtype = File
        self.assertIsInstance(actual_dtype, expected_dtype)
    
    def test_trace_source_datatype2(self) -> None:
        step_id = 'stp2'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = nextflow.trace.trace_source_datatype(src)
        expected_dtype = File
        self.assertIsInstance(actual_dtype, expected_dtype)
    
    def test_trace_source_datatype3(self) -> None:
        step_id = 'stp3'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = nextflow.trace.trace_source_datatype(src)
        expected_dtype = String
        self.assertIsInstance(actual_dtype, expected_dtype)
    
    def test_trace_source_datatype4(self) -> None:
        step_id = 'stp4'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = nextflow.trace.trace_source_datatype(src)
        expected_dtype = String
        self.assertIsInstance(actual_dtype, expected_dtype)
    
    def test_trace_source_datatype5(self) -> None:
        step_id = 'stp5'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = nextflow.trace.trace_source_datatype(src)
        expected_dtype = String
        self.assertIsInstance(actual_dtype, expected_dtype)

    def test_trace_source_scatter6(self) -> None:
        step_id = 'stp6'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_scatter = nextflow.trace.trace_source_scatter(src)
        expected_scatter = False
        self.assertEqual(actual_scatter, expected_scatter)

    def test_trace_source_scatter7(self) -> None:
        step_id = 'stp7'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_scatter = nextflow.trace.trace_source_scatter(src)
        expected_scatter = True
        self.assertEqual(actual_scatter, expected_scatter)


class TestPlumbingModule(unittest.TestCase):
    """tests the public functions in nfgen.plumbing."""

    def setUp(self) -> None:
        reset_globals()

    def test_get_array_depth(self):
        self.assertEqual(nextflow.plumbing.get_array_depth(String()), 0)
        self.assertEqual(nextflow.plumbing.get_array_depth(BamBai()), 0)
        self.assertEqual(nextflow.plumbing.get_array_depth(FastqGzPair()), 0)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(String())), 1)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(BamBai())), 1)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(FastqGzPair())), 1)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(Array(String()))), 2)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(Array(BamBai()))), 2)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(Array(FastqGzPair()))), 2)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(Array(Array(String())))), 3)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(Array(Array(BamBai())))), 3)
        self.assertEqual(nextflow.plumbing.get_array_depth(Array(Array(Array(FastqGzPair())))), 3)
    
    def test_is_array_depth_mismatch(self):
        # no mismatch
        self.assertFalse(nextflow.plumbing.is_array_depth_mismatch(String(), String(), False))
        self.assertFalse(nextflow.plumbing.is_array_depth_mismatch(Array(String()), Array(String()), False))
        self.assertFalse(nextflow.plumbing.is_array_depth_mismatch(Array(Array(String())), Array(Array(String())), False))
        # mismatch
        self.assertTrue(nextflow.plumbing.is_array_depth_mismatch(String(), Array(String()), False))
        self.assertTrue(nextflow.plumbing.is_array_depth_mismatch(Array(String()), String(), False))
        self.assertTrue(nextflow.plumbing.is_array_depth_mismatch(Array(Array(String())), Array(String()), False))

    # identifying common types ---
    def test_get_common_type_1(self):
        # non-union types, has intersection
        srctype = Bam()
        desttype = Bam()
        common_type = nextflow.plumbing.get_common_type(srctype, desttype)
        self.assertEqual(type(common_type), Bam)
        
    def test_get_common_type_2(self):
        # non-union types, no intersection
        srctype = Fasta()
        desttype = Bam()
        common_type = nextflow.plumbing.get_common_type(srctype, desttype)
        self.assertIsNone(common_type)
        
    def test_get_common_type_3(self):
        # single union type, has intersection
        srctype = UnionType(Bam, Sam, Cram)
        desttype = Bam()
        common_type = nextflow.plumbing.get_common_type(srctype, desttype)
        self.assertEqual(type(common_type), Bam)
        
    def test_get_common_type_4(self):
        # single union type, no intersection
        srctype = UnionType(Bam, Sam, Cram)
        desttype = Fasta()
        common_type = nextflow.plumbing.get_common_type(srctype, desttype)
        self.assertIsNone(common_type)
        
    def test_get_common_type_5(self):
        # both union types, has intersection
        srctype = UnionType(Bam, Sam, Cram)
        desttype = UnionType(Sam, Cram)
        common_type = nextflow.plumbing.get_common_type(srctype, desttype)
        self.assertEqual(type(common_type), Sam)
        
    def test_get_common_type_6(self):
        # both union types, no intersection
        srctype = UnionType(Bam, Sam, Cram)
        desttype = UnionType(Fasta, Fastq)
        common_type = nextflow.plumbing.get_common_type(srctype, desttype)
        self.assertIsNone(common_type)
    

    # array datatype mismatches ---
    # caused by datatypes involving arrays and/or scatter relationships
    def test_is_datatype_mismatch_arrays(self):
        """tests nfgen.plumibng.is_datatype_mismatch(srctype, desttype, srcscatter, destscatter)"""
        # secondary array (always considered mismatch)
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(Array(BamBai()), Array(BamBai()), False))

        # array depth - single single 
        self.assertFalse(nextflow.plumbing.is_datatype_mismatch(Bam(), Bam(), False)) # no mismatch
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(Bam(), Array(Bam()), False)) # mismatch
        
        # array depth - secondary secondary
        self.assertFalse(nextflow.plumbing.is_datatype_mismatch(BamBai(), BamBai(), False)) # no mismatch
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(Array(BamBai()), BamBai(), False)) # mismatch
        
        # array depth - secondary single 
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(BamBai(), Bam(), False)) # mismatch
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(FastaDict(), Bam(), False)) # mismatch
        
        # array depth - file pairs
        self.assertFalse(nextflow.plumbing.is_datatype_mismatch(FastqGzPair(), FastqGzPair(), False)) # no mismatch
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(FastqGzPair(), Array(FastqGzPair()), False)) # mismatch
        
    # def test_is_datatype_mismatch_scatter(self):
    #     # scatter - single single 
    #     self.assertFalse(nfgen.plumbing.is_datatype_mismatch(Array(Bam()), Bam(), True)) # no mismatch
    #     self.assertTrue(nfgen.plumbing.is_datatype_mismatch(Bam(), Bam(), True)) # mismatch
        
    #     # scatter - secondary secondary 
    #     self.assertFalse(nfgen.plumbing.is_datatype_mismatch(Array(BamBai()), BamBai(), True)) # no mismatch
    #     self.assertTrue(nfgen.plumbing.is_datatype_mismatch(BamBai(), BamBai(), True)) # mismatch

    #     # scatter - secondary single 
    #     self.assertTrue(nfgen.plumbing.is_datatype_mismatch(Array(BamBai()), Bam(), True)) # mismatch
    #     self.assertTrue(nfgen.plumbing.is_datatype_mismatch(FastaDict(), Array(Bam()), True)) # mismatch

    #     # scatter - file pairs
    #     self.assertFalse(nfgen.plumbing.is_datatype_mismatch(Array(FastqGzPair()), FastqGzPair(), True)) # no mismatch
    #     self.assertTrue(nfgen.plumbing.is_datatype_mismatch(FastqGzPair(), FastqGzPair(), True)) # mismatch
        
    def test_is_datatype_mismatch_basetype(self):
        # single single
        self.assertFalse(nextflow.plumbing.is_datatype_mismatch(Bam(), Bam(), False)) # no mismatch
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(Bam(), Fasta(), False)) # mismatch
        
        # secondary secondary
        self.assertFalse(nextflow.plumbing.is_datatype_mismatch(FastaDict(), FastaDict(), False)) # no mismatch
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(FastaDict(), BamBai(), False)) # mismatch
        
        # secondary single
        self.assertFalse(nextflow.plumbing.is_datatype_mismatch(Bam(), Bam(), False)) # no mismatch
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(Bam(), BamBai(), False)) # mismatch
        
        # file pairs
        self.assertFalse(nextflow.plumbing.is_datatype_mismatch(FastqGzPair(), FastqGzPair(), False)) # no mismatch
        self.assertTrue(nextflow.plumbing.is_datatype_mismatch(FastqGzPair(), BamBai(), False)) # mismatch

    def test_generate_datatype_mismatch_plumbing(self):
        # plumbing not required ---
        srctype, desttype, destscatter = Bam(), Bam(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '')
        
        srctype, desttype, destscatter = FastaWithIndexes(), FastaWithIndexes(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '')
        
        srctype, desttype, destscatter = FastqGzPair(), FastqGzPair(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '')

        srctype, desttype, destscatter = Array(Bam()), Array(Bam()), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '')
        
        # plumbing required ---
        
        # BASETYPES
        # secondary, secondary
        srctype, desttype, destscatter = FastaWithIndexes(), FastaDict(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.map{ tuple -> [tuple[0], tuple[4]] }')
        
        # secondary, single
        srctype, desttype, destscatter = FastaWithIndexes(), Fasta(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.map{ tuple -> tuple[0] }')

        # ARRAYS      
        # single, single  
        srctype, desttype, destscatter = Array(Bam()), Bam(), True
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten()')
        
        srctype, desttype, destscatter = Array(Fasta()), Fasta(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().first()')
        
        srctype, desttype, destscatter = Fasta(), Array(Fasta()), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.toList()')

        # any, secondary (always ends with .flatten().toList())
        srctype, desttype, destscatter = FastaWithIndexes(), Array(FastaWithIndexes()), True
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().toList()')
        
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Array(FastaWithIndexes()), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().toList()')
        
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Array(FastaWithIndexes()), True
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().toList()')

        # secondary, secondary 
        srctype, desttype, destscatter = Array(FastaWithIndexes()), FastaWithIndexes(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().collate( 8 ).first()')

        srctype, desttype, destscatter = FastaWithIndexes(), Array(FastaDict()), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.map{ tuple -> [tuple[0], tuple[4]] }.flatten().toList()')
        
        # secondary, single 
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Fasta(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.map{ tuple -> tuple[0] }.flatten().first()')

        # file pairs
        srctype, desttype, destscatter = Array(FastqGzPair()), FastqGzPair(), True
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().collate( 2 )')
        
        srctype, desttype, destscatter = Array(FastqGzPair()), FastqGzPair(), False
        plumbing = nextflow.plumbing.gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().collate( 2 ).first()')

        
        


class TestPlumbingTypeMismatch(unittest.TestCase):
    """
    This test group checks we can handle array / single datatype mismatches
    in a workflow. This occurs in wgsgermline. 
    """
    def setUp(self) -> None:
        self.wf = PlumbingTypeMismatchTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    def test_secondary_single_mismatch(self):
        step_id = 'bambai_to_bam'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_bam_bai.map{ tuple -> tuple[0] }'
        ])
        self.assertEqual(actual, expected)
    
    def test_secondary_single_array_mismatch(self):
        step_id = 'bambai_to_bam_array'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_bam_bai.map{ tuple -> tuple[0] }.toList()'
        ])
        self.assertEqual(actual, expected)
    
    def test_secondary_secondary_mismatch(self):
        step_id = 'fastawithindexes_to_fastadict'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_fasta_with_indexes.map{ tuple -> [tuple[0], tuple[4]] }'
        ])
        self.assertEqual(actual, expected)
    
    def test_secondary_secondary_array_mismatch(self):
        step_id = 'fastawithindexes_to_fastadict_array'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_fasta_with_indexes.map{ tuple -> [tuple[0], tuple[4]] }.flatten().toList()'
        ])
        self.assertEqual(actual, expected)

    def test_array_to_single(self):
        step_id = 'array_to_single'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_file_array.flatten().first()'
        ])
        self.assertEqual(actual, expected)
    
    def test_single_to_array(self):
        step_id = 'single_to_array'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_file.toList()'
        ])
        self.assertEqual(actual, expected)
    
    def test_secondary_array_to_secondary(self):
        step_id = 'secondary_array_to_secondary'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_bam_bai_array.flatten().collate( 2 ).first()'
        ])
        self.assertEqual(actual, expected)
    
    def test_to_secondary_array(self):
        step_id = 'secondary_array_to_secondary_array'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_bam_bai_array.flatten().toList()'
        ])
        self.assertEqual(actual, expected)



class TestPlumbingScatter(unittest.TestCase):
    """
    This test group checks whether data is being fed correctly
    to / between steps. 
    We need to ensure the following situations are handled:
        - secondary and secondary array types
        - scattering for each type
    """
    def setUp(self) -> None:
        self.wf = ComprehensiveScatterTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    def test_scatter_to_scatter(self):
        step_id = 'scatter_to_scatter'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        
        expected = set([
            'PRESTEP1.out.out'
        ])
        self.assertEqual(actual, expected)

    def test_scatter_to_array(self):
        step_id = 'scatter_to_array'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'PRESTEP1.out.out.toList()'
        ])
        self.assertEqual(actual, expected)

    def test_array_to_scatter1(self):
        step_id = 'prestep1'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_file_array.flatten()'
        ])
        self.assertEqual(actual, expected)
    
    def test_array_to_scatter2(self):
        step_id = 'array_to_scatter'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'PRESTEP2.out.out.flatten()'
        ])
        self.assertEqual(actual, expected)

    def test_scatter_secondary_to_scatter_secondary(self):
        step_id = 'scatter_secondary_to_scatter_secondary'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'PRESTEP3.out.out'
        ])
        self.assertEqual(actual, expected)

    def test_scatter_secondary_to_secondary_array(self):
        step_id = 'scatter_secondary_to_secondary_array'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'PRESTEP3.out.out.flatten().toList()'
        ])
        self.assertEqual(actual, expected)

    def test_secondary_array_to_scatter_secondary(self):
        step_id = 'secondary_array_to_scatter_secondary'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            'ch_in_bam_bai_array.flatten().collate( 2 )'
        ])
        self.assertEqual(actual, expected)

    # TODO use in future
    @unittest.skip('reimplement')
    def test_scatter_cross(self) -> None:
        wf = ScatterCrossTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))

        # cartesian cross channel manipulation in workflow
        operation = nextflow.channels.gen_scatter_cross_operation(step.sources, step.scatter)
        actual_op = operation.get_string()
        expected_op = """\
ch_in_file_array.flatten()
.combine(ch_in_str_array).flatten()
.multiMap { it ->
    in_file_array: it[0]
    in_str_array: it[1]
}
.set { ch_cartesian_cross }"""
        self.assertEqual(expected_op, actual_op)

        # step input values
        expected = [
            "ch_cartesian_cross.in_file_array",
            "ch_cartesian_cross.in_str_array",
        ]
        self.assertEqual(expected, actual)

    


class TestPlumbingBasic(unittest.TestCase):
    """
    This test group checks janis 'TInput' step inputs driven by workflow inputs
    or static values. 

    Ensures they are being handled correctly when parsed to nextflow.
    """

    def setUp(self) -> None:
        reset_globals()

    # workflow input step inputs
    def test_workflow_inputs(self):
        wf = StepInputsWFInputTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            "ch_in_file",
            "ch_in_file_opt",
        ])
        self.assertEqual(expected, actual)
        
    # static step inputs
    def test_static_inputs(self):
        wf = StepInputsTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        not_expected = {
            'pos_default',
            'pos_default',
            'pos_optional',
            'flag_true',
            'flag_false',
            'opt_basic',
            'opt_default',
            'opt_optional',
        }
        actual = nextflow.call.gen_task_call(step, scope, step_id)
        for tinput_name in not_expected:
            self.assertNotIn(tinput_name, actual)

    # connections
    def test_connections_files(self) -> None:
        wf = StepConnectionsTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        
        actual = split_task_call_to_lines(call)
        expected = ["STP1.out.out"]
        
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_filename_types(self) -> None:
        wf = FilenameTestWF1()
        do_preprocessing_workflow(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_file",
            "params.in_str",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
        
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_file",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_subworkflows1(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_file1",
            "ch_in_file_opt1",
            "params.in_str1",
            "params.in_str_opt1",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows2(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_file1",
            "ch_in_file_opt1",
            "'hello'",
            "'there'",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows3(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        
        step_id = 'stp3'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_file1",
            "'hello'",
            "null",
            "null",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows4(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        
        step_id = 'stp4'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_file1",
            "ch_in_file_opt1",
            "params.in_str1",
            "params.in_str_opt1",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows5(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        
        step_id = 'stp5'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_file1",
            "ch_in_file_opt1",
            "'hello'",
            "'there'",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows6(self) -> None:
        wf = DataSourceTestWF()
        do_preprocessing_workflow(wf)
        
        step_id = 'stp6'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_file1",
            "'hello'",
            "null",
            "null",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)




class TestPlumbingBasicArrays(unittest.TestCase):
    """
    This test group checks janis 'TInput' step inputs driven by workflow inputs
    or static values as above, but for Array situations. 

    Ensures they are being handled correctly when parsed to nextflow.
    """

    def setUp(self) -> None:
        reset_globals()

    def test_array_connections(self) -> None:
        wf = ArrayStepConnectionsTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            "STP1.out.out"
        ])
        self.assertEqual(expected, actual)
    
    def test_workflow_inputs_array(self) -> None:
        wf = ArrayStepInputsTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'stp1'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            "ch_in_file_array",
            "ch_in_file_array_opt",
        ])
        self.assertEqual(expected, actual)

    def test_static_step_inputs_array(self):
        wf = ArrayStepInputsTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'stp2'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        actual = nextflow.call.gen_task_call(step, scope, step_id)
        not_expected = {
            'pos_default',
            'pos_default',
            'pos_optional',
            'flag_true',
            'flag_false',
            'opt_basic',
            'opt_default',
            'opt_optional',
        }
        for tinput_name in not_expected:
            self.assertNotIn(tinput_name, actual)



# NOTE: previous secondaries / secondaries array format.
# class TestPlumbingSecondaries(unittest.TestCase):
#     # TODO test more features / edge cases?
#     """
#     This test group checks that janis File types with secondaries are being handled 
#     correctly to produce the desired nextflow workflow. 
    
#     SecondaryFiles are especially tricky when mapping to janis. 
#     Requires the remapping of tool inputs in some situations.
    
#     eg Array(BamBai) 
#         - Arrays where each item is a file with secondaries. 
#         - These gets remapped from a single tool input, to an array input
#           for each separate file in the datatype. 
#         - Works like transposing - instead of inp1=[[Bam, Bai], [Bam, Bai]],
#           results in inp1=[Bam, Bam], inp2=[Bai, Bai]

#     """

#     def setUp(self) -> None:
    
#     def test_secondaries_workflow_inputs(self) -> None:
#         wf = SecondariesTestWF()
#         refresh_workflow_inputs(wf)
#         # params created correctly?
#         actual_params = nfgen.params.serialize()
#         expected_params = {
#             'in_alignments_bam': 'null',
#             'in_alignments_bai': 'null',
#         }
#         self.assertDictContainsSubset(expected_params, actual_params)
        
#         # channels created correctly?
#         inp = wf.input_nodes['inAlignments']
#         ch_expected = 'ch_in_alignments = Channel.fromPath( [params.in_alignments_bam, params.in_alignments_bai] ).toList()'
#         ch_actual = nfgen.nextflow.channels.getall(inp.uuid)[0]
#         ch_actual = ch_actual.definition
#         self.assertEquals(ch_actual, ch_expected)
        
#         # step inputs created correctly?
#         step_id = 'stp1'
#         step = wf.step_nodes[step_id]
#         scope = nfgen.Scope()
#         scope.update(step)
#         actual_inputs = nfgen.call.get_args(step, scope)
#         expected_inputs = [
#             "ch_in_alignments"
#         ]
#         self.assertEquals(actual_inputs, expected_inputs)
    
#     def test_secondaries_array_workflow_inputs(self) -> None:
#         wf = SecondariesTestWF()
#         refresh_workflow_inputs(wf)
#         # params created correctly?
#         actual_params = nfgen.params.serialize()
#         expected_params = {
#             'in_alignments_arr_bams': '[]',
#             'in_alignments_arr_bais': '[]',
#         }
#         self.assertDictContainsSubset(expected_params, actual_params)
        
#         # channels created correctly?
#         inp = wf.input_nodes['inAlignmentsArr']
#         ch1_expected = 'ch_in_alignments_arr_bais = Channel.fromPath( params.in_alignments_arr_bais ).toList()'
#         ch2_expected = 'ch_in_alignments_arr_bams = Channel.fromPath( params.in_alignments_arr_bams ).toList()'
#         ch1_actual, ch2_actual = nfgen.nextflow.channels.getall(inp.uuid)
#         ch1_actual, ch2_actual = ch1_actual.definition, ch2_actual.definition
#         self.assertEquals(ch1_actual, ch1_expected)
#         self.assertEquals(ch2_actual, ch2_expected)
        
#         # step inputs created correctly?
#         step_id = 'stp4'
#         step = wf.step_nodes[step_id]
#         scope = nfgen.Scope()
#         scope.update(step)
#         actual_inputs = nfgen.call.get_args(step, scope)
#         expected_inputs = [
#             "ch_in_alignments_arr_bais",
#             "ch_in_alignments_arr_bams",
#         ]
#         self.assertEquals(actual_inputs, expected_inputs)
    
#     def test_secondaries_connections(self) -> None:
#         wf = SecondariesTestWF()
#         refresh_workflow_inputs(wf)
#         # step inputs created correctly?
#         step_id = 'stp2'
#         step = wf.step_nodes[step_id]
#         scope = nfgen.Scope()
#         scope.update(step)
#         actual_inputs = nfgen.call.get_args(step, scope)
#         expected_inputs = [
#             "STP1.out.out"
#         ]
#         self.assertEquals(actual_inputs, expected_inputs)




class TestPlumbingCombinations(unittest.TestCase):
    """
    Tests plumbing for most complex cases.
    Includes combinations of scatter, arrays, secondary files. 
    """
    def setUp(self) -> None:
        reset_globals()

    # scatter + secondaries
    @unittest.skip('not implemented')
    def test_scatter_secondaries_workflow_inputs(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_secondaries_dot(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_secondaries_cross(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_secondaries_connections(self) -> None:
        raise NotImplementedError

    # secondaries + arrays
    @unittest.skip('not implemented')
    def test_secondaries_workflow_inputs_array(self) -> None:
        # wf = ArraySecondariesTestWF()
        # actual_params = refresh_workflow_inputs(wf)
        # expected_params = {
        #     'in_alignments_bams': '[]',
        #     'in_alignments_bais': '[]',
        # }
        # self.assertEquals(actual_params, expected_params)
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_secondaries_connections_array(self) -> None:
        raise NotImplementedError
    
    # scatter + arrays
    @unittest.skip('not implemented')
    def test_scatter_basic_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_dot_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_cross_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_connection_array(self) -> None:
        raise NotImplementedError



class TestPlumbingEdgeCases(unittest.TestCase):
    """
    This test group checks we can handle array / single datatype mismatches
    in a workflow. This occurs in wgsgermline. 
    """
    def setUp(self) -> None:
        self.wf = PlumbingEdgeCaseTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    def test_pythontool_array_string_output(self) -> None:
        step_id = 'stp2'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            "STP1.out.out.filter{ it != '' }.map{ it -> it.split(', ') }.ifEmpty( null )"
        ])
        self.assertEqual(actual, expected)



class TestWorkflowOutputs(unittest.TestCase):
    """
    Tests workflow outputs being created correctly
    """
    def setUp(self) -> None:
        reset_globals()
        
    @unittest.skip('not implemented')
    def test_files(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_files_arrays(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_secondaries(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_secondaries_array(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_scatter_secondaries(self) -> None:
        raise NotImplementedError




class TestConfig(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()
    
    def test_outdir(self):
        wf = BasicIOTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        print(config)
        expected_values = {
            'outdir': './outputs',
        }
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_file(self):
        wf = BasicIOTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        print(config)
        expected_values = {
            'in_file': 'null',
            'in_str': 'null',
        }
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_nonfile(self):
        # string, int, bool
        wf = StepInputsTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        expected_values = {
            'in_file': 'null',
            'in_str': 'null',
            'in_int': 'null',
            'in_bool': 'null',
            # 'stp2_pos_default': '100',
            # 'stp2_pos_optional': "'static'",
            # 'stp2_flag_true': 'false',
            # 'stp2_flag_false': 'true',
            # 'stp2_opt_basic': "'static'",
            # 'stp2_opt_default': '100',
            # 'stp2_opt_optional': "''",
            # 'stp3_opt_basic': "'static'",
            # 'stp3_opt_default': '100',
            # 'stp3_opt_optional': "''",
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)

    def test_array(self):
        wf = ArrayIOExtrasTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        # I apologise for disgusting formatting below. 
        # This can probably be handled better.
        expected_values = {
"in_file_array": f"""[
    // list files here
]""",
"in_str_array": f"""[
    // list values here
]""",
"in_int_array": f"""[
    // list values here
]""",
# "in_float_array": f"""[
#     // list values here
# ]""",
# "stp4_inp": f"""[
#     'hello',
#     'there!'
# ]""",
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_file_pair(self):
        wf = FilePairsTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        print(config)
        expected_values = {
"in_reads": f"""[
    // read 1
    // read 2
]""",
        }
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_array_file_pair(self):
        wf = FilePairsTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        print(config)
        expected_values = {
"in_reads": f"""[
    // read 1
    // read 2
]""",
        }
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)

    def test_secondary(self):
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        expected_values = {
"in_alignments": f"""[
    // bam
    // bai
]""",
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)
    
    def test_array_secondary(self):
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        expected_values = {
"in_alignments_arr": f"""[
        [
            // bam
            // bai
        ],
    ]""",
        }
        print(config)
        for name, val in expected_values.items():
            self.assertIn(name, config)
            self.assertIn(val, config)
    
    def test_nonfile_array_workflow_inputs(self):
        # string, int, bool
        wf = ArrayStepInputsTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        print(config)
        expected_lines = {
            "outdir  = './outputs'",
            "in_file_array      = [",
            "    // list files here",
            "]",
            "in_file_array_opt  = [",
            "    // list files here",
            "]",
            "in_str_array       = []  // list values here",
            "in_int_array       = []  // list values here",
            # "stp2_pos_default   = [4, 5, 6]",
            # "stp2_pos_optional  = ['hi', 'there', 'friend']",
            # "stp2_opt_basic     = ['hi', 'there', 'friend']",
            # "stp2_opt_default   = [4, 5, 6]",
            # "stp2_opt_optional  = ['hi', 'there', 'friend']",
        }
        for ln in expected_lines:
            self.assertIn(ln, config)

    def test_workflow_config(self) -> None:
        wf = AssemblyTestWF()
        do_preprocessing_workflow(wf)
        # params = translator.build_inputs_dict(wf)
        config = translator.stringify_translated_inputs({})
        # basic structure
        self.assertIn('docker.enabled = true', config)
        self.assertIn('params {\n\n', config)
        self.assertIn('\n\n}', config)
        # expected params are present & correct
        expected_values = {
            'fastqc1_adapters': 'NO_FILE1',
            'fastqc1_contaminants': 'NO_FILE2',
            'fastqc1_limits': 'NO_FILE3',
            'fastqc2_adapters': 'NO_FILE4',
            'fastqc2_contaminants': 'NO_FILE5',
            'fastqc2_limits': 'NO_FILE6',
            'in_forward_reads': 'null',
            'in_long_reads': 'null',
            'in_reverse_reads': 'null',
            'test_input': 'null',
            # 'unicycler_kmers': "''",
            # 'unicycler_scores': "''",
            # 'unicycler_start_gene_cov': '95.0',
            # 'unicycler_start_gene_id': '90.0',
        }
        print(config)
        for name, val in expected_values.items():
            pattern = f'{name}.*?{val}'
            matches = re.findall(pattern, config)
            self.assertGreater(len(matches), 0)

    @unittest.skip('not implemented')
    def test_tool_config(self) -> None:
        raise NotImplementedError




class TestPlumbingExpressions(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_first_selector(self):
        wf = ConditionStepTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'print'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            "[params.mystring, GET_STRING.out.out].find{ it != null }"
        ])
        self.assertEqual(actual, expected)

    def test_with_expression(self):
        wf = StepInputExpressionTestWF()
        do_preprocessing_workflow(wf)
        step_id = 'print'
        step = wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([
            "params.mystring ? params.mystring : params.mystring_backup"
        ])
        self.assertEqual(actual, expected)



class TestUnwrapProcess(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()
        wf = UnwrapTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(step.tool, step.sources)
        self.prescript, self.script = nextflow.generate.process.gen_nf_process_script(
            scope=scope,
            tool=step.tool,
            variable_manager=variable_manager,
            sources=step.sources,
            stdout_filename='out'
        )

    # PROCESS RELATED
    def test_filename_generated(self) -> None:
        self.assertIn("--filenameGen generated.gz", self.script)
    
    def test_filename_reference(self) -> None:
        self.assertIn("--filenameRef ${in_file.simpleName}.fastq.gz", self.script)

    def test_input_selector_process_input(self) -> None:
        self.assertIn("--inputSelectorProcess ${in_file}", self.script)

    def test_input_selector_param_input(self) -> None:
        self.assertIn("--inputSelectorParam ${params.in_str}", self.script)
       
    def test_list(self) -> None:
        self.assertIn("--list [1, 2, 3, 4, 5]", self.script) # this is correct
       
    def test_two_value_operator(self) -> None:
        self.assertIn("--TwoValueOperator ${in_file + \".gz\"}", self.script)
    
    def test_first_operator(self) -> None:
        self.assertIn("--FirstOperator ${[params.in_str, []].find{ it != null }}", self.script)
    
    def test_index_operator(self) -> None:
        self.assertIn("--IndexOperator ${in_file_arr[0]}", self.script)
    
    def test_index_operator_secondaries(self) -> None:
        self.assertIn("--IndexOperatorSecondariesBam ${bam}", self.script)
        self.assertIn("--IndexOperatorSecondariesBai ${bai}", self.script)
    
    def test_index_operator_secondaries_array(self) -> None:
        print(self.script)
        self.assertIn("--IndexOperatorArraySecondariesBams ${indexed_bam_flat[0]}", self.script)
        self.assertIn("--IndexOperatorArraySecondariesBais ${indexed_bam_flat[1]}", self.script)


class TestUnwrapWorkflow(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_input_node_channel(self) -> None:
        # file input (channel)
        wf = StepConnectionsTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        node = wf.input_nodes['inFile']
        actual = nextflow.unwrap_expression(node, scope=scope, context='workflow')
        expected = 'ch_in_file'
        self.assertEqual(actual, expected)

    def test_input_node_param(self) -> None:
        # nonfile input (param)
        wf = StepConnectionsTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        node = wf.input_nodes['inStr']
        actual = nextflow.unwrap_expression(node, scope=scope, context='workflow')
        expected = 'params.in_str'
        self.assertEqual(actual, expected)

    def test_step_connection(self) -> None:
        wf = StepConnectionsTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        step_id = "stp2"
        inp_id = 'inp'
        sources = wf.step_nodes[step_id].sources
        src = sources[inp_id]
        actual = nextflow.unwrap_expression(src, scope=scope, context='workflow')
        expected = 'STP1.out.out'
        self.assertEqual(actual, expected)
    
    def test_alias_selector_wf(self) -> None:
        wf = AliasSelectorTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        step_id = "stp2"
        inp_id = 'inp'
        sources = wf.step_nodes[step_id].sources
        src = sources[inp_id]
        actual = nextflow.unwrap_expression(src, scope=scope, context='workflow')
        expected = 'STP1.out.out'
        self.assertEqual(actual, expected)
    
    def test_first_operator_wf(self) -> None:
        wf = ConditionStepTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        step_id = "print"
        inp_id = 'inp'
        sources = wf.step_nodes[step_id].sources
        src = sources[inp_id]
        actual = nextflow.unwrap_expression(src, scope=scope, context='workflow')
        expected = '[params.mystring, GET_STRING.out.out].find{ it != null }'
        self.assertEqual(actual, expected)
    
    def test_index_operator_wf(self) -> None:
        wf = IndexOperatorTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        step_id = "stp1"
        inp_id = 'inp'
        sources = wf.step_nodes[step_id].sources
        src = sources[inp_id]
        actual = nextflow.unwrap_expression(src, scope=scope, context='workflow')
        expected = 'ch_in_file_arr[0]'
        self.assertEqual(actual, expected)


class TestUnwrapStringFormatter(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()
    
    def test_string_formatter(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = BasicTestTool()
        do_preprocessing_tool(tool)
        scope = nextflow.Scope()
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(tool)
        sf = StringFormatter("no format")
        actual = nextflow.unwrap_expression(
            val=sf,
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual('"no format"', actual)

    def test_string_formatter_string(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = BasicTestTool()
        do_preprocessing_tool(tool)
        scope = nextflow.Scope()
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(tool)
        sf = StringFormatter("there's a {str_arg} arg", str_arg="string")
        actual = nextflow.unwrap_expression(
            val=sf,
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual('"there\'s a string arg"', actual)
    
    def test_string_formatter_inputselector_process_input(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = BasicTestTool()
        do_preprocessing_tool(tool)
        scope = nextflow.Scope()
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(tool)
        sf = StringFormatter("an input {arg}", arg=InputSelector("testtool"))
        actual = nextflow.unwrap_expression(
            val=sf,
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual('"an input ${testtool}"', actual)
    
    def test_string_formatter_inputselector_param_input(self):
        wf = StringFormatterTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(step.tool, step.sources)
        sf = StringFormatter("an input {arg}", arg=InputSelector("compressionLevel"))
        actual = nextflow.unwrap_expression(
            val=sf, 
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=step.tool,
            in_shell_script=True
        )
        expected = '"an input ${5}"'
        self.assertEqual(actual, expected)

    def test_string_formatter_two_param(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = InputQualityTestTool()
        do_preprocessing_tool(tool)
        scope = nextflow.Scope()
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(tool)
        sf = StringFormatter(
            "{username}:{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        actual = nextflow.unwrap_expression(
            val=sf,
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual('"${user}:${static}"', actual)

    def test_escaped_characters(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = InputQualityTestTool()
        do_preprocessing_tool(tool)
        scope = nextflow.Scope()
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(tool)
        sf = StringFormatter(
            "{username}\\t{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        actual_scripting = nextflow.unwrap_expression(
            val=sf,
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=False
        )
        actual_shell = nextflow.unwrap_expression(
            val=sf,
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual('user\\tstatic', actual_scripting)
        self.assertEqual('"${user}\\\\t${static}"', actual_shell)

    def test_expression_arg(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = BasicTestTool()
        do_preprocessing_tool(tool)
        scope = nextflow.Scope()
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(tool)
        sf = StringFormatter(
            "{name}:{items}",
            name=InputSelector("testtool"),
            items=JoinOperator(InputSelector("arrayInp"), separator=";"),
        )
        res = nextflow.unwrap_expression(
            val=sf,
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual('"${testtool}:${array_inp.join(\";\")}"', res)
    
    def test_string_formatter_advanced(self) -> None:
        wf = StringFormatterTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        variable_manager = VariableManager(scope)
        variable_manager.update_for_tool(step.tool, step.sources)
        arg = step.tool.arguments()[0]
        actual = nextflow.unwrap_expression(
            val=arg.value, 
            scope=scope,
            context='process_script',
            variable_manager=variable_manager,
            tool=step.tool,
            in_shell_script=True
        )
        expected = '"-Xmx${8 * 3 / 4}G ${5 ? "-Dsamjdk.compress_level=" + 5 : ""} ${["EXTRA_ARG1", []].find{ it != null }.join(" ")}"'
        self.assertEqual(actual, expected)



class TestOrdering(unittest.TestCase):

    def setUp(self) -> None:
        self.wf = OrderingTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    def test_process_call(self) -> None:
        # from workflow inputs
        step_id = 'stp1'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_fastq",
            "ch_in_fastq_array",
            "ch_in_file",
        ]
        self.assertEqual(expected, actual)
        
        # from process outputs
        step_id = 'stp3'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "STP1.out.outFastq",
            "STP1.out.outFastqArray",
            "STP1.out.outFile",
            "STP1.out.outInt",
            "STP1.out.outIntArray",
            "STP1.out.outStr",
        ]
        self.assertEqual(expected, actual)

        # from subworkflow outputs
        step_id = 'stp5'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "STP2.out.outFastq",
            "STP2.out.outFastqArray",
            "STP2.out.outFile",
            "STP2.out.outInt",
            "STP2.out.outIntArray",
            "STP2.out.outStr",
        ]
        self.assertEqual(expected, actual)

    def test_process_inputs(self) -> None:
        # from workflow inputs
        step = self.wf.step_nodes["stp1"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        expected_inputs = [
            "path in_fastq, stageAs: 'in_fastq.fastq'",
            "path in_fastq_array",
            "path in_file, stageAs: 'in_file'"
        ]

        actual_inputs = [inp.get_string() for inp in process.inputs]
        self.assertEqual(actual_inputs, expected_inputs)

        # from process outputs
        step = self.wf.step_nodes["stp3"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        expected_inputs = [
            "path in_fastq, stageAs: 'in_fastq.fastq'",
            "path in_fastq_array",
            "path in_file, stageAs: 'in_file'",
            "val in_int",
            "val in_int_array",
            "val in_str"
        ]
        actual_inputs = [inp.get_string() for inp in process.inputs]
        self.assertEqual(actual_inputs, expected_inputs)

        # from subworkflow outputs
        step = self.wf.step_nodes["stp5"]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool, step.sources)
        expected_inputs = [
            "path in_fastq, stageAs: 'in_fastq.fastq'",
            "path in_fastq_array",
            "path in_file, stageAs: 'in_file'",
            "val in_int",
            "val in_int_array",
            "val in_str"
        ]
        actual_inputs = [inp.get_string() for inp in process.inputs]
        self.assertEqual(actual_inputs, expected_inputs)

    def test_subworkflow_call(self) -> None:
        # from workflow inputs
        step_id = 'stp2'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "ch_in_fastq",
            "ch_in_fastq_array",
            "ch_in_file",
        ]
        self.assertEqual(expected, actual)
        
        # from process outputs
        step_id = 'stp4'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "STP1.out.outFastq",
            "STP1.out.outFastqArray",
            "STP1.out.outFile",
            "STP1.out.outInt",
            "STP1.out.outIntArray",
            "STP1.out.outStr",
        ]
        self.assertEqual(expected, actual)

        # from subworkflow outputs
        step_id = 'stp6'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = split_task_call_to_lines(call)
        expected = [
            "STP2.out.outFastq",
            "STP2.out.outFastqArray",
            "STP2.out.outFile",
            "STP2.out.outInt",
            "STP2.out.outIntArray",
            "STP2.out.outStr",
        ]
        self.assertEqual(expected, actual)
    

    def test_subworkflow_inputs(self) -> None:
        # from workflow inputs
        step_id = 'stp2'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        scope.update(step)
        
        nf_workflow = nextflow.generate.main.generate_workflow(
            name=step_id, 
            scope=scope,
            wf=step.tool,
            item_register=translator.item_register
        )
        expected_inputs = [
            "ch_in_fastq",
            "ch_in_fastq_array",
            "ch_in_file",
        ]
        actual_inputs = [ch.get_string() for ch in nf_workflow.take]
        self.assertEqual(actual_inputs, expected_inputs)

        # from process outputs
        step_id = 'stp4'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        scope.update(step)
        nf_workflow = nextflow.generate.main.generate_workflow(
            name=step_id, 
            scope=scope,
            wf=step.tool,
            item_register=translator.item_register
        )
        expected_inputs = [
            "ch_in_fastq",
            "ch_in_fastq_array",
            "ch_in_file",
            "ch_in_int",
            "ch_in_int_array",
            "ch_in_str"
        ]
        actual_inputs = [ch.get_string() for ch in nf_workflow.take]
        self.assertEqual(actual_inputs, expected_inputs)

        # from subworkflow outputs
        step_id = 'stp6'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        scope.update(step)
        nf_workflow = nextflow.generate.main.generate_workflow(
            name=step_id, 
            scope=scope,
            wf=step.tool,
            item_register=translator.item_register
        )
        expected_inputs = [
            "ch_in_fastq",
            "ch_in_fastq_array",
            "ch_in_file",
            "ch_in_int",
            "ch_in_int_array",
            "ch_in_str"
        ]
        actual_inputs = [ch.get_string() for ch in nf_workflow.take]
        self.assertEqual(actual_inputs, expected_inputs)

    @unittest.skip('not implemented')
    def test_workflow_imports(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_channel_definitions(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_config(self) -> None:
        raise NotImplementedError
    



class TestSubWorkflows(unittest.TestCase):
    # sometimes logic here a bit more complex or weird 
    # due to subworkflows. apologies. 
    
    def setUp(self) -> None:
        self.wf = SubworkflowTestWF()
        do_preprocessing_workflow(self.wf)
        reset_globals()

    @unittest.skip('not implemented')
    def test_param_system(self) -> None:
        # currently params system doesnt reach to subworkflows. 
        # will implement in future (time permitting)
        raise NotImplementedError
       
    def test_files_created(self) -> None:
        do_preprocessing_workflow(self.wf)
        _, substr_dict = translator.translate_workflow_internal(self.wf)
        expected_filepaths = set([
            'modules/file_tool',
            'modules/string_tool',
            'modules/int_tool',
            'modules/string_opt_tool',
            'subworkflows/oranges_subworkflow',
            'subworkflows/apples_subworkflow',
        ])
        actual_filepaths = set(substr_dict.keys())
        self.assertEqual(actual_filepaths, expected_filepaths)

    def test_structure(self) -> None:
        # take main emit
        mainstr, substr_dict = translator.translate_workflow_internal(self.wf)
        self.assertNotIn('take:', mainstr)
        self.assertNotIn('main:', mainstr)
        self.assertNotIn('emit:', mainstr)
        subwfstr = substr_dict['subworkflows/oranges_subworkflow']
        self.assertIn('take:', subwfstr)
        self.assertIn('main:', subwfstr)
        self.assertIn('emit:', subwfstr)

    def test_call(self) -> None:
        # translate workflow, building all nf items and files
        translator.translate_workflow_internal(self.wf)

        # call args are correct & in order
        step_id = 'apples_subworkflow'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()

        call = nextflow.call.gen_task_call(step, scope, step_id)
        actual = set(split_task_call_to_lines(call))
        expected = set([])
        self.assertSetEqual(actual, expected)

        # # subworkflow inputs are correct & in order
        # nf_workflow = nextflow.parsing.workflow.gen_workflow(
        #     name=step_id, 
        #     scope=scope,
        #     wf=step.tool,
        #     item_register=translator.item_register
        # )
        # # check the arg order matches the subworkflow input channel order
        # actual_channel_order = [t.get_string() for t in nf_workflow.take]
        # expected_channel_order = ['ch_in_int', 'ch_in_str', 'ch_in_str_opt']
        # self.assertEquals(actual_channel_order, expected_channel_order)

    def test_imports(self) -> None:
        # translate workflow, building all nf items and files
        wf = SubworkflowTestWF()
        translator.translate_workflow_internal(wf)

        # focusing in on specific subworkflow
        step = wf.step_nodes['apples_subworkflow']
        scope = nextflow.Scope()
        scope.update(step)
        subwf_file = translator.file_register.get(scope)
        actual_imports = subwf_file.items[0].imports
        expected_imports = [
            "include { STRING_TOOL } from '../modules/string_tool'",
            "include { STRING_OPT_TOOL } from '../modules/string_opt_tool'",
            "include { ORANGES_SUBWORKFLOW } from './oranges_subworkflow'"
        ]
        self.assertEquals(len(actual_imports), 3)
        self.assertEquals(actual_imports, expected_imports)
    
    @unittest.skip('not implemented')
    def test_nested(self) -> None:
        raise NotImplementedError
    

class TestNaming(unittest.TestCase):
    
    def setUp(self) -> None:
        self.wf = NamingTestWF()
        reset_globals()
        do_preprocessing_workflow(self.wf)

    def test_workflow(self) -> None:
        name = nextflow.naming.constructs.gen_varname_workflow(self.wf.id())
        self.assertEqual(name, 'NAMING_TEST_WF')
    
    def test_process(self) -> None:
        step = self.wf.step_nodes['stp1']
        name = nextflow.naming.constructs.gen_varname_process(step.id())
        self.assertEqual(name, 'STP1')

    def test_channels(self) -> None:
        scope = nextflow.Scope()
        channels = nextflow.channels.getall(scope)
        channel_names = {ch.name for ch in channels}
        expected_names = {
            'ch_process_input',
            'ch_process_input_array',
            'ch_secondary',
            'ch_panel_of_normals',
            'ch_reference',
            'ch_tumor_sample',
            'ch_germline_resource',
            'ch_normal_sample',
        }
        self.assertEqual(channel_names, expected_names)
    
    def test_params(self) -> None:
        params = nextflow.params.getall()
        param_names = {p.name for p in params}
        expected_names = {
            'outdir',
            'process_input',
            'process_input_array',
            'param_input',
            'param_input_array',
            'secondary',
            'normal_sample',
            'panel_of_normals',
            'reference',
            'germline_resource',
            'tumor_sample',
        }
        self.assertEqual(param_names, expected_names)
    
    def test_process_inputs(self) -> None:
        step = self.wf.step_nodes['stp1']
        scope = nextflow.Scope()
        scope.update(step)
        process_inputs = nextflow.generate.process.gen_nf_process_inputs(scope, step.tool)
        actual_input_names = [x.name for x in process_inputs]
        expected_input_names = [
            'secondary',
            'indexed_bam_flat',
            'process_input',
            'process_input_array',
        ]
        self.assertEqual(actual_input_names, expected_input_names)
        secondary_input = [x for x in process_inputs if x.name == 'secondary'][0]
        self.assertEqual(secondary_input.subnames, ['bam', 'bai'])

    def test_process_inputs_hard(self) -> None:
        step = self.wf.step_nodes['stp3']
        scope = nextflow.Scope()
        scope.update(step)
        process_inputs = nextflow.generate.process.gen_nf_process_inputs(scope, step.tool)
        
        # checking the names of process inputs are correct
        actual_input_names = set([x.name for x in process_inputs])
        expected_input_names = set([
            'normal_samples_flat',
            'tumor_samples_flat',
            'reference',
            'germlineResource',
            'panelOfNormals',
        ])
        
        # checking the subnames of process inputs are correct (for secondary files)
        self.assertEqual(actual_input_names, expected_input_names)

        # reference (FastaWithIndexes)
        process_input = [x for x in process_inputs if x.name == 'reference'][0]
        expected_subnames = ['fasta', 'amb', 'ann', 'bwt', 'dict', 'fai', 'pac', 'sa']
        self.assertEqual(process_input.subnames, expected_subnames)
        
        # panelOfNormals (VcfTabix) #1
        process_input = [x for x in process_inputs if x.name == 'panelOfNormals'][0]
        expected_subnames = ['panel_of_normals_vcf_gz', 'panel_of_normals_tbi']
        self.assertEqual(process_input.subnames, expected_subnames)
        
        # germlineResource (VcfTabix) #2
        process_input = [x for x in process_inputs if x.name == 'germlineResource'][0]
        expected_subnames = ['germline_resource_vcf_gz', 'germline_resource_tbi']
        self.assertEqual(process_input.subnames, expected_subnames)
        
        # normal_samples (Array(BamBai())) #1
        process_input = [x for x in process_inputs if x.name == 'normal_samples_flat'][0]
        expected_name = 'normal_samples_flat'
        self.assertEqual(process_input.name, expected_name)
        
        # tumor_samples (Array(BamBai())) #2
        process_input = [x for x in process_inputs if x.name == 'tumor_samples_flat'][0]
        expected_name = 'tumor_samples_flat'
        self.assertEqual(process_input.name, expected_name)

    def test_process_outputs(self) -> None:
        # TODO does not include secondaries array outputs 
        step_id = 'stp1'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_cmdtool(tool=step.tool, sources=step.sources, scope=scope)
        
        # process output names correct?
        actual_output_names = [x.name for x in process.outputs]
        expected_output_names = [
            'outProcessInput',
            'outParamInput',
            'outSecondary',
            'outProcessInputArray',
            'outParamInputArray',
        ]
        self.assertEqual(actual_output_names, expected_output_names)
        secondary_output = [x for x in process.outputs if x.name == 'outSecondary'][0]
        self.assertTrue(secondary_output)
        self.assertEqual(secondary_output.expressions, ['"*.bam"', '"*.bam.bai"'])

    def test_process_references(self) -> None:
        # singles & arrays
        step_id = 'stp1'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        scope.update(step)
        process = nextflow.generate.process.gen_process_from_cmdtool(tool=step.tool, sources=step.sources, scope=scope)
        print(process.get_string())
        
        # pre-script
        assert(process.pre_script)
        actual_prescript = process.pre_script.split('\n')
        expected_prescript = [
            'def process_input_array_joined = process_input_array.join(\' \')',
            'def param_input_array_joined = params.param_input_array.join(\' \')',
            "def secondary_array = get_primary_files(indexed_bam_flat)",
            'def secondary_array_joined = secondary_array.join(\' \')',
        ]
        self.assertEqual(actual_prescript, expected_prescript)

        # script body
        actual_script = process.script.split('\n')
        expected_script = [
            'echo \\',
            '--processInput ${process_input} \\',
            '--paramInput ${params.param_input} \\',
            '--secondary ${bam} \\',
            '--processInputArray ${process_input_array_joined} \\',
            '--paramInputArray ${param_input_array_joined} \\',
            '--secondaryArray ${secondary_array_joined} \\',
        ]
        self.assertEqual(actual_script, expected_script)

        # outputs
        actual_outputs = set([x.get_string() for x in process.outputs])
        expected_outputs = set([
            'path "${process_input + ".fastq"}", emit: outProcessInput',
            'path "param_input.txt", emit: outParamInput',
            'tuple path("*.bam"), path("*.bam.bai"), emit: outSecondary',
            'path "process_input_arr*", emit: outProcessInputArray',
            'path "param_input_arr*", emit: outParamInputArray'
        ])
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_connections(self) -> None:
        # ensure name references are correct between steps
        step_id = 'stp2'
        step = self.wf.step_nodes[step_id]
        scope = nextflow.Scope()
        task_call = nextflow.call.gen_task_call(step, scope, step_id)
        actual_inputs = set(split_task_call_to_lines(task_call))
        expected_inputs = set([
            'STP1.out.outProcessInput',
            'STP1.out.outProcessInputArray',
            'STP1.out.outSecondary',
            'ch_normal_sample.flatten().toList()',
            'STP1.out.outParamInput',
            'STP1.out.outParamInputArray',
        ])
        self.assertEqual(actual_inputs, expected_inputs)
