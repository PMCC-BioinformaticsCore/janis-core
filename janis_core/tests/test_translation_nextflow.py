import unittest

from typing import Any

from janis_core.tests.testtools import (
    InputQualityTestTool,
    BasicTestTool,
    FastqcTestTool,
    BwaMemTestTool,
    GridssTestTool,
    FileOutputPythonTestTool,
    MultiTypesInputPythonTool,
    SecondaryInputPythonTestTool,
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
    AllInputTypesTestWF,
    FilesDirectoriesToCreateTestWF
) 

from janis_core import (
    Workflow,
    CommandTool,
    PythonTool,
    InputSelector,
    StringFormatter,
    JoinOperator,
    TInput,
    Tool
)

from janis_core.translations import translate
from janis_core.translations import NextflowTranslator as translator
from janis_core.translations import nextflow

from janis_core.translations.nextflow.variables import VariableManager
from janis_core.translations.nextflow.variables import VariableType
from janis_core.translations.nextflow.variables import init_variable_manager_for_task

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
from janis_core import translation_utils as utils
from janis_core import settings

from janis_core.translations.nextflow.generate.workflow.common import get_common_type
from janis_core.translations.nextflow.generate.workflow.datatype_mismatch import (
    get_array_depth,
    is_datatype_mismatch,
    is_array_depth_mismatch,
    gen_datatype_mismatch_plumbing,
)


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
    nextflow.params.clear()

def do_preprocessing_workflow(wf: Workflow) -> None:
    reset_globals()
    nextflow.preprocessing.populate_task_inputs_workflowmode(wf, wf)

def do_preprocessing_tool(tool: CommandTool | PythonTool) -> None:
    reset_globals()
    nextflow.preprocessing.populate_task_inputs_toolmode(tool)

def split_to_lines(call: str) -> list[str]:
    lines = call.split('\n')                        # split
    lines = [ln.strip(' ') for ln in lines]         # strip indents
    lines = [ln for ln in lines if not ln == '']    # remove empty lines
    return lines

def simplify_call(textlines: list[str]) -> list[str]:
    lines = textlines[1:-1]     # get rid of task name & closing brace line
    lines = [ln.strip(' ') for ln in lines]     # strip indents
    lines = [ln.strip(' ,') for ln in lines]        # strip indents & commas
    lines = [ln for ln in lines if not ln == '']    # remove empty lines
    return lines

def simplify_file(text: str) -> list[str]:
    lines = text.split('\n')                    # split into lines
    lines = [ln.strip() for ln in lines]        # strip indents
    lines = [ln for ln in lines if ln != '']    # remove empty lines
    return lines

def _gen_call_lines_local(wf: Any, step: Any) -> list[str]:
    reset_globals()
    do_preprocessing_workflow(wf)
    processes = nextflow.generate.process.generate_processes(wf)
    workflows = nextflow.generate.workflow.generate_workflows(wf, processes)
    vmanager = init_variable_manager_for_task(wf)
    update_variables(wf, vmanager)
    
    if step.tool.id() in processes:
        task = processes[step.tool.id()]
    elif step.tool.id() in workflows:
        task = workflows[step.tool.id()]
    else:
        raise RuntimeError(f"Could not find task for step {step.id()}")

    call = nextflow.generate.workflow.call.gen_task_call(
        alias=step.id(),
        task=task,
        vmanager=vmanager,
        step=step
    )
    return simplify_call(call)




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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int2',
            'val in_int3',
            'val in_str2',
            'val in_str3',
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
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int1',
            'val in_int2',
            'val in_int3',
            'val in_int4',
            'val in_str1',
            'val in_str2',
            'val in_str3',
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_str1',
            'val in_str2',
            'val in_str3',
            'val in_int2',
            'val in_int4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)

        # sub wf process
        step = wf.step_nodes["stp2"]
        scope = nextflow.Scope()
        scope.update(step)
        step = step.tool.step_nodes["stp1"]
        scope.update(step)
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_str1',
            'val in_str2',
            'val in_str3',
            'val in_int2',
            'val in_int4',
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
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file, stageAs: \'in_file\'',
            'val in_int2',
            'val in_int3',
            'val in_int4',
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
        actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
        actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        
        expected_task_inputs = {
            'inFile',
            'inStr1',
            'inStr2',
            'inStr3',
            'inInt2',
            'inInt4',
        }
        expected_param_inputs = set()
        expected_static_inputs = {'inInt1'}
        expected_ignored_inputs = {
            'inStr4',
            'inInt3',
        }
        
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)
        
        # SubMinimalTaskInputsTestWF
        tool = wf.step_nodes['stp2'].tool

        actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
        actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
        actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
        actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        
        expected_task_inputs = {'inFile'}
        expected_param_inputs = {'inStr1', 'inInt2'}
        expected_static_inputs = set()
        expected_ignored_inputs = {
            'inStr2',
            'inStr3',
            'inInt1',
            'inInt3',
        }
        
        self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)





class TestVariableManager(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    @unittest.skip('implement in future')
    def test1(self) -> None:
        raise NotImplementedError





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




class TestFiles(unittest.TestCase):

    def setUp(self) -> None:
        self.maxDiff = None
        reset_globals()

    def test_files_created(self) -> None:
        wf = SubworkflowTestWF()
        _, subtask_dict = translator.translate_workflow_internal(wf)
        
        # check correct number of subtasks created
        self.assertEqual(len(subtask_dict), 6)
        expected_filepaths = set([
            'modules/file_test_tool',
            'modules/string_test_tool',
            'modules/int_test_tool',
            'modules/string_opt_test_tool',
            'subworkflows/oranges_workflow',
            'subworkflows/apples_workflow',
        ])
        actual_filepaths = set(subtask_dict.keys())
        self.assertSetEqual(actual_filepaths, expected_filepaths)

    def test_main_workflow_format(self) -> None:
        wf = AssemblyTestWF()
        mainstr, _ = translator.translate_workflow_internal(wf)
        actual = simplify_file(mainstr)
        expected = [
            "nextflow.enable.dsl=2",
            "include { FASTQC as FASTQC1 } from './modules/fastqc'",
            "include { FASTQC as FASTQC2 } from './modules/fastqc'",
            "include { FASTQC as FASTQC3 } from './modules/fastqc'",
            "include { CAT_TEST_TOOL } from './modules/cat_test_tool'",
            "include { UNICYCLER } from './modules/unicycler'",
            "// data which will be passed as channels",
            "ch_in_forward_reads  = Channel.fromPath( params.in_forward_reads )",
            "ch_in_long_reads     = Channel.fromPath( params.in_long_reads )",
            "ch_in_reverse_reads  = Channel.fromPath( params.in_reverse_reads )",
            "ch_test_input        = Channel.fromPath( params.test_input )",
            "// data which will be passed as optional files",
            "fastqc1_adapters      = file( params.fastqc1_adapters )",
            "fastqc1_contaminants  = file( params.fastqc1_contaminants )",
            "fastqc1_limits        = file( params.fastqc1_limits )",
            "fastqc2_adapters      = file( params.fastqc2_adapters )",
            "fastqc2_contaminants  = file( params.fastqc2_contaminants )",
            "fastqc2_limits        = file( params.fastqc2_limits )",
            "workflow {",
            "FASTQC1(",
            "ch_in_forward_reads,",
            "fastqc1_adapters,",
            "fastqc1_contaminants,",
            "fastqc1_limits",
            ")",
            "FASTQC2(",
            "ch_in_reverse_reads,",
            "fastqc2_adapters,",
            "fastqc2_contaminants,",
            "fastqc2_limits",
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
        self.assertEqual(len(actual), len(expected))
        for ln in actual:
            self.assertIn(ln, expected)

    def test_process_format(self) -> None:
        wf = AssemblyTestWF()
        _, substr_dict = translator.translate_workflow_internal(wf)
        actual_lines = simplify_file(substr_dict['modules/fastqc'])
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
            'def adapters = adapters != \'NO_FILE\' ? "--adapters ${adapters}" : ""',
            'def contaminants = contaminants != \'NO_FILE\' ? "--contaminants ${contaminants}" : ""',
            'def limits = limits != \'NO_FILE\' ? "--limits ${limits}" : ""',
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
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_subworkflow_format(self) -> None:
        wf = SubworkflowTestWF()
        _, substr_dict = translator.translate_workflow_internal(wf)
        actual_lines = simplify_file(substr_dict['subworkflows/apples_workflow'])
        expected_lines = [
            "nextflow.enable.dsl=2",
            "include { STRING_TEST_TOOL as STRING_TOOL } from '../modules/string_test_tool'",
            "include { STRING_OPT_TEST_TOOL as STRING_OPT_TOOL } from '../modules/string_opt_test_tool'",
            "include { ORANGES_WORKFLOW as ORANGES_SUBWORKFLOW } from './oranges_workflow'",
            "workflow APPLES_WORKFLOW {",
            "main:",
            "STRING_TOOL()",
            "STRING_OPT_TOOL()",
            "ORANGES_SUBWORKFLOW(",
            "STRING_TOOL.out.out",
            ")",
            "emit:",
            "outIntFile = ORANGES_SUBWORKFLOW.out.out",
            "outStringFile = STRING_TOOL.out.out",
            "}",
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_config_format(self) -> None:
        # TODO expand this to subworkflow with subworkflow, process specific params
        wf = SubworkflowTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        actual_lines = simplify_file(config)
        expected_lines = [
            "docker.enabled = true",
            "params {",
            "// OUTPUT DIRECTORY",
            "outdir  = './outputs'",
            "// INPUTS (MANDATORY)",
            "in_file  = null  // (generic file)",
            "in_int   = null  // (integer)",
            "in_str   = null  // (string)",
            "// INPUTS (OPTIONAL)",
            "in_str_opt  = null  // (optional string)",
            "// PROCESS: INT_TEST_TOOL",
            "int_test_tool.inp  = null  // (integer)",
            "// PROCESS: STRING_OPT_TEST_TOOL",
            "string_opt_test_tool.inp  = null  // (optional string)",
            "// PROCESS: STRING_TEST_TOOL",
            "string_test_tool.inp  = null  // (string)",
            "// SUBWORKFLOW: APPLES_WORKFLOW",
            "apples_workflow.in_int      = null  // (integer)",
            "apples_workflow.in_str      = null  // (string)",
            "apples_workflow.in_str_opt  = null  // (optional string)",
            "// SUBWORKFLOW: ORANGES_WORKFLOW",
            "oranges_workflow.in_int  = null  // (integer)",
            "}",
        ]
        for ln in actual_lines:
            print(ln)
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_config_params(self) -> None:
        # test_nonfile
        # string, int, bool
        wf = AllInputTypesTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        actual_lines = simplify_file(config)
        expected_lines = [
            "docker.enabled = true",
            "params {",
            "// OUTPUT DIRECTORY",
            "outdir  = './outputs'",
            "// INPUTS (MANDATORY)",
            "in_file               = null  // (generic file)",
            "in_file_array         = []    // (array)         eg. [file1, ...]",
            "in_secondaries        = []    // (indexedbam)    eg. [bam, bai]",
            "in_secondaries_array  = [[]]  // (array)         eg. [[bam, bai]]",
            "in_filepair           = []    // (fastqpair)     eg. [pair1, pair2]",
            "in_filepair_array     = [[]]  // (array)         eg. [[pair1, pair2]]",
            "in_nonfile            = null  // (integer)",
            "in_nonfile_array      = []    // (array)         eg. [integer1, ...]",
            "// INPUTS (OPTIONAL)",
            "in_file_array_optional         = ['NO_FILE']",
            "in_file_optional               = 'NO_FILE'",
            "in_secondaries_array_optional  = [['NO_FILE', 'NO_FILE']]",
            "in_secondaries_optional        = ['NO_FILE', 'NO_FILE']",
            "in_filepair_array_optional     = [['NO_FILE', 'NO_FILE']]",
            "in_filepair_optional           = ['NO_FILE', 'NO_FILE']",
            "in_nonfile_array_optional      = []    // (optional array)    eg. [integer1, ...]",
            "in_nonfile_optional            = null  // (optional integer)",
            "// PROCESS: NON_FILE_TEST_TOOL",
            "non_file_test_tool.nonfile                 = null  // (integer)",
            "non_file_test_tool.nonfile_array           = []    // (array)             eg. [integer1, ...]",
            "non_file_test_tool.nonfile_array_optional  = []    // (optional array)    eg. [integer1, ...]",
            "non_file_test_tool.nonfile_optional        = null  // (optional integer)",
            "}",
        ]
        for ln in actual_lines:
            print(ln)
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_config_params_pythontool(self) -> None:
        # test_nonfile
        # string, int, bool
        wf = InputsPythonToolTestWF()
        do_preprocessing_workflow(wf)
        config = translator.stringify_translated_inputs({})
        actual_lines = simplify_file(config)
        expected_lines = [
            "docker.enabled = true",
            "params {",
            "// OUTPUT DIRECTORY",
            "outdir  = './outputs'",
            "// INPUTS (MANDATORY)",
            "in_file            = null  // (generic file)",
            "in_secondary_type  = []    // (generic file)  eg. [txt, txt]",
            "in_int             = null  // (integer)",
            "in_str             = null  // (string)",
            "in_str_arr         = []    // (array)         eg. [string1, ...]",
            "// PROCESS: JOIN_ARRAY_PYTHON_TEST_TOOL",
            "join_array_python_test_tool.code_file  = '/home/grace/work/pp/translation/janis-core/templates/JoinArrayPythonTestTool.py'",
            "join_array_python_test_tool.inp        = []  // (array)  eg. [string1, ...]",
            "// PROCESS: MULTI_TYPES_INPUT_PYTHON_TOOL",
            "multi_types_input_python_tool.code_file  = '/home/grace/work/pp/translation/janis-core/templates/MultiTypesInputPythonTool.py'",
            "multi_types_input_python_tool.inp2       = null  // (string)",
            "multi_types_input_python_tool.inp3       = null  // (integer)",
            "// PROCESS: SECONDARY_INPUT_PYTHON_TEST_TOOL",
            "secondary_input_python_test_tool.code_file  = '/home/grace/work/pp/translation/janis-core/templates/SecondaryInputPythonTestTool.py'",
            "}",
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

    def test_basic(self) -> None:
        wf = AssemblyTestWF()
        maintask, subtask_dict = translator.translate_workflow_internal(wf)
        actual_subtasks = len(subtask_dict)
        expected_subtasks = 3
        self.assertEqual(actual_subtasks, expected_subtasks)

        actual_paths = set(subtask_dict.keys())
        expected_paths = set([
            'modules/fastqc',
            'modules/unicycler',
            'modules/cat_test_tool',
        ])
        self.assertSetEqual(actual_paths, expected_paths)
        print(maintask)
        print()
        subtask = subtask_dict['modules/fastqc']
        print(subtask)
        print()
        subtask = subtask_dict['modules/unicycler']
        print(subtask)
        print()
        subtask = subtask_dict['modules/cat_test_tool']
        print(subtask)
        print()
        print()

    def test_imports(self) -> None:
        # translate workflow, building all nf items and files
        # TODO improve this - more complex workflow imports
        wf = SubworkflowTestWF()
        _, substr_dict = translator.translate_workflow_internal(wf)
        actual_lines = simplify_file(substr_dict['subworkflows/apples_workflow'])
        actual_lines = [ln for ln in actual_lines if ln.startswith('include {')]
        expected_lines = [
            "include { STRING_TEST_TOOL as STRING_TOOL } from '../modules/string_test_tool'",
            "include { STRING_OPT_TEST_TOOL as STRING_OPT_TOOL } from '../modules/string_opt_test_tool'",
            "include { ORANGES_WORKFLOW as ORANGES_SUBWORKFLOW } from './oranges_workflow'"
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_channel_declarations(self) -> None:
        """
        Every Optional(File) type wf input should have a channel. 
        '.ifEmpty(null)' should appear in the channel string definition.
        """
        wf = AllInputTypesTestWF()
        mainstr, _ = translator.translate_workflow_internal(wf)
        actual_lines = simplify_file(mainstr)
        actual_lines = [ln for ln in actual_lines if 'Channel.' in ln]
        expected_lines = [
            'ch_in_file               = Channel.fromPath( params.in_file )',
            'ch_in_file_array         = Channel.fromPath( params.in_file_array ).toList()',
            'ch_in_secondaries        = Channel.fromPath( params.in_secondaries ).toList()',
            'ch_in_secondaries_array  = Channel.fromPath( params.in_secondaries_array.flatten() ).collate( 2 )',
            'ch_in_filepair           = Channel.of( params.in_filepair ).toList()',
            'ch_in_filepair_array     = Channel.of( params.in_filepair_array ).toList()',
            'ch_in_nonfile            = Channel.of( params.in_nonfile )',
            'ch_in_nonfile_array      = Channel.of( params.in_nonfile_array ).toList()',
        ]
        print(mainstr)
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

    def test_variable_declarations(self) -> None:
        wf = AllInputTypesTestWF()
        mainstr, _ = translator.translate_workflow_internal(wf)
        actual_lines = simplify_file(mainstr)
        actual_lines = [ln for ln in actual_lines if 'file(' in ln]
        expected_lines = [
            'in_file_optional               = file( params.in_file_optional )',
            'in_filepair_optional           = params.in_filepair_optional.each { value -> file(value) }',
            'in_secondaries_optional        = params.in_secondaries_optional.each { value -> file(value) }',
            'in_file_array_optional         = params.in_file_array_optional.each { value -> file(value) }',
            'in_filepair_array_optional     = params.in_filepair_array_optional.each { outer -> outer.each { inner -> file(inner) } }',
            'in_secondaries_array_optional  = params.in_secondaries_array_optional.each { outer -> outer.each { inner -> file(inner) } }',
        ]
        print(mainstr)
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

    def test_duplicate_tool_usage(self) -> None:
        wf = DuplicateTasksTestWF()
        maintask, subtask_dict = translator.translate_workflow_internal(wf)

        # main workflow fmt
        print(maintask)
        actual_lines = simplify_file(maintask)
        expected_lines = [
            "nextflow.enable.dsl=2",
            "include { ECHO_TEST_TOOL as STP1 } from './modules/echo_test_tool'",
            "include { ECHO_TEST_TOOL as STP2 } from './modules/echo_test_tool'",
            "include { ECHO_TEST_WORKFLOW1 as STP3 } from './subworkflows/echo_test_workflow1'",
            "include { ECHO_TEST_WORKFLOW1 as STP4 } from './subworkflows/echo_test_workflow1'",
            "include { ECHO_TEST_WORKFLOW2 as STP5 } from './subworkflows/echo_test_workflow2'",
            "// data which will be passed as channels",
            "ch_in_file  = Channel.fromPath( params.in_file )",
            "workflow {",
            "STP1(",
            "ch_in_file,",
            "5,",
            "'hello'",
            ")",
            "STP2(",
            "ch_in_file,",
            "10,",
            "'there'",
            ")",
            "STP3(",
            "ch_in_file,",
            "15,",
            "'my'",
            ")",
            "STP4(",
            "ch_in_file,",
            "20,",
            "'friend'",
            ")",
            "STP5(",
            "ch_in_file",
            ")",
            "}",
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

        # echo test tool fmt
        subtask = subtask_dict['modules/echo_test_tool']
        print(subtask)
        actual_lines = simplify_file(subtask)
        expected_lines = [
            "nextflow.enable.dsl=2",
            "process ECHO_TEST_TOOL {",
            "debug true",
            'container "ubuntu:latest"',
            'publishDir "${params.outdir}/echo_test_tool"',
            "input:",
            "path in_file, stageAs: 'in_file'",
            "val in_int",
            "val in_str",
            "output:",
            "stdout, emit: out",
            "script:",
            '"""',
            "echo \\",
            "${in_file} \\",
            "${in_str} \\",
            "${in_int} \\",
            '"""',
            "}",
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

        # echo_test_workflow1 fmt
        subtask = subtask_dict['subworkflows/echo_test_workflow1']
        print(subtask)
        actual_lines = simplify_file(subtask)
        expected_lines = [
            "nextflow.enable.dsl=2",
            "include { ECHO_TEST_TOOL as STP1 } from '../modules/echo_test_tool'",
            "include { ECHO_TEST_TOOL as STP2 } from '../modules/echo_test_tool'",
            "workflow ECHO_TEST_WORKFLOW1 {",
            "take:",
            "ch_in_file",
            "ch_in_int",
            "ch_in_str",
            "main:",
            "STP1(",
            "ch_in_file,",
            "25,",
            "'good'",
            ")",
            "STP2(",
            "ch_in_file,",
            "30,",
            "'night'",
            ")",
            "emit:",
            "out = STP2.out.out",
            "}",
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

        # echo_test_workflow2 fmt
        subtask = subtask_dict['subworkflows/echo_test_workflow2']
        print(subtask)
        actual_lines = simplify_file(subtask)
        expected_lines = [
            "nextflow.enable.dsl=2",
            "include { ECHO_TEST_WORKFLOW1 as LOL } from './echo_test_workflow1'",
            "workflow ECHO_TEST_WORKFLOW2 {",
            "take:",
            "ch_in_file",
            "main:",
            "LOL(",
            "ch_in_file,",
            "25,",
            "'good'",
            ")",
            "}",
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

    def test_duplicate_subworkflow_usage(self) -> None:
        raise NotImplementedError
    



class TestTranslateInterface(unittest.TestCase):
    """
    Tests janis CommandTool can be parsed to nextflow process (end-to-end).
    """

    def setUp(self) -> None:
        reset_globals() 
    
    def test_fastqc_tool(self) -> None:
        tool = FastqcTestTool()
        process = translate(tool, 'nextflow')
        print()
    
    def test_bwamem_tool(self) -> None:
        tool = BwaMemTestTool()
        process = translate(tool, 'nextflow')
        print()
    
    def test_gridss_tool(self) -> None:
        tool = GridssTestTool()
        process = translate(tool, 'nextflow')
        print()
 
    def test_assembly_workflow(self) -> None:
        wf = AssemblyTestWF()
        maintask, inputs_dict, subtask_dict = translate(wf, 'nextflow')
        print()

class TestCmdtoolProcessDirectives(unittest.TestCase):
    """
    Tests identifying tool inputs which should be process inputs.
    Need a process input for each tool input in step sources

    INCLUDES WORKFLOW OUTPUTS as publishDir
    """

    def setUp(self) -> None:
        reset_globals()

    def test_directives(self) -> None:
        wf = DirectivesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path reads'
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
        
        # filepair array
        step = wf.step_nodes["stp3"]
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        actual_inputs = {inp.get_string() for inp in process.inputs}
        expected_inputs = {
            "path inp1, stageAs: 'inp1'",
            "val inp2",
        }
        self.assertEqual(actual_inputs, expected_inputs)  

        step = wf.step_nodes["stp2"]
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'stdout, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)  

    def test_wildcard(self) -> None:
        step = self.wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "myfile.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_wildcard_array(self) -> None:
        wf = WildcardSelectorOutputTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp2"]
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "*.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_input_selector(self) -> None:
        wf = InputSelectorTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path inp, emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_input_selector_param(self) -> None:
        step = self.wf.step_nodes["stp4"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "myfile.txt", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
        
    def test_input_selector_array(self) -> None:
        wf = InputSelectorTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp3']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path inp, emit: out'}
        print(process.get_string())
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_filenames_generated(self) -> None:
        wf = FilenameTestWF2()
        do_preprocessing_workflow(wf)
        
        # inp1 is in step.sources
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        raise NotImplementedError
        # eg read1.fastq, read2.fastq
        # collection method is list, len(list) == 2.
        step = self.wf.step_nodes['stp6']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "[${inp.simpleName + "-R1.fastq"}, ${inp.simpleName + "-R2.fastq"}]", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_secondaries(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'tuple path("*.bam"), path("*.bam.bai"), emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_secondaries_replaced(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp3']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'tuple path("*.bam"), path("*.bai"), emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_secondaries_edge_basename(self) -> None:
        wf = SecondariesTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp5']
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "${inp.simpleName + ".gz"}", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_edge_markduplicates_metrics(self) -> None:
        step = self.wf.step_nodes['stp8']
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
        
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
        process = nextflow.generate.process.generate_process(step.tool)
        
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
        
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
        process = nextflow.generate.process.generate_process(step.tool)
        
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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
        process = nextflow.generate.process.generate_process(step.tool)
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




class TestTranslateWorkflowInternal(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()

    def test_assembly(self) -> None:
        wf = AssemblyTestWF()
        maintask, _ = translator.translate_workflow_internal(wf)




class TestTranslateToolInternal(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()

    def test_fastqc_tool(self) -> None:
        tool = FastqcTestTool()
        toolstr = translator.translate_tool_internal(tool)
    
    def test_bwamem_tool(self) -> None:
        tool = BwaMemTestTool()
        toolstr = translator.translate_tool_internal(tool)
    
    def test_gridss_tool(self) -> None:
        tool = GridssTestTool()
        toolstr = translator.translate_tool_internal(tool)




class TestTranslateCodeToolInternal(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()

    def test1(self) -> None:
        tool = FileOutputPythonTestTool()
        toolstr = translator.translate_code_tool_internal(tool)
        print(toolstr)
        print()
    
    def test2(self) -> None:
        tool = MultiTypesInputPythonTool()
        toolstr = translator.translate_code_tool_internal(tool)
        print(toolstr)
        print()
    
    def test3(self) -> None:
        tool = SecondaryInputPythonTestTool()
        toolstr = translator.translate_code_tool_internal(tool)
        print(toolstr)
        print()
        



class TestTranslateHelperFiles(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_python_tool_helpers(self) -> None:
        # first test wf
        wf = InputsPythonToolTestWF()
        helper_files = translator.translate_helper_files(wf)
        actual_paths = list(helper_files.keys())
        expected_paths = [
            'templates/MultiTypesInputPythonTool.py',
            'templates/JoinArrayPythonTestTool.py',
            'templates/SecondaryInputPythonTestTool.py'
        ]
        for filepath, filecontents in helper_files.items():
            print(f'\n--- {filepath} ---')
            print(filecontents)
        self.assertEqual(len(actual_paths), len(expected_paths))
        for path in actual_paths:
            self.assertIn(path, expected_paths)
        
        # second test wf
        wf = OutputsPythonToolTestWF()
        helper_files = translator.translate_helper_files(wf)
        actual_paths = list(helper_files.keys())
        expected_paths = [
            'templates/FileOutputPythonTestTool.py',
            'templates/FileInputPythonTestTool.py',
            'templates/SplitTextPythonTestTool.py'
        ]
        self.assertEqual(len(actual_paths), len(expected_paths))
        for path in actual_paths:
            self.assertIn(path, expected_paths)
 
    def test_template_helpers(self) -> None:
        wf = FilesDirectoriesToCreateTestWF()
        helper_files = translator.translate_helper_files(wf)
        actual_paths = list(helper_files.keys())
        expected_paths = [
            'templates/myscript.sh',
            'templates/check_value.js',
            'templates/MultiTypesInputPythonTool.py'
        ]
        self.assertEqual(len(actual_paths), len(expected_paths))
        for path in actual_paths:
            self.assertIn(path, expected_paths)
        



class TestSubWorkflows(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()







class TestPythontoolProcessInputs(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()
        self.wf = InputsPythonToolTestWF()
        do_preprocessing_workflow(self.wf)
        self.processes = nextflow.generate.process.generate_processes(self.wf)
        self.workflows = nextflow.generate.workflow.generate_workflows(self.wf, self.processes)
        self.vmanager = init_variable_manager_for_task(self.wf)
        update_variables(self.wf, self.vmanager)

    def test_input_generation(self) -> None:
        # File, String, Int input types
        task = self.processes['MultiTypesInputPythonTool']
        actual_inputs = [inp.get_string() for inp in task.ordered_inputs]
        expected_inputs = {
            "path code_file",
            "path inp1, stageAs: 'inp1'"
        }
        self.assertEqual(len(actual_inputs), len(expected_inputs))
        for ln in expected_inputs:
            self.assertIn(ln, actual_inputs)
        
        # Array(String) input type
        task = self.processes['JoinArrayPythonTestTool']
        actual_inputs = [inp.get_string() for inp in task.ordered_inputs]
        expected_inputs = {
            "path code_file",
        }
        self.assertEqual(len(actual_inputs), len(expected_inputs))
        for ln in expected_inputs:
            self.assertIn(ln, actual_inputs)
        
        # File (secondaries) input type
        task = self.processes['SecondaryInputPythonTestTool']
        actual_inputs = {inp.get_string() for inp in task.ordered_inputs}
        expected_inputs = {
            "path code_file",
            'tuple path(bam), path(bai)',
        }
        self.assertEqual(len(actual_inputs), len(expected_inputs))
        for ln in expected_inputs:
            self.assertIn(ln, actual_inputs)




class TestPythontoolProcessOutputs(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()
        self.wf = OutputsPythonToolTestWF()
        do_preprocessing_workflow(self.wf)
        self.processes = nextflow.generate.process.generate_processes(self.wf)
        self.workflows = nextflow.generate.workflow.generate_workflows(self.wf, self.processes)
        self.vmanager = init_variable_manager_for_task(self.wf)
        update_variables(self.wf, self.vmanager)
    
    def test_output_generation(self) -> None:
        # file output
        task = self.processes['FileOutputPythonTestTool']
        actual_outputs = [out.get_string() for out in task.ordered_outputs]
        expected_outputs = [
            'val "${file("${task.workDir}/" + file("${task.workDir}/out_out").text.replace(\'"\', \'\'))}", emit: out'
        ]
        for ln in expected_outputs:
            self.assertIn(ln, actual_outputs)
        
        # String output
        task = self.processes['FileInputPythonTestTool']
        actual_outputs = [out.get_string() for out in task.ordered_outputs]
        expected_outputs = {
            'val "${file("${task.workDir}/out_out").text}", emit: out'
        }
        for ln in expected_outputs:
            self.assertIn(ln, actual_outputs)
        
        # Array(String) output
        task = self.processes['SplitTextPythonTestTool']
        actual_outputs = [out.get_string() for out in task.ordered_outputs]
        expected_outputs = {
            'val "${file("${task.workDir}/out_out").text.replace(\'[\', \'\').replace(\']\', \'\')}", emit: out'
        }
        for ln in expected_outputs:
            self.assertIn(ln, actual_outputs)



class TestPythontoolProcess(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()
        self.wf = InputsPythonToolTestWF()
        do_preprocessing_workflow(self.wf)
        self.processes = nextflow.generate.process.generate_processes(self.wf)
        self.workflows = nextflow.generate.workflow.generate_workflows(self.wf, self.processes)
        self.vmanager = init_variable_manager_for_task(self.wf)
        update_variables(self.wf, self.vmanager)

    def test_format(self) -> None:
        task = self.processes['MultiTypesInputPythonTool']
        actual_lines = simplify_file(task.get_string())
        expected_lines = [
            'process MULTI_TYPES_INPUT_PYTHON_TOOL {',
            'debug true',
            'container "python:3.8.1"',
            'publishDir "${params.outdir}/multi_types_input_python_tool"',
            'input:',
            'path code_file',
            'path inp1, stageAs: \'inp1\'',
            'output:',
            'val "${file("${task.workDir}/" + file("${task.workDir}/out_out").text.replace(\'"\', \'\'))}", emit: out',
            'exec:',
            'script:',
            '"""',
            '#!/usr/bin/env python',
            'from ${code_file} import code_block',
            'import os',
            'import json',
            'result = code_block(',
            'inp1="${inp1}",',
            'inp2="${params.multi_types_input_python_tool.inp2}",',
            'inp3=${params.multi_types_input_python_tool.inp3}',
            ')',
            'work_dir = os.getcwd()',
            'for key in result:',
            'with open(os.path.join(work_dir, f"out_{key}"), "w") as fp:',
            'fp.write(json.dumps(result[key]))',
            '"""',
            '}',
        ]
        self.assertEqual(actual_lines, expected_lines)




class TestPythontoolProcessScript(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()
        self.wf = InputsPythonToolTestWF()
        do_preprocessing_workflow(self.wf)
        self.processes = nextflow.generate.process.generate_processes(self.wf)
        self.workflows = nextflow.generate.workflow.generate_workflows(self.wf, self.processes)
        self.vmanager = init_variable_manager_for_task(self.wf)
        update_variables(self.wf, self.vmanager)

    def test_input_references(self) -> None:
        # File, String, Int input types
        task = self.processes['MultiTypesInputPythonTool']
        actual_script = task.script
        expected_lines = {
            'result = code_block(',
            'inp1="${inp1}",',
            'inp2="${params.multi_types_input_python_tool.inp2}",',
            'inp3=${params.multi_types_input_python_tool.inp3}',
            ')'
        }
        print(task.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

        # Array(String) input type
        task = self.processes['JoinArrayPythonTestTool']
        actual_script = task.script
        expected_lines = {
            'result = code_block(inp="${params.join_array_python_test_tool.inp}".split(" "))',
        }
        print(task.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

        # File (secondaries) input type
        task = self.processes['SecondaryInputPythonTestTool']
        actual_script = task.script
        expected_lines = {
            'result = code_block(inp="${bam}")',
        }
        print(task.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)




class TestEntityTracing(unittest.TestCase):
    """
    This test group checks we can handle array / single datatype mismatches
    in a workflow. This occurs in wgsgermline. 
    """
    def setUp(self) -> None:
        reset_globals()
        self.maxDiff = None
        self.wf = EntityTraceTestWF()
        do_preprocessing_workflow(self.wf)

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
        self.assertEqual(get_array_depth(String()), 0)
        self.assertEqual(get_array_depth(BamBai()), 0)
        self.assertEqual(get_array_depth(FastqGzPair()), 0)
        self.assertEqual(get_array_depth(Array(String())), 1)
        self.assertEqual(get_array_depth(Array(BamBai())), 1)
        self.assertEqual(get_array_depth(Array(FastqGzPair())), 1)
        self.assertEqual(get_array_depth(Array(Array(String()))), 2)
        self.assertEqual(get_array_depth(Array(Array(BamBai()))), 2)
        self.assertEqual(get_array_depth(Array(Array(FastqGzPair()))), 2)
        self.assertEqual(get_array_depth(Array(Array(Array(String())))), 3)
        self.assertEqual(get_array_depth(Array(Array(Array(BamBai())))), 3)
        self.assertEqual(get_array_depth(Array(Array(Array(FastqGzPair())))), 3)
    
    def test_is_array_depth_mismatch(self):
        # no mismatch
        self.assertFalse(is_array_depth_mismatch(String(), String(), False))
        self.assertFalse(is_array_depth_mismatch(Array(String()), Array(String()), False))
        self.assertFalse(is_array_depth_mismatch(Array(Array(String())), Array(Array(String())), False))
        # mismatch
        self.assertTrue(is_array_depth_mismatch(String(), Array(String()), False))
        self.assertTrue(is_array_depth_mismatch(Array(String()), String(), False))
        self.assertTrue(is_array_depth_mismatch(Array(Array(String())), Array(String()), False))

    # identifying common types ---
    def test_get_common_type_1(self):
        # non-union types, has intersection
        srctype = Bam()
        desttype = Bam()
        common_type = get_common_type(srctype, desttype)
        self.assertEqual(type(common_type), Bam)
        
    def test_get_common_type_2(self):
        # non-union types, no intersection
        srctype = Fasta()
        desttype = Bam()
        common_type = get_common_type(srctype, desttype)
        self.assertIsNone(common_type)
        
    def test_get_common_type_3(self):
        # single union type, has intersection
        srctype = UnionType(Bam, Sam, Cram)
        desttype = Bam()
        common_type = get_common_type(srctype, desttype)
        self.assertEqual(type(common_type), Bam)
        
    def test_get_common_type_4(self):
        # single union type, no intersection
        srctype = UnionType(Bam, Sam, Cram)
        desttype = Fasta()
        common_type = get_common_type(srctype, desttype)
        self.assertIsNone(common_type)
        
    def test_get_common_type_5(self):
        # both union types, has intersection
        srctype = UnionType(Bam, Sam, Cram)
        desttype = UnionType(Sam, Cram)
        common_type = get_common_type(srctype, desttype)
        self.assertEqual(type(common_type), Sam)
        
    def test_get_common_type_6(self):
        # both union types, no intersection
        srctype = UnionType(Bam, Sam, Cram)
        desttype = UnionType(Fasta, Fastq)
        common_type = get_common_type(srctype, desttype)
        self.assertIsNone(common_type)
    

    # array datatype mismatches ---
    # caused by datatypes involving arrays and/or scatter relationships
    def test_is_datatype_mismatch_arrays(self):
        """tests nfgen.plumibng.is_datatype_mismatch(srctype, desttype, srcscatter, destscatter)"""
        # secondary array (always considered mismatch)
        self.assertTrue(is_datatype_mismatch(Array(BamBai()), Array(BamBai()), False))

        # array depth - single single 
        self.assertFalse(is_datatype_mismatch(Bam(), Bam(), False)) # no mismatch
        self.assertTrue(is_datatype_mismatch(Bam(), Array(Bam()), False)) # mismatch
        
        # array depth - secondary secondary
        self.assertFalse(is_datatype_mismatch(BamBai(), BamBai(), False)) # no mismatch
        self.assertTrue(is_datatype_mismatch(Array(BamBai()), BamBai(), False)) # mismatch
        
        # array depth - secondary single 
        self.assertTrue(is_datatype_mismatch(BamBai(), Bam(), False)) # mismatch
        self.assertTrue(is_datatype_mismatch(FastaDict(), Bam(), False)) # mismatch
        
        # array depth - file pairs
        self.assertFalse(is_datatype_mismatch(FastqGzPair(), FastqGzPair(), False)) # no mismatch
        self.assertTrue(is_datatype_mismatch(FastqGzPair(), Array(FastqGzPair()), False)) # mismatch
        
    def test_is_datatype_mismatch_basetype(self):
        # single single
        self.assertFalse(is_datatype_mismatch(Bam(), Bam(), False)) # no mismatch
        self.assertTrue(is_datatype_mismatch(Bam(), Fasta(), False)) # mismatch
        
        # secondary secondary
        self.assertFalse(is_datatype_mismatch(FastaDict(), FastaDict(), False)) # no mismatch
        self.assertTrue(is_datatype_mismatch(FastaDict(), BamBai(), False)) # mismatch
        
        # secondary single
        self.assertFalse(is_datatype_mismatch(Bam(), Bam(), False)) # no mismatch
        self.assertTrue(is_datatype_mismatch(Bam(), BamBai(), False)) # mismatch
        
        # file pairs
        self.assertFalse(is_datatype_mismatch(FastqGzPair(), FastqGzPair(), False)) # no mismatch
        self.assertTrue(is_datatype_mismatch(FastqGzPair(), BamBai(), False)) # mismatch

    def test_generate_datatype_mismatch_plumbing(self):
        # plumbing not required ---
        srctype, desttype, destscatter = Bam(), Bam(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '')
        
        srctype, desttype, destscatter = FastaWithIndexes(), FastaWithIndexes(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '')
        
        srctype, desttype, destscatter = FastqGzPair(), FastqGzPair(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '')

        srctype, desttype, destscatter = Array(Bam()), Array(Bam()), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '')
        
        # plumbing required ---
        
        # BASETYPES
        # secondary, secondary
        srctype, desttype, destscatter = FastaWithIndexes(), FastaDict(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.map{ tuple -> [tuple[0], tuple[4]] }')
        
        # secondary, single
        srctype, desttype, destscatter = FastaWithIndexes(), Fasta(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.map{ tuple -> tuple[0] }')

        # ARRAYS      
        # single, single  
        srctype, desttype, destscatter = Array(Bam()), Bam(), True
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten()')
        
        srctype, desttype, destscatter = Array(Fasta()), Fasta(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().first()')
        
        srctype, desttype, destscatter = Fasta(), Array(Fasta()), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.toList()')

        # any, secondary (always ends with .flatten().toList())
        srctype, desttype, destscatter = FastaWithIndexes(), Array(FastaWithIndexes()), True
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().toList()')
        
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Array(FastaWithIndexes()), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().toList()')
        
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Array(FastaWithIndexes()), True
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().toList()')

        # secondary, secondary 
        srctype, desttype, destscatter = Array(FastaWithIndexes()), FastaWithIndexes(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().collate( 8 ).first()')

        srctype, desttype, destscatter = FastaWithIndexes(), Array(FastaDict()), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.map{ tuple -> [tuple[0], tuple[4]] }.flatten().toList()')
        
        # secondary, single 
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Fasta(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.map{ tuple -> tuple[0] }.flatten().first()')

        # file pairs
        srctype, desttype, destscatter = Array(FastqGzPair()), FastqGzPair(), True
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().collate( 2 )')
        
        srctype, desttype, destscatter = Array(FastqGzPair()), FastqGzPair(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter)
        self.assertEqual(plumbing, '.flatten().collate( 2 ).first()')

        
        


class TestPlumbingTypeMismatch(unittest.TestCase):
    """
    This test group checks we can handle array / single datatype mismatches
    in a workflow. This occurs in wgsgermline. 
    """
    def setUp(self) -> None:
        reset_globals()
        self.wf = PlumbingTypeMismatchTestWF()

    def test_secondary_single_mismatch(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['bambai_to_bam'])
        expected = ['ch_in_bam_bai.map{ tuple -> tuple[0] }']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_single_array_mismatch(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['bambai_to_bam_array'])
        expected = ['ch_in_bam_bai.map{ tuple -> tuple[0] }.toList()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_secondary_mismatch(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['fastawithindexes_to_fastadict'])
        expected = ['ch_in_fasta_with_indexes.map{ tuple -> [tuple[0], tuple[4]] }']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_secondary_array_mismatch(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['fastawithindexes_to_fastadict_array'])
        expected = ['ch_in_fasta_with_indexes.map{ tuple -> [tuple[0], tuple[4]] }.flatten().toList()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_array_to_single(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['array_to_single'])
        expected = ['ch_in_file_array.flatten().first()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_single_to_array(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['single_to_array'])
        expected = ['ch_in_file.toList()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_array_to_secondary(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['secondary_array_to_secondary'])
        expected = ['ch_in_bam_bai_array.flatten().collate( 2 ).first()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_to_secondary_array(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['secondary_array_to_secondary_array'])
        expected = ['ch_in_bam_bai_array.flatten().toList()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)




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
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['scatter_to_scatter'])
        expected = ['PRESTEP1.out.out']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_scatter_to_array(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['scatter_to_array'])
        expected = ['PRESTEP1.out.out.toList()']
        self.assertEqual(actual, expected)

    def test_array_to_scatter1(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['prestep1'])
        expected = ['ch_in_file_array.flatten()']
        self.assertEqual(actual, expected)
    
    def test_array_to_scatter2(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['array_to_scatter'])
        expected = ['PRESTEP2.out.out.flatten()']
        self.assertEqual(actual, expected)

    def test_scatter_secondary_to_scatter_secondary(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['scatter_secondary_to_scatter_secondary'])
        expected = ['PRESTEP3.out.out']
        self.assertEqual(actual, expected)

    def test_scatter_secondary_to_secondary_array(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['scatter_secondary_to_secondary_array'])
        expected = ['PRESTEP3.out.out.flatten().toList()']
        self.assertEqual(actual, expected)

    def test_secondary_array_to_scatter_secondary(self):
        actual = _gen_call_lines_local(self.wf, step=self.wf.step_nodes['secondary_array_to_scatter_secondary'])
        expected = ['ch_in_bam_bai_array.flatten().collate( 2 )']
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
        actual = set(simplify_call(call))

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
    def test_python_tool_inputs(self):
        wf = InputsPythonToolTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp0'])
        expected = [
            "params.multi_types_input_python_tool.code_file",
            "ch_in_file",
        ]
        self.assertListEqual(expected, actual)

    def test_workflow_inputs(self):
        wf = StepInputsWFInputTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp1'])
        expected = [
            "ch_in_file",
            "in_file_opt",
        ]
        self.assertListEqual(expected, actual)
        
    # static step inputs
    def test_static_inputs(self):
        wf = StepInputsStaticTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = ['ch_in_file']
        self.assertListEqual(actual, expected)

    # connections
    def test_connections_files(self) -> None:
        # setup
        wf = StepConnectionsTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = ["STP1.out.out"]
        self.assertListEqual(actual, expected)

    def test_filename_types(self) -> None:
        wf = FilenameTestWF1()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp1'])
        expected = ["ch_in_file", "params.in_str"]
        self.assertListEqual(actual, expected)
        
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = ["ch_in_file"]
        self.assertListEqual(actual, expected)

    def test_subworkflows1(self) -> None:
        wf = DataSourceTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp1'])
        expected = [
            "ch_in_file1",
            "in_file_opt1",
            "params.in_str1",
            "params.in_str_opt1",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows2(self) -> None:
        wf = DataSourceTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = [
            "ch_in_file1",
            "in_file_opt1",
            "'hello'",
            "'there'",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows3(self) -> None:
        wf = DataSourceTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp3'])
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
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp4'])
        expected = [
            "ch_in_file1",
            "in_file_opt1",
            "params.in_str1",
            "params.in_str_opt1",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows5(self) -> None:
        wf = DataSourceTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp5'])
        expected = [
            "ch_in_file1",
            "in_file_opt1",
            "'hello'",
            "'there'",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_subworkflows6(self) -> None:
        wf = DataSourceTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp6'])
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
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = ["STP1.out.out"]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_workflow_inputs_array(self) -> None:
        wf = ArrayStepInputsTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp1'])
        expected = [
            "ch_in_file_array",
            "in_file_array_opt",
            "params.in_str_array",
            "params.in_int_array",
            "params.in_str_array",
            "params.in_int_array",
            "params.in_str_array",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_static_step_inputs_array(self):
        wf = ArrayStepInputsTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = [
            "ch_in_file_array",
            "null",
            "['hi', 'there', 'friend']",
            "[4, 5, 6]",
            "['hi', 'there', 'friend']",
            "[4, 5, 6]",
            "['hi', 'there', 'friend']",
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)




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
        reset_globals()

    def test_pythontool_array_string_output(self) -> None:
        wf = PlumbingEdgeCaseTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = [
            "STP1.out.out.filter{ it != '' }.map{ it -> it.split(', ') }.ifEmpty( null )"
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)




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





class TestPlumbingExpressions(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_first_selector(self):
        wf = ConditionStepTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['print'])
        expected = [
            "[params.mystring, GET_STRING.out.out].find{ it != null }"
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_with_expression(self):
        wf = StepInputExpressionTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['print'])
        expected = [
            "params.mystring ? params.mystring : params.mystring_backup"
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
        




class TestUnwrapProcess(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()
        wf = UnwrapTestWF()
        do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        self.prescript = process.pre_script
        self.script = process.script
        print(self.prescript)
        print(self.script)

    # PROCESS RELATED
    def test_filename_generated(self) -> None:
        self.assertIn("--filenameGen generated.gz", self.script)
    
    def test_filename_reference(self) -> None:
        self.assertIn("--filenameRef ${in_file.simpleName}.fastq.gz", self.script)

    def test_input_selector_process_input(self) -> None:
        self.assertIn("--inputSelectorProcess ${in_file}", self.script)

    def test_input_selector_param_input(self) -> None:
        self.assertIn("--inputSelectorParam ${params.unwrap_test_tool.in_str}", self.script)
       
    def test_list(self) -> None:
        self.assertIn("--list [1, 2, 3, 4, 5]", self.script) # this is correct
       
    def test_two_value_operator(self) -> None:
        self.assertIn("--TwoValueOperator ${in_file + \".gz\"}", self.script)
    
    def test_first_operator(self) -> None:
        self.assertIn("--FirstOperator ${[params.unwrap_test_tool.in_str, []].find{ it != null }}", self.script)
    
    def test_index_operator(self) -> None:
        self.assertIn("--IndexOperator ${in_file_arr[0]}", self.script)
    
    def test_index_operator_secondaries(self) -> None:
        self.assertIn("--IndexOperatorSecondariesBam ${bam}", self.script)
        self.assertIn("--IndexOperatorSecondariesBai ${bai}", self.script)
    
    def test_index_operator_secondaries_array(self) -> None:
        print(self.script)
        self.assertIn("--IndexOperatorArraySecondariesBams ${indexed_bam_flat[0]}", self.script)
        self.assertIn("--IndexOperatorArraySecondariesBais ${indexed_bam_flat[1]}", self.script)



def update_variables(tool: Tool, vmanager: VariableManager) -> None:
    for tinput in _get_param_inputs(tool, vmanager):
        if utils.is_file_type(tinput.intype):
            if tinput.intype.optional: 
                f_name = nextflow.naming.constructs.gen_varname_file(tinput.id(), dtype=tinput.intype)
                vmanager.update(tinput.id(), 'local', f_name)
            else:
                ch_name = nextflow.naming.constructs.gen_varname_channel(tinput.id(), dtype=tinput.intype)
                vmanager.update(tinput.id(), 'channel', ch_name)

def _get_param_inputs(tool: Tool, vmanager: VariableManager) -> list[TInput]:
    out: list[TInput] = []
    for tinput in tool.tool_inputs():
        var = vmanager.get(tinput.id()).current
        if var.vtype == VariableType.PARAM:
            out.append(tinput)
    return out


class TestUnwrapWorkflow(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_input_node_channel(self) -> None:
        # file input (channel)
        wf = StepConnectionsTestWF()
        do_preprocessing_workflow(wf)
        node = wf.input_nodes['inFile']
        variable_manager = init_variable_manager_for_task(wf)
        update_variables(wf, variable_manager)
        actual = nextflow.unwrap_expression(
            val=node,
            context='workflow',
            variable_manager=variable_manager,
            quote_strings=True
        )
        expected = 'ch_in_file'
        self.assertEqual(actual, expected)

    def test_input_node_param(self) -> None:
        # nonfile input (param)
        wf = StepConnectionsTestWF()
        do_preprocessing_workflow(wf)
        node = wf.input_nodes['inStr']
        variable_manager = init_variable_manager_for_task(wf)
        update_variables(wf, variable_manager)
        actual = nextflow.unwrap_expression(
            val=node,
            context='workflow',
            variable_manager=variable_manager,
            quote_strings=True
        )
        expected = 'params.in_str'
        self.assertEqual(actual, expected)

    def test_step_connection(self) -> None:
        wf = StepConnectionsTestWF()
        do_preprocessing_workflow(wf)
        sources = wf.step_nodes["stp2"].sources
        src = sources['inp']
        variable_manager = init_variable_manager_for_task(wf)
        update_variables(wf, variable_manager)
        actual = nextflow.unwrap_expression(
            val=src,
            context='workflow',
            variable_manager=variable_manager,
            quote_strings=True
        )
        expected = 'STP1.out.out'
        self.assertEqual(actual, expected)
    
    def test_alias_selector_wf(self) -> None:
        wf = AliasSelectorTestWF()
        do_preprocessing_workflow(wf)
        sources = wf.step_nodes["stp2"].sources
        src = sources['inp']
        variable_manager = init_variable_manager_for_task(wf)
        update_variables(wf, variable_manager)
        actual = nextflow.unwrap_expression(
            val=src,
            context='workflow',
            variable_manager=variable_manager,
            quote_strings=True
        )
        expected = 'STP1.out.out'
        self.assertEqual(actual, expected)
    
    def test_first_operator_wf(self) -> None:
        wf = ConditionStepTestWF()
        do_preprocessing_workflow(wf)
        scope = nextflow.Scope()
        sources = wf.step_nodes["print"].sources
        src = sources['inp']
        variable_manager = init_variable_manager_for_task(wf)
        update_variables(wf, variable_manager)
        actual = nextflow.unwrap_expression(
            val=src,
            context='workflow',
            variable_manager=variable_manager,
            quote_strings=True
        )
        expected = '[params.mystring, GET_STRING.out.out].find{ it != null }'
        self.assertEqual(actual, expected)
    
    def test_index_operator_wf(self) -> None:
        wf = IndexOperatorTestWF()
        do_preprocessing_workflow(wf)
        sources = wf.step_nodes["stp1"].sources
        src = sources['inp']
        variable_manager = init_variable_manager_for_task(wf)
        update_variables(wf, variable_manager)
        actual = nextflow.unwrap_expression(
            val=src,
            context='workflow',
            variable_manager=variable_manager,
            quote_strings=True
        )
        expected = 'ch_in_file_arr[0]'
        self.assertEqual(actual, expected)


class TestUnwrapStringFormatter(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()
    
    def test_string_formatter(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = BasicTestTool()
        do_preprocessing_tool(tool)
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter("no format")
        actual = nextflow.unwrap_expression(
            val=sf,
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
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter("there's a {str_arg} arg", str_arg="string")
        actual = nextflow.unwrap_expression(
            val=sf,
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
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter("an input {arg}", arg=InputSelector("testtool"))
        actual = nextflow.unwrap_expression(
            val=sf,
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
        variable_manager = init_variable_manager_for_task(step.tool)
        sf = StringFormatter("an input {arg}", arg=InputSelector("compressionLevel"))
        actual = nextflow.unwrap_expression(
            val=sf, 
            context='process_script',
            variable_manager=variable_manager,
            tool=step.tool,
            in_shell_script=True
        )
        expected = '"an input ${compression_level}"'
        self.assertEqual(actual, expected)

    def test_string_formatter_two_param(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = InputQualityTestTool()
        do_preprocessing_tool(tool)
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter(
            "{username}:{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        actual = nextflow.unwrap_expression(
            val=sf,
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
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter(
            "{username}\\t{password}",
            username=InputSelector("user"),
            password=InputSelector("static"),
        )
        actual_scripting = nextflow.unwrap_expression(
            val=sf,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=False
        )
        actual_shell = nextflow.unwrap_expression(
            val=sf,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual('"user\\tstatic"', actual_scripting)
        self.assertEqual('"${user}\\\\t${static}"', actual_shell)

    def test_expression_arg(self):
        settings.translate.nextflow.MODE = 'tool'
        tool = BasicTestTool()
        do_preprocessing_tool(tool)
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter(
            "{name}:{items}",
            name=InputSelector("testtool"),
            items=JoinOperator(InputSelector("arrayInp"), separator=";"),
        )
        res = nextflow.unwrap_expression(
            val=sf,
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
        variable_manager = init_variable_manager_for_task(step.tool)
        arg = step.tool.arguments()[0]
        actual = nextflow.unwrap_expression(
            val=arg.value, 
            context='process_script',
            variable_manager=variable_manager,
            tool=step.tool,
            in_shell_script=True
        )
        expected = '"-Xmx${8 * 3 / 4}G ${compression_level ? "-Dsamjdk.compress_level=" + compression_level : ""} ${[java_options, []].find{ it != null }.join(" ")}"'
        self.assertEqual(actual, expected)



class TestOrdering(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()
        self.wf = OrderingTestWF()
        do_preprocessing_workflow(self.wf)
        self.processes = nextflow.generate.process.generate_processes(self.wf)
        self.workflows = nextflow.generate.workflow.generate_workflows(self.wf, self.processes)

    def test_process_call(self) -> None:
        # from workflow inputs
        vmanager = init_variable_manager_for_task(self.wf)
        update_variables(self.wf, vmanager)
        
        step = self.wf.step_nodes['stp1']
        task = self.processes[step.tool.id()]
        call = nextflow.generate.workflow.call.gen_task_call(
            alias=step.id(),
            task=task,
            vmanager=vmanager,
            step=step
        )

        actual = simplify_call(call)
        expected = [
            "ch_in_fastq",
            "ch_in_fastq_array",
            "ch_in_file",
            "params.in_int",
            "params.in_int_array",
            "params.in_str",
        ]
        self.assertEqual(expected, actual)
        
        # from process outputs
        step = self.wf.step_nodes['stp3']
        task = self.processes[step.tool.id()]
        call = nextflow.generate.workflow.call.gen_task_call(
            alias=step.id(),
            task=task,
            vmanager=vmanager,
            step=step
        )

        actual = simplify_call(call)
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
        step = self.wf.step_nodes['stp5']
        task = self.processes[step.tool.id()]
        call = nextflow.generate.workflow.call.gen_task_call(
            alias=step.id(),
            task=task,
            vmanager=vmanager,
            step=step
        )

        actual = simplify_call(call)
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
        process = self.processes['MultiTypeTestTool']
        actual_inputs = [i.get_string() for i in process.ordered_inputs]
        expected_inputs = [
            "path in_fastq, stageAs: 'in_fastq.fastq'",
            "path in_fastq_array",
            "path in_file, stageAs: 'in_file'",
            "val in_int",
            "val in_int_array",
            "val in_str",
        ]
        self.assertEqual(actual_inputs, expected_inputs)

    def test_process_directives(self) -> None:
        wf = DirectivesTestWF()
        do_preprocessing_workflow(wf)
        processes = nextflow.generate.process.generate_processes(wf)
        vmanager = init_variable_manager_for_task(wf)
        update_variables(wf, vmanager)
        
        step = wf.step_nodes["stp1"]
        process = processes[step.tool.id()]
        actual_order = [type(x).__name__ for x in process.ordered_directives]
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

    def test_subworkflow_call(self) -> None:
        vmanager = init_variable_manager_for_task(self.wf)
        update_variables(self.wf, vmanager)
        
        # from workflow inputs
        step = self.wf.step_nodes['stp2']
        task = self.workflows[step.tool.id()]
        call = nextflow.generate.workflow.call.gen_task_call(
            alias=step.id(),
            task=task,
            vmanager=vmanager,
            step=step
        )
        actual = simplify_call(call)
        expected = [
            "ch_in_fastq",
            "ch_in_file",
            "ch_in_fastq_array",
            "params.in_int",
            "params.in_int_array",
            "params.in_str",
        ]
        self.assertEqual(expected, actual)
        
        # from process outputs
        step = self.wf.step_nodes['stp4']
        task = self.workflows[step.tool.id()]
        call = nextflow.generate.workflow.call.gen_task_call(
            alias=step.id(),
            task=task,
            vmanager=vmanager,
            step=step
        )

        actual = simplify_call(call)
        expected = [
            "STP1.out.outFastq",
            "STP1.out.outFile",
            "STP1.out.outFastqArray",
            "STP1.out.outInt",
            "STP1.out.outIntArray",
            "STP1.out.outStr",
        ]
        self.assertEqual(expected, actual)

        # from subworkflow outputs
        step = self.wf.step_nodes['stp6']
        task = self.workflows[step.tool.id()]
        call = nextflow.generate.workflow.call.gen_task_call(
            alias=step.id(),
            task=task,
            vmanager=vmanager,
            step=step
        )

        actual = simplify_call(call)
        expected = [
            "STP2.out.outFastq",
            "STP2.out.outFile",
            "STP2.out.outFastqArray",
            "STP2.out.outInt",
            "STP2.out.outIntArray",
            "STP2.out.outStr",
        ]
        self.assertEqual(expected, actual)
    
    def test_subworkflow_inputs(self) -> None:
        workflow = self.workflows['MultiTypeTestWF']
        actual_inputs = [i.get_string() for i in workflow.ordered_take]
        expected_inputs = [
            'ch_in_fastq',
            'ch_in_file',
            'ch_in_fastq_array',
            'ch_in_int',
            'ch_in_int_array',
            'ch_in_str',
        ]
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
        name1 = nextflow.naming.constructs.gen_varname_channel('in_fastq')
        name2 = nextflow.naming.constructs.gen_varname_channel('in_fastq', name_override='in_reads')
        name3 = nextflow.naming.constructs.gen_varname_channel('inFastq')
        self.assertEqual(name1, 'ch_in_fastq')
        self.assertEqual(name2, 'ch_in_reads')
        self.assertEqual(name3, 'ch_in_fastq')
    
