
import unittest
import os
import regex as re
from typing import Any, Optional

from janis_core.translations.common import trace
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
    ComponentsMandatoryTestWF,
    ComponentsOptionalTestWF,
    ComponentsMandatoryArrayTestWF,
    ComponentsOptionalArrayTestWF,
    CallWFInputTestWF,
    CallDuplicateUseageWFInputTestWF,
    CallStaticTestWF,
    CallPartialStaticTestWF,
    CallMinimalTestWF,
    CallConnectionsTestWF,
    CallArrayInputsTestWF,
    CallArrayConnectionsTestWF,
    DirectivesTestWF,

    # arrays
    ArrayIOTestWF,
    ArrayIOExtrasTestWF,

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
    AllFilePairsTestWF,
    FilePairsTestWF0,
    FilePairsTestWF1,
    FilePairsTestWF2,
    FilePairsTestWF3,
    FilePairsOptionalTestWF0,
    FilePairsOptionalTestWF1,
    FilePairsOptionalTestWF2,
    FilePairsOptionalTestWF3,
    FilePairsArrayTestWF,
    FilePairsArrayOptionalTestWF,
    ProcessInputsTestWF,
    OrderingTestWF,
    PlumbingEdgeCaseTestWF,
    MandatoryInputTypesTestWF,
    OptionalInputTypesTestWF,
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
    WorkflowBase,
    CommandTool,
    PythonTool,
    InputSelector,
    StringFormatter,
    JoinOperator,
    TInput,
    Tool,
    WorkflowBuilder,
    Workflow,
    CommandToolBuilder
)

from janis_core.translations import translate
from janis_core.translations import NextflowTranslator as NextflowTranslator
translator = NextflowTranslator()
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

from janis_core.redefinitions.types import (
    Fasta, 
    FastaDict, 
    FastaWithIndexes,
    Fastq, 
    FastqGzPair,
    Bam, 
    BamBai,
    Sam,
    Cram, 
)

from janis_core.translations.nextflow.nfgen_utils import to_groovy
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType
from janis_core import settings

from janis_core.translations.nextflow.generate.workflow.common import get_common_type
from janis_core.translations.nextflow.generate.workflow.datatype_mismatch import (
    get_array_depth,
    is_datatype_mismatch,
    is_array_depth_mismatch,
    gen_datatype_mismatch_plumbing,
)


### helper functions

TRANSLATED_DIR = os.path.join(os.getcwd(), 'translated')

def reset_globals() -> None:
    # reset the translator
    translator = NextflowTranslator()

    # nextflow specific
    settings.translate.MODE = 'extended'
    settings.translate.nextflow.ENTITY = 'workflow'
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
    settings.translate.nextflow.BASE_OUTDIR = os.getcwd()
    nextflow.task_inputs.clear()
    nextflow.params.clear()

def do_preprocessing_workflow(wf: Workflow, ignore_task_inputs: bool=False) -> WorkflowBuilder:
    from janis_core.translations.common import prune_workflow
    from janis_core.translations.common import to_builders
    wf = to_builders(wf)
    if settings.translate.MODE in ['skeleton', 'regular'] and isinstance(wf, WorkflowBuilder):
        assert(isinstance(wf, WorkflowBuilder))
        prune_workflow(wf)
    if not ignore_task_inputs:
        nextflow.preprocessing.populate_task_inputs_workflowmode(wf, wf)
    assert(isinstance(wf, WorkflowBuilder))
    return wf

def do_preprocessing_tool(tool: CommandTool | PythonTool) -> CommandToolBuilder | PythonTool:
    from janis_core.translations.common import to_builders
    tool = to_builders(tool)
    nextflow.preprocessing.populate_task_inputs_toolmode(tool)
    assert(isinstance(tool, CommandToolBuilder | PythonTool))
    return tool

def split_to_lines(call: str) -> list[str]:
    lines = call.split('\n')                        # split
    lines = [ln.strip(' ') for ln in lines]         # strip indents
    lines = [ln for ln in lines if not ln == '']    # remove empty lines
    return lines

def simplify_call(textlines: list[str]) -> list[str]:
    lines = textlines[1:-1]     # get rid of task name & closing brace line
    lines = [ln.split('//')[0] for ln in lines]     # get rid of comments
    lines = [ln.strip(' ') for ln in lines]         # strip indents
    lines = [ln.strip(' ,') for ln in lines]        # strip indents & commas
    lines = [ln for ln in lines if not ln == '']    # remove empty lines
    return lines

def simplify_prescript(text: Optional[str]) -> list[str]:
    if text is None:
        return []
    lines = text.split('\n')                        # split into lines
    lines = [ln.strip(' ') for ln in lines]         # strip indents
    lines = [ln for ln in lines if not ln == '']    # remove empty lines
    return lines

def simplify_script(text: Optional[str]) -> list[str]:
    if text is None:
        return []
    lines = text.split('\n')                        # split into lines
    lines = [ln.strip(' ') for ln in lines]         # strip indents
    lines = [ln.rstrip(' \\') for ln in lines]      # remove escape characters
    lines = [ln for ln in lines if not ln == '']    # remove empty lines
    return lines

def simplify_file(text: str) -> list[str]:
    lines = text.split('\n')                    # split into lines
    lines = [ln.split('//')[0] for ln in lines] # get rid of comments
    lines = [ln.strip() for ln in lines]        # strip indents
    lines = [ln for ln in lines if ln != '']    # remove empty lines
    return lines

def _get_process_directive_lines(text: str) -> list[str]:
    sec_start = r'process.*?{'
    sec_end = r'input:'
    return _lines_within_section(text, sec_start, sec_end)

def _get_process_input_lines(text: str) -> list[str]:
    sec_start = r'input:'
    sec_end = r'(output:)|(script:)|(exec:)'
    return _lines_within_section(text, sec_start, sec_end)

def _get_process_output_lines(text: str) -> list[str]:
    sec_start = r'output:'
    sec_end = r'(script:)|(exec:)'
    return _lines_within_section(text, sec_start, sec_end)

# def _get_process_prescript_lines(text: str) -> list[str]:
#     sec_start = r'script:'
#     sec_end = r'"""'
#     return _lines_within_section(text, sec_start, sec_end)

# def _get_process_script_lines(text: str) -> list[str]:
#     sec_start = r'"""'
#     sec_start = r'"""'
#     return _lines_within_section(text, sec_start, sec_end)

def _get_task_call_lines(text: str, task_name: str) -> list[str]:
    out: list[str] =[]
    lines = simplify_file(text)
    within_section: bool = False
    for ln in lines:
        if ln == f'{task_name}(':
            within_section = True
            continue
        elif ln == ')' and within_section:
            within_section = False
            continue
        elif within_section:
            out.append(ln)
    return out
        
def _lines_within_section(text: str, sec_start: str, sec_end: str) -> list[str]:
    out: list[str] = []
    lines = text.split('\n')                        # split into lines
    lines = [ln.strip(' ') for ln in lines]         # strip indents
    lines = [ln for ln in lines if not ln == '']    # remove empty lines
    lines = [ln for ln in lines if not ln == '"""'] # remove script """ start end
    within_section: bool = False
    for line in lines:
        if re.findall(sec_start, line) and not within_section:
            within_section = True
            continue
        elif re.findall(sec_end, line) and within_section:
            break
        if within_section:
            out.append(line)
    return out

def _gen_call_lines_local(wf: Any, step: Any) -> list[str]:
    reset_globals()
    wf = do_preprocessing_workflow(wf)
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




class TestDatatypeUtils(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()

    def test_get_dtt(self) -> None:
        # secondary array
        testtype = Array(BamBai())
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.SECONDARY_ARRAY)
        
        # secondary 
        testtype = BamBai()
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.SECONDARY)

        # file pair array
        testtype = Array(FastqGzPair())
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.FILE_PAIR_ARRAY)

        # file pair 
        testtype = FastqGzPair()
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.FILE_PAIR)

        # file array 1
        testtype = Array(File())
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.FILE_ARRAY)
        
        # file array 2
        testtype = Array(Fasta())
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.FILE_ARRAY)

        # file 1
        testtype = File()
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.FILE)
        
        # file 2
        testtype = Fasta()
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.FILE)

        # flag array
        testtype = Array(Boolean())
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.FLAG_ARRAY)

        # flag
        testtype = Boolean()
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.FLAG)

        # generic array 1
        testtype = Array(String())
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.GENERIC_ARRAY)
        
        # generic array 2
        testtype = Array(Int())
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.GENERIC_ARRAY)
       
        # generic 1 
        testtype = String()
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.GENERIC)
        
        # generic 2 
        testtype = Int()
        dtt = utils.get_dtt(testtype)
        self.assertEqual(dtt, DTypeType.GENERIC)

        
    def test_is_secondary_array_type(self) -> None:
        # secondary array
        testtype = Array(BamBai())
        verdict = utils.is_secondary_array_type(testtype)
        self.assertTrue(verdict)
        
        # secondary 
        testtype = BamBai()
        verdict = utils.is_secondary_array_type(testtype)
        self.assertFalse(verdict)
        
        # file pair array
        testtype = Array(FastqGzPair())
        verdict = utils.is_secondary_array_type(testtype)
        self.assertFalse(verdict)

        # file array 1
        testtype = Array(File())
        verdict = utils.is_secondary_array_type(testtype)
        self.assertFalse(verdict)
   
    def test_is_secondary_type(self) -> None:
        # secondary array
        testtype = Array(BamBai())
        verdict = utils.is_secondary_type(testtype)
        self.assertTrue(verdict)
        
        # secondary 
        testtype = BamBai()
        verdict = utils.is_secondary_type(testtype)
        self.assertTrue(verdict)
        
        # file pair array
        testtype = Array(FastqGzPair())
        verdict = utils.is_secondary_type(testtype)
        self.assertFalse(verdict)

        # file array 1
        testtype = Array(File())
        verdict = utils.is_secondary_type(testtype)
        self.assertFalse(verdict)

    def test_is_file_pair_array_type(self) -> None:
        # secondary array
        testtype = Array(BamBai())
        verdict = utils.is_file_pair_array_type(testtype)
        self.assertFalse(verdict)
        
        # secondary 
        testtype = BamBai()
        verdict = utils.is_file_pair_array_type(testtype)
        self.assertFalse(verdict)
        
        # file pair array
        testtype = Array(FastqGzPair())
        verdict = utils.is_file_pair_array_type(testtype)
        self.assertTrue(verdict)
    
        # file pair 
        testtype = FastqGzPair()
        verdict = utils.is_file_pair_array_type(testtype)
        self.assertFalse(verdict)

        # file array 1
        testtype = Array(File())
        verdict = utils.is_file_pair_array_type(testtype)
        self.assertFalse(verdict)

    def test_is_file_pair_type(self) -> None:
        # file pair array
        testtype = Array(FastqGzPair())
        verdict = utils.is_file_pair_type(testtype)
        self.assertTrue(verdict)
        
        # file pair 
        testtype = FastqGzPair()
        verdict = utils.is_file_pair_type(testtype)
        self.assertTrue(verdict)

        # file array 1
        testtype = Array(File())
        verdict = utils.is_file_pair_type(testtype)
        self.assertFalse(verdict)

    def test_is_file_array_type(self) -> None:
        # file pair array
        testtype = Array(FastqGzPair())
        verdict = utils.is_file_array_type(testtype)
        self.assertTrue(verdict)
        
        # file pair
        testtype = FastqGzPair()
        verdict = utils.is_file_array_type(testtype)
        self.assertFalse(verdict)

        # file array
        testtype = Array(File())
        verdict = utils.is_file_array_type(testtype)
        self.assertTrue(verdict)

    def test_is_file_type(self) -> None:
        # secondary array
        testtype = Array(BamBai())
        verdict = utils.is_file_type(testtype)
        self.assertTrue(verdict)
        
        # secondary 
        testtype = BamBai()
        verdict = utils.is_file_type(testtype)
        self.assertTrue(verdict)

        # file pair array
        testtype = Array(FastqGzPair())
        verdict = utils.is_file_type(testtype)
        self.assertTrue(verdict)
        
        # file pair
        testtype = FastqGzPair()
        verdict = utils.is_file_type(testtype)
        self.assertTrue(verdict)

        # file array
        testtype = Array(File())
        verdict = utils.is_file_type(testtype)
        self.assertTrue(verdict)
        
        # file 
        testtype = File()
        verdict = utils.is_file_type(testtype)
        self.assertTrue(verdict)

        # generic array
        testtype = Array(String())
        verdict = utils.is_file_type(testtype)
        self.assertFalse(verdict)
       
        # generic 1 
        testtype = Int()
        verdict = utils.is_file_type(testtype)
        self.assertFalse(verdict)

    def test_is_array_type(self) -> None:
        # secondary array
        testtype = Array(BamBai())
        verdict = utils.is_array_type(testtype)
        self.assertTrue(verdict)
        
        # secondary 
        testtype = BamBai()
        verdict = utils.is_array_type(testtype)
        self.assertFalse(verdict)

        # file pair array
        testtype = Array(FastqGzPair())
        verdict = utils.is_array_type(testtype)
        self.assertTrue(verdict)
        
        # file pair
        testtype = FastqGzPair()
        verdict = utils.is_array_type(testtype)
        self.assertFalse(verdict)

        # file array
        testtype = Array(File())
        verdict = utils.is_array_type(testtype)
        self.assertTrue(verdict)
        
        # file 
        testtype = File()
        verdict = utils.is_array_type(testtype)
        self.assertFalse(verdict)

        # generic array
        testtype = Array(String())
        verdict = utils.is_array_type(testtype)
        self.assertTrue(verdict)
       
        # generic 1 
        testtype = Int()
        verdict = utils.is_array_type(testtype)
        self.assertFalse(verdict)

        # TODO
        # utils.is_flag_array_type()
        # utils.is_flag_type()
        # utils.get_base_type()
        # utils.ensure_single_type()
        # utils.get_extensions()
        
        

class TestTaskInputs(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()
        settings.translate.MODE = 'regular'

    # no subworkflows
    def test_one_call(self) -> None:
        wf = MinimalTaskInputsTestWF1()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file',
            'val in_int1',
            'val in_int2',
            'val in_str1',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    def test_two_calls(self) -> None:
        wf = MinimalTaskInputsTestWF2()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file',
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
    
    def test_three_calls(self) -> None:
        wf = MinimalTaskInputsTestWF3()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file',
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
        wf = MinimalTaskInputsTestWF4()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        step = step.tool.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file',
            'val in_str2',
            'val in_str3',
            'val in_int1',
            'val in_int2',
            'val in_int4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)

    def test_two_calls_sub(self) -> None:
        # TODO improve? 
        # main wf process
        wf = MinimalTaskInputsTestWF5()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file',
            'val in_str1',
            'val in_str2',
            'val in_str3',
            'val in_int1',
            'val in_int2',
            'val in_int4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)

        # sub wf process - just tests that the inputs in the 2nd usage 
        # are the exact same as the first usage
        step = wf.step_nodes["stp2"]
        step = step.tool.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file',
            'val in_str1',
            'val in_str2',
            'val in_str3',
            'val in_int1',
            'val in_int2',
            'val in_int4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
    
    def test_three_calls_sub(self) -> None:
        # TODO improve? 
        wf = MinimalTaskInputsTestWF6()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        expected_inputs = {
            'path in_file',
            'val in_int1',
            'val in_int2',
            'val in_int3',
            'val in_int4',
            'val in_str2',
            'val in_str3',
            'val in_str4',
        }
        actual_inputs = {inp.get_string() for inp in process.inputs}
        self.assertEqual(actual_inputs, expected_inputs)
        


# class TestPreprocessingTaskInputs(unittest.TestCase):

#     def setUp(self) -> None:
#         reset_globals()
#         settings.translate.MODE = 'regular'

#     def test_main_wf(self) -> None:
#         wf = MinimalTaskInputsTestWF1()
#         wf = do_preprocessing_workflow(wf)

#         # main wf
#         actual_task_inputs = nextflow.task_inputs.task_inputs(wf.id())
#         expected_task_inputs: set[str] = set()
#         self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
#         actual_param_inputs = nextflow.task_inputs.param_inputs(wf.id())
#         expected_param_inputs = {'inFile', 'inStr1', 'inInt1'}
#         self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
#         actual_static_inputs = nextflow.task_inputs.static_inputs(wf.id())
#         expected_static_inputs: set[str] = set()
#         self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
#         actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(wf.id())
#         expected_ignored_inputs = set()
#         self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)

#     # no subworkflows
#     def test_one_call(self) -> None:
#         wf = MinimalTaskInputsTestWF1()
#         wf = do_preprocessing_workflow(wf)
#         tool = wf.step_nodes['stp1'].tool

#         actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
#         expected_task_inputs = {'inFile', 'inStr1', 'inInt2'}
#         self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
#         actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
#         expected_param_inputs = set()
#         self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
#         actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
#         expected_static_inputs = {'inInt1'}
#         self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
#         actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
#         expected_ignored_inputs = set()
#         self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)
    
#     def test_two_calls(self) -> None:
#         wf = MinimalTaskInputsTestWF2()
#         wf = do_preprocessing_workflow(wf)
#         tool = wf.step_nodes['stp1'].tool

#         actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
#         expected_task_inputs = {'inFile', 'inStr1', 'inStr2', 'inStr3', 'inStr4', 'inInt1', 'inInt2', 'inInt3'}
#         self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
#         actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
#         expected_param_inputs = set()
#         self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
#         actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
#         expected_static_inputs = set()
#         self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
#         actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
#         expected_ignored_inputs = set()
#         self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)
    
#     def test_three_calls(self) -> None:
#         wf = MinimalTaskInputsTestWF3()
#         wf = do_preprocessing_workflow(wf)
#         tool = wf.step_nodes['stp1'].tool

#         actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
#         expected_task_inputs = {
#             'inFile', 
#             'inStr1', 
#             'inStr2', 
#             'inStr3', 
#             'inStr4', 
#             'inInt1', 
#             'inInt2', 
#             'inInt3',
#             'inInt4',
#         }
#         self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
#         actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
#         expected_param_inputs = set()
#         self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
#         actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
#         expected_static_inputs = set()
#         self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
#         actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
#         expected_ignored_inputs = set()
#         self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)


#     # subworkflows
#     def test_one_call_sub(self) -> None:
#         wf = MinimalTaskInputsTestWF4()
#         wf = do_preprocessing_workflow(wf)
#         subwf = wf.step_nodes['stp1'].tool
#         tool = subwf.step_nodes['stp1'].tool

#         actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
#         actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
#         actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
#         actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        
#         expected_task_inputs = {'inFile', 'inStr1', 'inInt2'}
#         # self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
#         expected_param_inputs = set()
#         self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
#         expected_static_inputs = {'inInt1'}
#         self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
#         expected_ignored_inputs = set()
#         self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)
        
#         # actual_task_inputs = nextflow.task_inputs.task_inputs(subwf.id())
#         # expected_task_inputs = {'inFile', 'inStr1', 'inInt2'}
#         # self.assertSetEqual(actual_task_inputs, expected_task_inputs)
        
#         # actual_param_inputs = nextflow.task_inputs.param_inputs(subwf.id())
#         # expected_param_inputs = set()
#         # self.assertSetEqual(actual_param_inputs, expected_param_inputs)
        
#         # actual_static_inputs = nextflow.task_inputs.static_inputs(subwf.id())
#         # expected_static_inputs = {'inInt1'}
#         # self.assertSetEqual(actual_static_inputs, expected_static_inputs)
        
#         # actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(subwf.id())
#         # expected_ignored_inputs = set()
#         # self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)

#     def test_two_calls_sub(self) -> None:
#         wf = MinimalTaskInputsTestWF5()
#         wf = do_preprocessing_workflow(wf)
        
#         # TaskInputsTestTool1
#         tool = wf.step_nodes['stp1'].tool
#         actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
#         actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
#         actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
#         actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        
#         expected_task_inputs = {
#             'inFile',
#             'inStr1',
#             'inStr2',
#             'inStr3',
#             'inInt2',
#             'inInt4',
#         }
#         expected_param_inputs = set()
#         expected_static_inputs = {'inInt1'}
#         expected_ignored_inputs = {
#             'inStr4',
#             'inInt3',
#         }
        
#         self.assertSetEqual(actual_task_inputs, expected_task_inputs)
#         self.assertSetEqual(actual_param_inputs, expected_param_inputs)
#         self.assertSetEqual(actual_static_inputs, expected_static_inputs)
#         self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)
        
#         # SubMinimalTaskInputsTestWF
#         tool = wf.step_nodes['stp2'].tool

#         actual_task_inputs = nextflow.task_inputs.task_inputs(tool.id())
#         actual_param_inputs = nextflow.task_inputs.param_inputs(tool.id())
#         actual_static_inputs = nextflow.task_inputs.static_inputs(tool.id())
#         actual_ignored_inputs = nextflow.task_inputs.ignored_inputs(tool.id())
        
#         expected_task_inputs = {'inFile'}
#         expected_param_inputs = {'inStr1', 'inInt2'}
#         expected_static_inputs = set()
#         expected_ignored_inputs = {
#             'inStr2',
#             'inStr3',
#             'inInt1',
#             'inInt3',
#         }
        
#         self.assertSetEqual(actual_task_inputs, expected_task_inputs)
#         self.assertSetEqual(actual_param_inputs, expected_param_inputs)
#         self.assertSetEqual(actual_static_inputs, expected_static_inputs)
#         self.assertSetEqual(actual_ignored_inputs, expected_ignored_inputs)





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
        expected = '"Hello"'
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
        # self.assertEquals(settings.translate.nextflow.LIB_FILENAME, 'lib.nf')
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
        mainstr, _, subtask_dict = translate(wf, dest_fmt='nextflow')
        
        # check correct number of subtasks created
        self.assertEqual(len(subtask_dict), 6)
        expected_filepaths = set([
            'modules/file_test_tool.nf',
            'modules/string_test_tool.nf',
            'modules/int_test_tool.nf',
            'modules/string_opt_test_tool.nf',
            'subworkflows/oranges_workflow.nf',
            'subworkflows/apples_workflow.nf',
        ])
        actual_filepaths = set([x[0] for x in subtask_dict])
        self.assertSetEqual(actual_filepaths, expected_filepaths)

    def test_main_workflow_format(self) -> None:
        wf = AssemblyTestWF()
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)
        expected = [
            "nextflow.enable.dsl=2",
            "include { FASTQC as FASTQC1 } from './modules/fastqc'",
            "include { FASTQC as FASTQC2 } from './modules/fastqc'",
            "include { FASTQC as FASTQC3 } from './modules/fastqc'",
            "include { CAT_TEST_TOOL } from './modules/cat_test_tool'",
            "include { UNICYCLER } from './modules/unicycler'",
            "fastqc1_adapters      = file( params.fastqc1_adapters )",
            "fastqc1_contaminants  = file( params.fastqc1_contaminants )",
            "fastqc1_limits        = file( params.fastqc1_limits )",
            "fastqc2_adapters      = file( params.fastqc2_adapters )",
            "fastqc2_contaminants  = file( params.fastqc2_contaminants )",
            "fastqc2_limits        = file( params.fastqc2_limits )",
            "in_forward_reads      = file( params.in_forward_reads )",
            "in_long_reads         = file( params.in_long_reads )",
            "in_reverse_reads      = file( params.in_reverse_reads )",
            "test_input            = file( params.test_input )",
            "workflow {",
        ]
        actual = simplify_file(mainstr)
        actual = actual[:len(expected)]
        for ln in actual:
            self.assertIn(ln, expected)

    @unittest.skip("Nextflow translation rewrites expected output, need a new test")
    def test_process_format(self) -> None:
        wf = AssemblyTestWF()
        mainstr, _, subtasks = translate(wf, dest_fmt='nextflow', to_console=False)
        process = [x[1] for x in subtasks if x[0] == 'modules/fastqc.nf'][0]
        print(process)
        actual_lines = simplify_file(process)
        expected_lines = [
            'nextflow.enable.dsl=2',
            'process FASTQC {',
            'container "quay.io/biocontainers/fastqc:0.11.8--2"',
            'publishDir "${params.outdir}/fastqc"',
            'input:',
            'path input_file',
            'path adapters, stageAs: \'adapters/*\'',
            'path contaminants, stageAs: \'contaminants/*\'',
            'path limits, stageAs: \'limits/*\'',
            'val extract',
            'val nogroup',
            'val quiet',
            'val kmers',
            'val min_length',
            'val option_f',
            'val outdir',
            'output:',
            'path "output.html", emit: outHtmlFile',
            'path "output.txt", emit: outTextFile',
            'script:',
            'def adapters = adapters.simpleName != params.NULL_VALUE ? "--adapters ${adapters}" : ""',
            'def contaminants = contaminants.simpleName != params.NULL_VALUE ? "--contaminants ${contaminants}" : ""',
            'def limits = limits.simpleName != params.NULL_VALUE ? "--limits ${limits}" : ""',
            'def extract = extract ? "--extract" : ""',
            'def nogroup = nogroup ? "--nogroup" : ""',
            'def quiet = quiet ? "--quiet" : ""',
            'def kmers = kmers != params.NULL_VALUE ? kmers : 7',
            'def min_length = min_length != params.NULL_VALUE ? "--min_length ${min_length}" : ""',
            'def option_f = option_f != params.NULL_VALUE ? "-f ${option_f}" : ""',
            'def outdir = outdir != params.NULL_VALUE ? "--outdir ${outdir}" : ""',
            '"""',
            'fastqc \\',
            '${adapters} \\',
            '${contaminants} \\',
            '${limits} \\',
            '--kmers ${kmers} \\',
            '${min_length} \\',
            '${outdir} \\',
            '${option_f} \\',
            '${extract} \\',
            '${nogroup} \\',
            '${quiet} \\',
            '${input_file}',
            '"""',
            '}',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_subworkflow_format(self) -> None:
        wf = SubworkflowTestWF()
        mainstr, _, subtasks = translate(wf, dest_fmt='nextflow', to_console=False)
        process = [x[1] for x in subtasks if x[0] == 'subworkflows/apples_workflow.nf'][0]
        print(process)
        actual_lines = simplify_file(process)
        expected_lines = [
            "nextflow.enable.dsl=2",
            "include { STRING_TEST_TOOL as STRING_TOOL } from '../modules/string_test_tool'",
            "include { STRING_OPT_TEST_TOOL as STRING_OPT_TOOL } from '../modules/string_opt_test_tool'",
            "include { ORANGES_WORKFLOW as ORANGES_SUBWORKFLOW } from './oranges_workflow'",
            "workflow APPLES_WORKFLOW {",
            "take:",
            "ch_in_int",
            "ch_in_str",
            "ch_in_str_opt",
            "main:",
            "STRING_TOOL(",
            "ch_in_str",
            ")",
            "STRING_OPT_TOOL(",
            "ch_in_str_opt",
            ")",
            "ORANGES_SUBWORKFLOW(",
            "STRING_TOOL.out.out,",
            "ch_in_int",
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
        _, config, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        actual_lines = split_to_lines(config)
        expected_lines = [
            'nextflow.enable.dsl = 2',
            'singularity.enabled = true',
            'singularity.autoMounts = true',
            'singularity.cacheDir = "$HOME/.singularity/cache"',
            'params {',
            '// Placeholder for null values.',
            '// Do not alter unless you know what you are doing.',
            'NULL_VALUE = \'NULL\'',
            '// WORKFLOW OUTPUT DIRECTORY',
            'outdir  = \'./outputs\'',
            '// INPUTS (MANDATORY)',
            'in_file  = NULL_VALUE  // (MANDATORY generic file)',
            'in_int   = NULL_VALUE  // (MANDATORY integer)',
            'in_str   = NULL_VALUE  // (MANDATORY string)',
            '// INPUTS (OPTIONAL)',
            'in_str_opt  = NULL_VALUE  // (optional string)',
            '}',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_config_params(self) -> None:
        # test_nonfile
        # string, int, bool
        wf = AllInputTypesTestWF()
        _, config, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        actual_lines = split_to_lines(config)
        for ln in actual_lines:
            print(ln)
        expected_lines = [
            'nextflow.enable.dsl = 2',
            'singularity.enabled = true',
            'singularity.autoMounts = true',
            'singularity.cacheDir = "$HOME/.singularity/cache"',
            'params {',
            '// Placeholder for null values.',
            '// Do not alter unless you know what you are doing.',
            'NULL_VALUE = \'NULL\'',
            '// WORKFLOW OUTPUT DIRECTORY',
            'outdir  = \'./outputs\'',
            '// INPUTS (MANDATORY)',
            'in_file               = NULL_VALUE  // (MANDATORY generic file)',
            'in_file_array         = []          // (MANDATORY array)         eg. [file1, ...]',
            'in_secondaries        = []          // (MANDATORY indexedbam)    eg. [bam, bai]',
            'in_secondaries_array  = [[]]        // (MANDATORY array)         eg. [[bam, bai]]',
            'in_filepair           = []          // (MANDATORY fastqpair)     eg. [pair1, pair2]',
            'in_filepair_array     = [[]]        // (MANDATORY array)         eg. [[pair1, pair2]]',
            'in_nonfile            = NULL_VALUE  // (MANDATORY integer)',
            'in_nonfile_array      = NULL_VALUE  // (MANDATORY array)         eg. [integer1, ...]',
            '// INPUTS (OPTIONAL)',
            'in_file_array_optional         = []          // (optional array)         eg. [file1, ...]',
            'in_file_optional               = NULL_VALUE  // (optional generic file)',
            'in_secondaries_array_optional  = [[]]        // (optional array)         eg. [[bam, bai]]',
            'in_secondaries_optional        = []          // (optional indexedbam)    eg. [bam, bai]',
            'in_filepair_array_optional     = [[]]        // (optional array)         eg. [[pair1, pair2]]',
            'in_filepair_optional           = []          // (optional fastqpair)     eg. [pair1, pair2]',
            'in_nonfile_array_optional      = NULL_VALUE  // (optional array)         eg. [integer1, ...]',
            'in_nonfile_optional            = NULL_VALUE  // (optional integer)',
            '}',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)
    
    def test_config_params_pythontool(self) -> None:
        # test_nonfile
        # string, int, bool
        
        wf = InputsPythonToolTestWF()
        _, config, _ = translate(wf, dest_fmt='nextflow', to_console=False, export_path='translated')
        actual_lines = split_to_lines(config)
        for ln in actual_lines:
            print(ln)
        expected_lines = [
            'nextflow.enable.dsl = 2',
            'singularity.enabled = true',
            'singularity.autoMounts = true',
            'singularity.cacheDir = "$HOME/.singularity/cache"',
            'params {',
            '// Placeholder for null values.',
            '// Do not alter unless you know what you are doing.',
            'NULL_VALUE = \'NULL\'',
            '// WORKFLOW OUTPUT DIRECTORY',
            'outdir  = \'./outputs\'',
            '// INPUTS (MANDATORY)',
            'in_file            = NULL_VALUE  // (MANDATORY generic file)',
            'in_secondary_type  = []          // (MANDATORY generic file)  eg. [txt, txt]',
            'in_int             = NULL_VALUE  // (MANDATORY integer)',
            'in_str             = NULL_VALUE  // (MANDATORY string)',
            'in_str_arr         = NULL_VALUE  // (MANDATORY array)         eg. [string1, ...]',
            '// PROCESS: JOIN_ARRAY_PYTHON_TEST_TOOL',
            f'join_array_python_test_tool.code_file  = "{TRANSLATED_DIR}/templates/JoinArrayPythonTestTool.py"',
            '// PROCESS: MULTI_TYPES_INPUT_PYTHON_TOOL',
            f'multi_types_input_python_tool.code_file  = "{TRANSLATED_DIR}/templates/MultiTypesInputPythonTool.py"',
            '// PROCESS: SECONDARY_INPUT_PYTHON_TEST_TOOL',
            f'secondary_input_python_test_tool.code_file  = "{TRANSLATED_DIR}/templates/SecondaryInputPythonTestTool.py"',
            '}',
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
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        actual_lines = simplify_file(mainstr)
        actual_lines = [ln for ln in actual_lines if 'Channel.' in ln]
        expected_lines = [
            'ch_in_file_array         = Channel.fromPath( params.in_file_array ).toList()',
            'ch_in_secondaries_array  = Channel.fromPath( params.in_secondaries_array.flatten() ).collate( 2 )',
            'ch_in_filepair_array     = Channel.of( params.in_filepair_array ).toList()',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

    def test_variable_declarations(self) -> None:
        wf = AllInputTypesTestWF()
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        actual_lines = simplify_file(mainstr)
        actual_lines = [ln for ln in actual_lines if 'file(' in ln]
        expected_lines = [
            'in_file                        = file( params.in_file )',
            'in_file_array_optional         = params.in_file_array_optional.collect{ file(it) }',
            'in_file_optional               = file( params.in_file_optional )',
            'in_filepair                    = params.in_filepair.collect{ file(it) }',
            'in_filepair_array_optional     = params.in_filepair_array_optional.collect{ it.collect{ file(it) } }',
            'in_filepair_optional           = params.in_filepair_optional.collect{ file(it) }',
            'in_secondaries                 = params.in_secondaries.collect{ file(it) }',
            'in_secondaries_array_optional  = params.in_secondaries_array_optional.collect{ it.collect{ file(it) } }',
            'in_secondaries_optional        = params.in_secondaries_optional.collect{ file(it) }',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in actual_lines:
            self.assertIn(ln, expected_lines)

    def test_duplicate_tool_usage(self) -> None:
        wf = DuplicateTasksTestWF()
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)

        # main workflow fmt
        actual_lines = simplify_file(mainstr)
        expected_lines = [
            "include { ECHO_TEST_TOOL as STP1 } from './modules/echo_test_tool'",
            "include { ECHO_TEST_TOOL as STP2 } from './modules/echo_test_tool'",
            "include { ECHO_TEST_WORKFLOW1 as STP3 } from './subworkflows/echo_test_workflow1'",
            "include { ECHO_TEST_WORKFLOW1 as STP4 } from './subworkflows/echo_test_workflow1'",
            "include { ECHO_TEST_WORKFLOW2 as STP5 } from './subworkflows/echo_test_workflow2'",
            "STP1(",
            "STP2(",
            "STP3(",
            "STP4(",
            "STP5(",
        ]
        for ln in expected_lines:
            self.assertIn(ln, actual_lines)

    @unittest.skip('TODO')
    def test_duplicate_subworkflow_usage(self) -> None:
        raise NotImplementedError
    



class TestCmdtoolProcessDirectives(unittest.TestCase):
    """
    Tests identifying tool inputs which should be process inputs.
    Need a process input for each tool input in step sources

    INCLUDES WORKFLOW OUTPUTS as publishDir
    """

    def setUp(self) -> None:
        self.destfmt = 'nextflow'
        reset_globals()

    def test_directives(self) -> None:
        wf = DirectivesTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.destfmt, export_path='./translated')
        subtasks.sort(key=lambda x: x[0])
        actual_lines = _get_process_directive_lines(subtasks[1][1])
        expected_directives = {
            'container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"',
            'publishDir "${params.outdir}/resources_test_tool"',
            'cpus "${params.resources_test_tool.cpus}"',
            'disk "${params.resources_test_tool.disk}"',
            'memory "${params.resources_test_tool.memory}"',
            'time "${params.resources_test_tool.time}"'
        }
        for direc in expected_directives:
            self.assertIn(direc, actual_lines)
    
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
        self.dest = 'nextflow'
        reset_globals()
    
    def test_secondaries(self) -> None:
        wf = SecondariesTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/secondaries_test_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            'path bam1',
            'path bam2',
            'path bam3',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
    
    def test_secondaries_optional(self) -> None:
        wf = SecondariesTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/secondaries_optional_test_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            "path bam1, stageAs: 'bam1/*'"
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
    
    def test_secondaries_array(self) -> None:
        wf = SecondariesTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/secondaries_array_test_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            'path bams1_flat',
            'path bams2_flat',
            'path bams3_flat',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)

    def test_secondaries_array_optional(self) -> None:
        wf = SecondariesTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/secondaries_array_optional_test_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            "path bams1_flat",
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)

    def test_file_pair(self) -> None:
        wf = FilePairsTestWF0()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/file_pair_test_tool0.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            'tuple path(reads1), path(reads2)',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)

    def test_file_pair_optional(self) -> None:
        wf = FilePairsOptionalTestWF0()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/file_pair_optional_test_tool0.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            'tuple path(reads1), path(reads2)',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
    
    def test_file_pair_array(self) -> None:
        wf = FilePairsArrayTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/file_pair_array_test_tool1.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            'path read_pairs_flat',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
    
    def test_file_pair_array_optional(self) -> None:
        wf = FilePairsArrayOptionalTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/file_pair_array_optional_test_tool1.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            'path read_pairs_flat',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
        
    # a rather weak test
    def test_generics(self) -> None:
        wf = ComponentsMandatoryTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/components_mandatory_test_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            'path pos_basic',
            'val flag_false',
            'val flag_true',
            'val opt_basic',
            'val opt_default',
            'val pos_default',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
    
    # @unittest.skip('filenames are scuffed')
    def test_filenames(self) -> None:
        wf = FilenameTestWF1()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        
        # basic
        process = [x[1] for x in subtasks if x[0] == 'modules/filename_test_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            "path inp1",
            "val inp2",
        ]
        print(process)
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
        
        # filename using input selector
        process = [x[1] for x in subtasks if x[0] == 'modules/filename_input_selector_test_tool.nf'][0]
        print(process)
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            "path inp1",
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)

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
        self.dest = 'nextflow'

    def test_stdout(self):
        wf = BasicIOTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/file_test_tool.nf'][0]
        actual_lines = _get_process_output_lines(process)
        expected_lines = [
            'path "out", emit: out',
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)

    def test_wildcard(self) -> None:
        wf = OutputCollectionTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/wildcard_selector_test_tool.nf'][0]
        actual_lines = _get_process_output_lines(process)
        expected_lines = [
            'path "myfile.txt", emit: out'
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
    
    def test_wildcard_array(self) -> None:
        wf = WildcardSelectorOutputTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/array_wildcard_selector_test_tool.nf'][0]
        actual_lines = _get_process_output_lines(process)
        expected_lines = [
            'path "*.txt", emit: out'
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
        
    def test_input_selector(self) -> None:
        wf = InputSelectorTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/file_input_selector_test_tool.nf'][0]
        actual_lines = _get_process_output_lines(process)
        expected_lines = [
            'path inp, emit: out'
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)

    def test_input_selector_param(self) -> None:
        wf = OutputCollectionTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/input_selector_test_tool.nf'][0]
        actual_lines = _get_process_output_lines(process)
        expected_lines = [
            'path output_filename, emit: out'
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
        
    def test_input_selector_array(self) -> None:
        wf = InputSelectorTestWF()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/array_input_selector_test_tool.nf'][0]
        actual_lines = _get_process_output_lines(process)
        expected_lines = [
            'path inp, emit: out'
        ]
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)

    def test_filenames_generated(self) -> None:
        wf = FilenameTestWF2()
        _, _, subtasks = translate(wf, dest_fmt=self.dest, export_path='./translated')
        process = [x[1] for x in subtasks if x[0] == 'modules/filename_collection_test_tool.nf'][0]
        actual_lines = _get_process_output_lines(process)
        expected_lines = [
            'path "${inp1.baseName + ".csv"}", emit: out3',
            'path "${"generated" + ".csv"}", emit: out4',
            'path "generated.csv", emit: out5',
            'path "generated.merged.csv", emit: out6'
        ]
        print(process)
        self.assertEqual(len(actual_lines), len(expected_lines))
        for inp in expected_lines:
            self.assertIn(inp, actual_lines)
        
    def test_filenames_referenced(self) -> None:
        # inp1, inp1_1, inp4, inp5, inp6 are in step.sources
        wf = FilenameTestWF1()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp5']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {
            'path "${inp1.baseName + ".csv"}", emit: out3',
            'path "${inp4 + ".csv"}", emit: out4',
            'path "${inp5 + ".csv"}", emit: out5',
            'path "${inp6 + ".merged.csv"}", emit: out6'
        }
        print(process.get_string())
        self.assertEqual(actual_outputs, expected_outputs)

    def test_file_pair(self) -> None:
        # eg read1.fastq, read2.fastq
        # collection method is list, len(list) == 2.
        wf = OutputCollectionTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp6']
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "[${inp.baseName + "-R1.fastq"}, ${inp.baseName + "-R2.fastq"}]", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_secondaries(self) -> None:
        wf = SecondariesTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'tuple path("*.bam"), path("*.bam.bai"), emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_secondaries_replaced(self) -> None:
        wf = SecondariesTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp5']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'tuple path("*.bam"), path("*.bai"), emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)

    def test_secondaries_edge_basename(self) -> None:
        wf = SecondariesTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp6']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        print(process.outputs[0].get_string())
        expected_outputs = {
            'tuple path("${bam[0].name}"), path("*.bai"), emit: out'
        }
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_secondaries_edge_no_secondaries_present_as(self) -> None:
        wf = SecondariesTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp7']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        print(process.outputs[0].get_string())
        expected_outputs = {
            'tuple path(inp), path("*.tbi"), emit: out'
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
        wf = OutputCollectionTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp5']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        expected_outputs = {'path "${inp.baseName + ".gz"}", emit: out'}
        self.assertEqual(actual_outputs, expected_outputs)
    
    def test_edge_markduplicates_metrics(self) -> None:
        wf = OutputCollectionTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp8']
        process = nextflow.generate.process.generate_process(step.tool)
        actual_outputs = {out.get_string() for out in process.outputs}
        print(process.get_string())
        expected_outputs = {
            'path "${[output_prefix, "generated"].find{ it != null } + ".metrics.txt"}", emit: metrics'
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

    @unittest.skip('not implemented')
    def test_get_ordered_inputs_arguments(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_prescript_presence(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_script_presence(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_autofill(self) -> None:
        raise NotImplementedError
    
    @unittest.skip('not implemented')
    def test_multiple_statements(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_directories_to_create(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_files_to_create_cmdtool(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_files_to_create_codetool(self) -> None:
        raise NotImplementedError

    @unittest.skip('not implemented')
    def test_files_to_create_cmdtool_exprtool(self) -> None:
        raise NotImplementedError

    def test_variables_defined1(self) -> None:
        # inputs referencing undefined inputs
        wf = EntityTraceTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp8"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        assert(actual_prescript)
        expected_lines = {
            'def java_options_joined = java_options != params.NULL_VALUE ? java_options.join(\' \') : ""',
            'def compression_level = compression_level != params.NULL_VALUE ? compression_level : ""'
        }

        for ln in expected_lines:
            self.assertIn(ln, actual_prescript)
        
    def test_components_mandatory(self) -> None:
        wf = ComponentsMandatoryTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())

        # pre-script
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def pos_default = pos_default != params.NULL_VALUE ? pos_default : 95',
            'def flag_true = flag_true == false ? "" : "--flag-true"',
            'def flag_false = flag_false ? "--flag-false" : ""',
            'def opt_default = opt_default != params.NULL_VALUE ? opt_default : 5',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in actual_prescript:
            self.assertIn(ln, expected_prescript)
        
        # script
        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${pos_basic}',
            '${pos_default}',
            '${flag_true}',
            '${flag_false}',
            '--opt-basic ${opt_basic}',
            '--opt-default ${opt_default}',
            '> out'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)
    
    def test_components_optional(self) -> None:
        wf = ComponentsOptionalTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        
        # pre-script
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def pos_optional = pos_optional != params.NULL_VALUE ? pos_optional : ""',
            'def flag_true = flag_true == false ? "" : "--flag-true"',
            'def flag_false = flag_false ? "--flag-false" : ""',
            'def opt_optional = opt_optional != params.NULL_VALUE ? "--opt-optional ${opt_optional}" : ""',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in actual_prescript:
            self.assertIn(ln, expected_prescript)
        
        # script
        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${pos_optional}',
            '${flag_true}',
            '${flag_false}',
            '${opt_optional}',
            '> out'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)
    
    def test_components_array_mandatory(self) -> None:
        wf = ComponentsMandatoryArrayTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        
        # pre-script
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def pos_basic_arr_joined = pos_basic_arr.join(\' \')',
            'def opt_basic_arr_joined = opt_basic_arr.join(\' \')',
            'def opt_basic_arr_prefixeach_joined = opt_basic_arr_prefixeach.collect{ "--opt-basic-prefixeach ${it}" }.join(\' \')',
            'def opt_default_arr_joined = opt_default_arr != params.NULL_VALUE ? opt_default_arr.join(\' \') : "100 200 300"',
            'def opt_default_arr_prefixeach_joined = opt_default_arr_prefixeach != params.NULL_VALUE ? opt_default_arr_prefixeach.collect{ "--opt-default-prefixeach ${it}" }.join(\' \') : "--opt-default-prefixeach "hi" --opt-default-prefixeach "there""',
            'def pos_default_arr_joined = pos_default_arr != params.NULL_VALUE ? pos_default_arr.join(\' \') : "1 2"',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in actual_prescript:
            self.assertIn(ln, expected_prescript)
        
        # script
        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${pos_basic_arr_joined}',
            '${pos_default_arr_joined}',
            '--opt-basic ${opt_basic_arr_joined}',
            '--opt-default ${opt_default_arr_joined}',
            '${opt_basic_arr_prefixeach_joined}',
            '${opt_default_arr_prefixeach_joined}',
            '> out'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)
    
    def test_components_array_optional(self) -> None:
        wf = ComponentsOptionalArrayTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes['stp1']
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        
        # pre-script
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def pos_optional_arr_joined = pos_optional_arr[0] != null ? pos_optional_arr.join(\' \') : ""',
            'def opt_optional_arr_joined = opt_optional_arr != params.NULL_VALUE ? "--opt-optional-arr " + opt_optional_arr.join(\' \') : ""',
            'def opt_optional_arr_prefixeach_joined = opt_optional_arr_prefixeach != params.NULL_VALUE ? opt_optional_arr_prefixeach.collect{ "--opt-optional-arr-prefixeach ${it}" }.join(\' \') : ""',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in actual_prescript:
            self.assertIn(ln, expected_prescript)
        
        # script
        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${pos_optional_arr_joined}',
            '${opt_optional_arr_joined}',
            '${opt_optional_arr_prefixeach_joined}',
            '> out'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)
        
    def test_secondaries(self) -> None:
        wf = SecondariesTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        # pre-script
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def bam1 = bam1[0]',
            'def bam2 = bam2[0]',
            'def bam3 = bam3[0]',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in actual_prescript:
            self.assertIn(ln, expected_prescript)
        
        # script
        actual_script = simplify_script(process.script)
        expected_lines = [
            'echo',
            '${bam1}',
            '${bam2}',
            '${bam3}',
            '--arg1 ${bam1}',
            '--arg2 ${bam1}',
            '--arg3 ${bam1}',
        ]
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)
    
    def test_secondaries_optional(self) -> None:
        wf = SecondariesTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp2"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        # pre-script
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def bam1 = bam1[0] != null ? bam1[0] : ""',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in actual_prescript:
            self.assertIn(ln, expected_prescript)
        
        # script
        actual_script = simplify_script(process.script)
        expected_lines = [
            'echo',
            '${bam1}',
        ]
        print(process.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

    def test_secondaries_array(self) -> None:
        wf = SecondariesTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp3"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        
        # pre-script
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = {
            'def bams1 = get_primary_files(bams1_flat, 2)',
            'def bams1_joined = bams1.join(\' \')',
            'def bams2 = get_primary_files(bams2_flat, 2)',
            'def bams2_joined = bams2.join(\' \')',
            'def bams3 = get_primary_files(bams3_flat, 2)',
            'def bams3_joined = bams3.collect{ "--bams3 ${it}" }.join(\' \')',
        }
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)
        
        # script
        actual_script = simplify_script(process.script)
        expected_script = {
            'echo',
            '${bams1_joined}',
            '--bams2 ${bams2_joined}',
            '${bams3_joined}',
            '--bams-arg1 ${bams1_joined}',
            '--bams-arg2 ${bams1[0]}',
            '--bams-arg3 ${bams1[1]}',
            '> out'
        }
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)
    
    def test_secondaries_array_optional(self) -> None:
        wf = SecondariesTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp4"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())
        
        # pre-script
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = {
            'def bams1 = get_primary_files(bams1_flat, 2)',
            'def bams1_joined = bams1[0] != null ? bams1.join(\' \') : ""',
        }
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)
        
        # script
        actual_script = simplify_script(process.script)
        expected_script = {
            'echo',
            '${bams1_joined}',
            '> out'
        }
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_file_pair0(self) -> None:
        wf = FilePairsTestWF0()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def reads_joined = reads1 + \' \' + reads2'
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${reads_joined}',
            '> stdout'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)
    
    def test_file_pair1(self) -> None:
        wf = FilePairsTestWF1()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def reads_a_joined = "--prefix ${reads_a1} ${reads_a2}"',
            'def reads_b_joined = "--prefixeach ${reads_b1} --prefixeach ${reads_b2}"',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${reads_a_joined}',
            '${reads_b_joined}',
            '> stdout'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_file_pair2(self) -> None:
        wf = FilePairsTestWF2()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def reads_joined = reads1 + \' \' + reads2',
            'def read1 = read1.simpleName != params.NULL_VALUE ? read1 : ${reads1}',
            'def read2 = read2.simpleName != params.NULL_VALUE ? read2 : ${reads2}',

        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${reads_joined}',
            '--reads-index-0 ${read1}',
            '--reads-index-1 ${read2}',
            '> stdout'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_file_pair3(self) -> None:
        wf = FilePairsTestWF3()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def reads_joined = reads1 + \' \' + reads2',
            'def read1 = read1.simpleName != params.NULL_VALUE ? read1 : ${reads1}',
            'def read2 = read2.simpleName != params.NULL_VALUE ? read2 : ${reads2}',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${reads_joined}',
            '--reads-index-0 ${read1}',
            '--reads-index-1 ${read2}',
            '> stdout'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_file_pair_optional0(self) -> None:
        # name accession should be different?
        wf = FilePairsOptionalTestWF0()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def reads_joined = reads1.simpleName != params.NULL_VALUE ? reads1 + \' \' + reads2 : ""',
            'def reads1 = reads1.simpleName != params.NULL_VALUE ? reads1 : ""',
            'def reads2 = reads2.simpleName != params.NULL_VALUE ? reads2 : ""',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${reads_joined}',
            '> stdout'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_file_pair_optional1(self) -> None:
        # name accession should be different?
        wf = FilePairsOptionalTestWF1()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def reads_a_joined = reads_a1.simpleName != params.NULL_VALUE ? "--prefix ${reads_a1} ${reads_a2}" : ""',
            'def reads_a1 = reads_a1.simpleName != params.NULL_VALUE ? "--prefix " + reads_a1 : ""',
            'def reads_a2 = reads_a2.simpleName != params.NULL_VALUE ? "--prefix " + reads_a2 : ""',
            'def reads_b_joined = reads_b1.simpleName != params.NULL_VALUE ? "--prefixeach ${reads_b1} --prefixeach ${reads_b2}" : ""',
            'def reads_b1 = reads_b1.simpleName != params.NULL_VALUE ? "--prefixeach " + reads_b1 : ""',
            'def reads_b2 = reads_b2.simpleName != params.NULL_VALUE ? "--prefixeach " + reads_b2 : ""',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${reads_a_joined}',
            '${reads_b_joined}',
            '> stdout'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    # @unittest.skip('known bug but rare case')
    def test_file_pair_optional2(self) -> None:
        # name accession should be different?
        wf = FilePairsOptionalTestWF2()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def reads_joined = reads1.simpleName != params.NULL_VALUE ? reads1 + \' \' + reads2 : ""',
            'def reads1 = reads1.simpleName != params.NULL_VALUE ? reads1 : ""',
            'def reads2 = reads2.simpleName != params.NULL_VALUE ? reads2 : ""',
            'def read1 = read1.simpleName != params.NULL_VALUE ? read1 : ${reads1}',
            'def read2 = read2.simpleName != params.NULL_VALUE ? read2 : ${reads2}',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${reads_joined}',
            '--reads-index-0 ${read1}',
            '--reads-index-1 ${read2}',
            '> stdout',
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_file_pair_array(self) -> None:
        wf = FilePairsArrayTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def read_pairs = read_pairs_flat.collate(2, 2)',
            'def read_pairs_joined = read_pairs.collect{ it.join(\' \') }',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${read_pairs_joined}',
            '> stdout'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_file_pair_array_optional(self) -> None:
        wf = FilePairsArrayOptionalTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        
        print(process.get_string())
        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = [
            'def read_pairs = read_pairs_flat.collate(2, 2)',
            'def read_pairs_joined = read_pairs.collect{ it[0] != null ? it.join(\' \') : "" }',
        ]
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${read_pairs_joined}',
            '> stdout'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_filename_generated_tool(self):
        wf = FilenameTestWF1()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp3"]
        process = nextflow.generate.process.generate_process(step.tool)
        actual_process = simplify_file(process.get_string())
        print(process.get_string())
        expected_process = [
            'process FILENAME_GENERATED_TOOL {',
            'publishDir "${params.outdir}/filename_generated_tool"',
            'input:',
            'path file_inp',
            'path file_inp_optional, stageAs: \'file_inp_optional/*\'',
            'val inp',
            'val inp_optional',
            'output:',
            'val "*", emit: out',
            'script:',
            'def file_inp_optional = file_inp_optional.simpleName != params.NULL_VALUE ? file_inp_optional : ""',
            'def inp_optional = inp_optional != params.NULL_VALUE ? inp_optional : ""',
            '\"\"\"',
            'echo \\',
            '${inp} \\',
            '${inp_optional} \\',
            '${file_inp.baseName}.transformed.fnp \\',
            '${file_inp_optional.baseName}.optional.txt',
            '\"\"\"',
            '}'
        ]
        self.assertEqual(len(expected_process), len(actual_process))
        for ln in actual_process:
            self.assertIn(ln, expected_process)
   
    def test_filename_types1(self) -> None:
        wf = FilenameTestWF1()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())

        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = []
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${inp1}',
            '${inp2}',
            '> out'
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

    def test_filename_types2(self) -> None:
        wf = FilenameTestWF1()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp2"]
        process = nextflow.generate.process.generate_process(step.tool)
        print(process.get_string())

        actual_prescript = simplify_prescript(process.pre_script)
        expected_prescript = []
        self.assertEqual(len(actual_prescript), len(expected_prescript))
        for ln in expected_prescript:
            self.assertIn(ln, actual_prescript)

        actual_script = simplify_script(process.script)
        expected_script = [
            'echo',
            '${inp1}',
            '${inp1.baseName}.processed.txt',
            '> out',
        ]
        self.assertEqual(len(actual_script), len(expected_script))
        for ln in expected_script:
            self.assertIn(ln, actual_script)

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
        wf = do_preprocessing_workflow(wf, ignore_task_inputs=True)
        maintask, _ = translator.translate_workflow_internal(wf)



class TestTranslateToolInternal(unittest.TestCase):
    
    def setUp(self) -> None:
        reset_globals()

    def test_fastqc_tool(self) -> None:
        tool = FastqcTestTool()
        tool = tool.to_command_tool_builder()
        toolstr = translator.translate_tool_internal(tool)
    
    def test_bwamem_tool(self) -> None:
        tool = BwaMemTestTool()
        tool = tool.to_command_tool_builder()
        toolstr = translator.translate_tool_internal(tool)
    
    def test_gridss_tool(self) -> None:
        tool = GridssTestTool()
        tool = tool.to_command_tool_builder()
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
        wf = InputsPythonToolTestWF()
        self.mainstr, _, self.subtasks = translate(wf, dest_fmt='nextflow')

    def test_file_inputs(self) -> None:
        # File, String, Int input types
        process = [x[1] for x in self.subtasks if x[0] == 'modules/multi_types_input_python_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = {
            "path code_file",
            "path inp1",
            "val inp2",
            "val inp3",
        }
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in expected_lines:
            self.assertIn(ln, actual_lines)

    def test_generic_array_inputs(self) -> None:
        # Array(String) input type
        process = [x[1] for x in self.subtasks if x[0] == 'modules/join_array_python_test_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = {
            "path code_file",
            "val inp",
        }
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in expected_lines:
            self.assertIn(ln, actual_lines)
        
    def test_secondaries_inputs(self) -> None:
        # File (secondaries) input type
        process = [x[1] for x in self.subtasks if x[0] == 'modules/secondary_input_python_test_tool.nf'][0]
        actual_lines = _get_process_input_lines(process)
        expected_lines = {
            "path code_file",
            'path inp',
        }
        self.assertEqual(len(actual_lines), len(expected_lines))
        for ln in expected_lines:
            self.assertIn(ln, actual_lines)

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
        wf = InputsPythonToolTestWF()
        self.mainstr, _, self.subtasks = translate(wf, dest_fmt='nextflow', to_console=False)

    def test_format(self) -> None:
        process = [x[1] for x in self.subtasks if x[0] == 'modules/multi_types_input_python_tool.nf'][0]
        actual_lines = simplify_file(process)
        expected_lines = [
            'nextflow.enable.dsl=2',
            'process MULTI_TYPES_INPUT_PYTHON_TOOL {',
            'container "python:3.8.1"',
            'publishDir "${params.outdir}/multi_types_input_python_tool"',
            'input:',
            'path code_file',
            'path inp1',
            'val inp2',
            'val inp3',
            'output:',
            'val "${file("${task.workDir}/" + file("${task.workDir}/out_out").text.replace(\'"\', \'\'))}", emit: out',
            'exec:',
            'script:',
            '"""',
            '#!/usr/bin/env python',
            'from ${code_file.simpleName} import code_block',
            'import os',
            'import json',
            'result = code_block(',
            'inp1="${inp1}",',
            'inp2="${inp2}",',
            'inp3=${inp3}',
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
            'inp2="${inp2}",',
            'inp3=${inp3}',
            ')'
        }
        print(task.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

        # Array(String) input type
        task = self.processes['JoinArrayPythonTestTool']
        actual_script = task.script
        expected_lines = {
            'result = code_block(inp="${inp}".split(" "))',
        }
        print(task.get_string())
        for ln in expected_lines:
            self.assertIn(ln, actual_script)

        # File (secondaries) input type
        task = self.processes['SecondaryInputPythonTestTool']
        actual_script = task.script
        expected_lines = {
            'result = code_block(inp="${inp}")',
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
        actual_counts = trace.trace_entity_counts(src)
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
        actual_counts = trace.trace_entity_counts(src)
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
        actual_dtype = trace.trace_source_datatype(src)
        expected_dtype = File
        self.assertIsInstance(actual_dtype, expected_dtype)
    
    def test_trace_source_datatype2(self) -> None:
        step_id = 'stp2'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = trace.trace_source_datatype(src)
        expected_dtype = File
        self.assertIsInstance(actual_dtype, expected_dtype)
    
    def test_trace_source_datatype3(self) -> None:
        step_id = 'stp3'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = trace.trace_source_datatype(src)
        expected_dtype = String
        self.assertIsInstance(actual_dtype, expected_dtype)
    
    def test_trace_source_datatype4(self) -> None:
        step_id = 'stp4'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = trace.trace_source_datatype(src)
        expected_dtype = String
        self.assertIsInstance(actual_dtype, expected_dtype)
    
    def test_trace_source_datatype5(self) -> None:
        step_id = 'stp5'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_dtype = trace.trace_source_datatype(src)
        expected_dtype = String
        self.assertIsInstance(actual_dtype, expected_dtype)

    def test_trace_source_scatter6(self) -> None:
        step_id = 'stp6'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_scatter = trace.trace_source_scatter(src)
        expected_scatter = False
        self.assertEqual(actual_scatter, expected_scatter)

    def test_trace_source_scatter7(self) -> None:
        step_id = 'stp7'
        step = self.wf.step_nodes[step_id]
        src = step.sources['inp']
        actual_scatter = trace.trace_source_scatter(src)
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
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '')
        
        srctype, desttype, destscatter = FastaWithIndexes(), FastaWithIndexes(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '')
        
        srctype, desttype, destscatter = FastqGzPair(), FastqGzPair(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '')

        srctype, desttype, destscatter = Array(Bam()), Array(Bam()), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '')
        
        # plumbing required ---
        
        # BASETYPES
        # secondary, secondary
        srctype, desttype, destscatter = FastaWithIndexes(), FastaDict(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.map{ tuple -> [tuple[0], tuple[4]] }')
        
        # secondary, single
        srctype, desttype, destscatter = FastaWithIndexes(), Fasta(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.map{ tuple -> tuple[0] }')

        # ARRAYS      
        # single, single  
        srctype, desttype, destscatter = Array(Bam()), Bam(), True
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.flatten()')
        
        srctype, desttype, destscatter = Array(Fasta()), Fasta(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.flatten().first()')
        
        srctype, desttype, destscatter = Fasta(), Array(Fasta()), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.toList()')

        # any, secondary (always ends with .flatten().toList())
        srctype, desttype, destscatter = FastaWithIndexes(), Array(FastaWithIndexes()), True
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.flatten().toList()')
        
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Array(FastaWithIndexes()), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.flatten().toList()')
        
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Array(FastaWithIndexes()), True
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.flatten().toList()')

        # secondary, secondary 
        srctype, desttype, destscatter = Array(FastaWithIndexes()), FastaWithIndexes(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.flatten().collate( 8 ).first()')

        srctype, desttype, destscatter = FastaWithIndexes(), Array(FastaDict()), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.map{ tuple -> [tuple[0], tuple[4]] }.flatten().toList()')
        
        # secondary, single 
        srctype, desttype, destscatter = Array(FastaWithIndexes()), Fasta(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.map{ tuple -> tuple[0] }.flatten().first()')

        # file pairs
        srctype, desttype, destscatter = Array(FastqGzPair()), FastqGzPair(), True
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.flatten().collate( 2 )')
        
        srctype, desttype, destscatter = Array(FastqGzPair()), FastqGzPair(), False
        plumbing = gen_datatype_mismatch_plumbing(srctype, desttype, destscatter, is_connection=True)
        self.assertEqual(plumbing, '.flatten().collate( 2 ).first()')

        
        


class TestPlumbingTypeMismatch(unittest.TestCase):
    """
    This test group checks we can handle array / single datatype mismatches
    in a workflow. This occurs in wgsgermline. 
    """
    def setUp(self) -> None:
        reset_globals()
        wf = PlumbingTypeMismatchTestWF()
        self.mainstr, _, _ = translate(wf, dest_fmt='nextflow')
        print(self.mainstr)

    def test_secondary_single_mismatch(self):
        actual = _get_task_call_lines(self.mainstr, 'BAMBAI_TO_BAM')
        expected = ['in_bam_bai[0]']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_single_array_mismatch(self):
        actual = _get_task_call_lines(self.mainstr, 'BAMBAI_TO_BAM_ARRAY')
        expected = ['in_bam_bai[0].toList()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_secondary_mismatch(self):
        actual = _get_task_call_lines(self.mainstr, 'FASTAWITHINDEXES_TO_FASTADICT')
        expected = ['in_fasta_with_indexes[0, 4]']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_secondary_array_mismatch(self):
        actual = _get_task_call_lines(self.mainstr, 'FASTAWITHINDEXES_TO_FASTADICT_ARRAY')
        expected = ['in_fasta_with_indexes[0, 4].flatten().toList()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_array_to_single(self):
        actual = _get_task_call_lines(self.mainstr, 'ARRAY_TO_SINGLE')
        expected = ['ch_in_file_array.flatten().first()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_single_to_array(self):
        actual = _get_task_call_lines(self.mainstr, 'SINGLE_TO_ARRAY')
        expected = ['in_file.toList()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_array_to_secondary(self):
        actual = _get_task_call_lines(self.mainstr, 'SECONDARY_ARRAY_TO_SECONDARY')
        expected = ['ch_in_bam_bai_array.flatten().collate( 2 ).first()']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)
    
    def test_secondary_array_to_secondary_array(self):
        actual = _get_task_call_lines(self.mainstr, 'SECONDARY_ARRAY_TO_SECONDARY_ARRAY')
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
        reset_globals()
        wf = ComprehensiveScatterTestWF()
        self.mainstr, _, _ = translate(wf, dest_fmt='nextflow')
        print(self.mainstr)

    def test_scatter_to_scatter(self):
        actual = _get_task_call_lines(self.mainstr, 'SCATTER_TO_SCATTER')
        expected = ['PRESTEP1.out.out']
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_scatter_to_array(self):
        actual = _get_task_call_lines(self.mainstr, 'SCATTER_TO_ARRAY')
        expected = ['PRESTEP1.out.out.toList()']
        self.assertEqual(actual, expected)

    def test_array_to_scatter1(self):
        actual = _get_task_call_lines(self.mainstr, 'PRESTEP1')
        expected = ['ch_in_file_array.flatten()']
        self.assertEqual(actual, expected)
    
    def test_array_to_scatter2(self):
        actual = _get_task_call_lines(self.mainstr, 'ARRAY_TO_SCATTER')
        expected = ['PRESTEP2.out.out.flatten()']
        self.assertEqual(actual, expected)

    def test_scatter_secondary_to_scatter_secondary(self):
        actual = _get_task_call_lines(self.mainstr, 'SCATTER_SECONDARY_TO_SCATTER_SECONDARY')
        expected = ['PRESTEP3.out.out']
        self.assertEqual(actual, expected)

    def test_scatter_secondary_to_secondary_array(self):
        actual = _get_task_call_lines(self.mainstr, 'SCATTER_SECONDARY_TO_SECONDARY_ARRAY')
        expected = ['PRESTEP3.out.out.flatten().toList()']
        self.assertEqual(actual, expected)

    def test_secondary_array_to_scatter_secondary(self):
        actual = _get_task_call_lines(self.mainstr, 'SECONDARY_ARRAY_TO_SCATTER_SECONDARY')
        expected = ['ch_in_bam_bai_array.flatten().collate( 2 )']
        self.assertEqual(actual, expected)

    # TODO use in future
    @unittest.skip('reimplement')
    def test_scatter_cross(self) -> None:
        wf = ScatterCrossTestWF()
        wf = do_preprocessing_workflow(wf)
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
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)
        actual = _get_task_call_lines(mainstr, 'STP0')
        expected = [
            "params.multi_types_input_python_tool.code_file,",
            "in_file,",
            "params.in_str,",
            "params.in_int",
        ]
        self.assertListEqual(expected, actual)

    def test_workflow_inputs(self):
        wf = CallWFInputTestWF()
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)

        actual = _get_task_call_lines(mainstr, 'STP1')
        expected = [
            "in_file,",
            "params.in_bool,",
            "params.in_bool,",
            "params.in_str,",
            "params.in_int,",
            "params.in_int",
        ]
        self.assertListEqual(expected, actual)
        
        actual = _get_task_call_lines(mainstr, 'STP2')
        expected = [
            "params.in_bool,",
            "params.in_bool,",
            "params.in_str,",
            "params.in_str",
        ]
        self.assertListEqual(expected, actual)
        
        actual = _get_task_call_lines(mainstr, 'STP3')
        expected = [
            "ch_in_file_arr,",
            "params.in_str_arr,",
            "params.in_str_arr,",
            "params.in_int_arr,",
            "params.in_str_arr,",
            "params.in_int_arr",
        ]
        self.assertListEqual(expected, actual)
        
        actual = _get_task_call_lines(mainstr, 'STP4')
        expected = [
            "in_file_arr_opt,",
            "params.in_str_arr_opt,",
            "params.in_str_arr_opt",
        ]
        self.assertListEqual(expected, actual)
    
    def test_workflow_inputs_duplicates(self):
        wf = CallDuplicateUseageWFInputTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp1'])
        expected = [
            "ch_in_file",
            "params.in_bool",
            "params.in_bool",
            "params.in_str",
            "params.in_int",
            "params.in_int",
        ]
        self.assertListEqual(expected, actual)
        
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp1_supp'])
        expected = [
            "ch_in_file",
            "true",
            "true",
            '"hi"',
            "10",
            "10",
        ]
        self.assertListEqual(expected, actual)
        
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = [
            'params.in_bool',
            'params.in_bool',
            'params.in_str',
            'params.in_str',
        ]
        self.assertListEqual(expected, actual)
        
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2_supp'])
        # this commented out section is a TODO. 
        # expected = [
        #     'params.NULL_VALUE',
        #     'params.NULL_VALUE',
        #     'params.NULL_VALUE',
        #     'params.NULL_VALUE',
        # ]
        expected = [
            'params.stp2_supp_flag_false',
            'params.stp2_supp_flag_true',
            'params.stp2_supp_opt_optional',
            'params.stp2_supp_pos_optional',
        ]
        self.assertListEqual(expected, actual)
        
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp3'])
        expected = [
            "ch_in_file_arr",
            "params.in_str_arr",
            "params.in_str_arr",
            "params.in_int_arr",
            "params.in_str_arr",
            "params.in_int_arr",
        ]
        self.assertListEqual(expected, actual)
        
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp3_supp'])
        expected = [
            "ch_in_file_arr",
            "['hi', 'there']",
            "['hi', 'there']",
            "[1, 2, 3]",
            "['hi', 'there']",
            "[1, 2, 3]",
        ]
        self.assertListEqual(expected, actual)
        
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp4'])
        expected = [
            "in_file_arr_opt",
            "params.in_str_arr_opt",
            "params.in_str_arr_opt",
        ]
        self.assertListEqual(expected, actual)
        
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp4_supp'])
        expected = [
            'stp4_supp_pos_optional_arr',
            'params.stp4_supp_opt_optional_arr',
            'params.stp4_supp_opt_optional_arr_prefixeach',
        ]
        self.assertListEqual(expected, actual)
        
    # static step inputs
    def test_static_inputs(self):
        wf = CallStaticTestWF()
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)
        
        actual = _get_task_call_lines(mainstr, 'STP1')
        expected = [
            'in_file,',
            'true,',
            'false,',
            '"static",',
            '100,',
            '100',
        ]
        self.assertListEqual(actual, expected)
        
        actual = _get_task_call_lines(mainstr, 'STP2')
        expected = [
            'params.stp2_flag_false,',
            'params.stp2_flag_true,',
            'params.stp2_opt_optional,',
            'params.stp2_pos_optional',
        ]
        self.assertListEqual(actual, expected)
        
        actual = _get_task_call_lines(mainstr, 'STP3')
        expected = [
            "ch_in_file_arr,",
            "['hi', 'there'],",
            "['hi', 'there'],",
            "[1, 2, 3],",
            "['hi', 'there'],",
            "[1, 2, 3]",
        ]
        self.assertListEqual(actual, expected)
        
        actual = _get_task_call_lines(mainstr, 'STP4')
        expected = [
            'stp4_pos_optional_arr,',
            'params.stp4_opt_optional_arr,',
            'params.stp4_opt_optional_arr_prefixeach',
        ]
        self.assertListEqual(actual, expected)

    # connections
    def test_connections_files(self) -> None:
        # singles
        wf = CallConnectionsTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = ["STP1.out.out"]
        self.assertListEqual(actual, expected)
        
        # arrays
        wf = CallArrayConnectionsTestWF()
        actual = _gen_call_lines_local(wf, step=wf.step_nodes['stp2'])
        expected = ["STP1.out.out"]
        self.assertListEqual(actual, expected)

    def test_filename_types(self) -> None:
        wf = FilenameTestWF1()
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)
        
        actual = _get_task_call_lines(mainstr, 'STP1')
        expected = ["in_file,", "params.in_str"]
        self.assertListEqual(actual, expected)
        
        actual = _get_task_call_lines(mainstr, 'STP2')
        expected = ["in_file"]
        self.assertListEqual(actual, expected)

    def test_subworkflow(self) -> None:
        wf = DataSourceTestWF()
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)
        
        actual = _get_task_call_lines(mainstr, 'STP1')
        expected = [
            "in_file1,",
            "in_file_opt1,",
            "params.in_str1,",
            "params.in_str_opt1",
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
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)
        actual = _get_task_call_lines(mainstr, 'STP2')
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
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)
        actual = _get_task_call_lines(mainstr, 'PRINT')
        expected = [
            "[params.mystring, GET_STRING.out.out].find{ it != null }"
        ]
        self.assertEqual(len(actual), len(expected))
        for arg in expected:
            self.assertIn(arg, actual)

    def test_with_expression(self):
        wf = StepInputExpressionTestWF()
        mainstr, _, _ = translate(wf, dest_fmt='nextflow', to_console=False)
        print(mainstr)
        actual = _get_task_call_lines(mainstr, 'PRINT')
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
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        process = nextflow.generate.process.generate_process(step.tool)
        self.prescript = process.pre_script
        self.script = process.script
        print(process.get_string())

    # PROCESS RELATED
    @unittest.skip('filenames are scuffed')
    def test_filename_generated(self) -> None:
        self.assertIn("--filenameGen generated.gz", self.script)
    @unittest.skip('filenames are scuffed')
    def test_filename_reference(self) -> None:
        self.assertIn("--filenameRef ${in_file.simpleName}.fastq.gz", self.script)

    def test_input_selector_process_input(self) -> None:
        self.assertIn("--inputSelectorProcess ${in_file}", self.script)

    def test_input_selector_param_input(self) -> None:
        self.assertIn("--inputSelectorParam ${in_str}", self.script)
    
    def test_input_selector_array(self) -> None:
        self.assertIn("--InputSelectorArray ${in_file_arr_joined}", self.script)
       
    def test_list(self) -> None:
        self.assertIn('--list "1 2 3 4 5"', self.script) # this is correct
       
    def test_two_value_operator(self) -> None:
        self.assertIn("--TwoValueOperator ${in_file + \".gz\"}", self.script)
    
    def test_first_operator(self) -> None:
        self.assertIn("--FirstOperator ${[in_str, []].find{ it != null }}", self.script)
    
    def test_index_operator_array(self) -> None:
        self.assertIn("--IndexOperatorArray ${in_file_arr[0]}", self.script)
    
    def test_index_operator_secondaries(self) -> None:
        self.assertIn("--IndexOperatorSecondariesBam ${in_bam_bai}", self.script)
        self.assertIn("--IndexOperatorSecondariesBai ${in_bam_bai}", self.script)
    
    def test_index_operator_secondaries_array(self) -> None:
        print(self.script)
        self.assertIn("--IndexOperatorArraySecondariesBams ${in_bam_bai_arr[0]}", self.script)
        self.assertIn("--IndexOperatorArraySecondariesBais ${in_bam_bai_arr[1]}", self.script)



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
        wf = CallConnectionsTestWF()
        wf = do_preprocessing_workflow(wf)
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
        wf = CallConnectionsTestWF()
        wf = do_preprocessing_workflow(wf)
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
        wf = CallConnectionsTestWF()
        wf = do_preprocessing_workflow(wf)
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
        wf = do_preprocessing_workflow(wf)
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
        wf = do_preprocessing_workflow(wf)
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
        wf = do_preprocessing_workflow(wf)
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
    
    def test_basic(self):
        settings.translate.nextflow.ENTITY = 'tool'
        tool = BasicTestTool()
        tool = do_preprocessing_tool(tool)
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter("no format")
        actual = nextflow.unwrap_expression(
            val=sf,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual("no format", actual)

    def test_string(self):
        settings.translate.nextflow.ENTITY = 'tool'
        tool = BasicTestTool()
        tool = do_preprocessing_tool(tool)
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter("there's a {str_arg} arg", str_arg="string")
        actual = nextflow.unwrap_expression(
            val=sf,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True
        )
        self.assertEqual("there\'s a string arg", actual)
    
    def test_input_selector_process_input(self):
        settings.translate.nextflow.ENTITY = 'tool'
        tool = BasicTestTool()
        tool = do_preprocessing_tool(tool)
        variable_manager = init_variable_manager_for_task(tool)
        sf = StringFormatter("an input {arg}", arg=InputSelector("testtool"))
        actual = nextflow.unwrap_expression(
            val=sf,
            context='process_script',
            variable_manager=variable_manager,
            tool=tool, 
            in_shell_script=True,
        )
        self.assertEqual("an input ${testtool}", actual)
    
    def test_input_selector_param_input(self):
        wf = StringFormatterTestWF()
        wf = do_preprocessing_workflow(wf)
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
        expected = "an input ${compression_level}"
        self.assertEqual(actual, expected)

    def test_two_params(self):
        settings.translate.nextflow.ENTITY = 'tool'
        tool = InputQualityTestTool()
        tool = do_preprocessing_tool(tool)
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
        self.assertEqual("${user}:${static}", actual)

    def test_escaped_characters(self):
        settings.translate.nextflow.ENTITY = 'tool'
        tool = InputQualityTestTool()
        tool = do_preprocessing_tool(tool)
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
        print(actual_shell)
        self.assertEqual("user\\tstatic", actual_scripting)
        self.assertEqual("${user}\\\\t${static}", actual_shell)

    def test_expression(self):
        settings.translate.nextflow.ENTITY = 'tool'
        tool = BasicTestTool()
        tool = do_preprocessing_tool(tool)
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
        self.assertEqual("${testtool}:${array_inp.join(\";\")}", res)
    
    def test_argument(self) -> None:
        wf = StringFormatterTestWF()
        wf = do_preprocessing_workflow(wf)
        step = wf.step_nodes["stp1"]
        vmanager = init_variable_manager_for_task(step.tool)
        arg = step.tool.arguments()[0]
        actual = nextflow.generate.process.script.common.eval_cmdline_targ(
            arg=arg,
            tool=step.tool,
            vmanager=vmanager,
            shell_quote=True
        )
        expected = '--java-options "-Xmx${8 * 3 / 4}G ${compression_level ? "-Dsamjdk.compress_level=" + compression_level : ""} ${[java_options, []].find{ it != null }.join(" ")}"'
        self.assertEqual(actual, expected)


class TestOrdering(unittest.TestCase):

    def setUp(self) -> None:
        reset_globals()
        wf = OrderingTestWF()
        self.mainstr, _, self.subtasks = translate(wf, dest_fmt='nextflow', to_console=False)

    def test_process_call(self) -> None:
        # from workflow inputs
        actual = _get_task_call_lines(self.mainstr, 'STP1')
        expected = [
            "ch_in_fastq_array,",
            "in_fastq,",
            "in_file,",
            "params.in_int_array,",
            "params.in_int,",
            "params.in_str",
        ]
        self.assertEqual(expected, actual)
        
        # from process outputs
        actual = _get_task_call_lines(self.mainstr, 'STP3')
        expected = [
            "STP1.out.outFastqArray,",
            "STP1.out.outFastq,",
            "STP1.out.outFile,",
            "STP1.out.outIntArray,",
            "STP1.out.outInt,",
            "STP1.out.outStr",
        ]
        self.assertEqual(expected, actual)

        # from subworkflow outputs
        actual = _get_task_call_lines(self.mainstr, 'STP5')
        expected = [
            "STP2.out.outFastqArray,",
            "STP2.out.outFastq,",
            "STP2.out.outFile,",
            "STP2.out.outIntArray,",
            "STP2.out.outInt,",
            "STP2.out.outStr",
        ]
        self.assertEqual(expected, actual)

    def test_process_inputs(self) -> None:
        process = [x[1] for x in self.subtasks if x[0] == 'modules/multi_type_test_tool.nf'][0]
        print(process)
        actual_lines = _get_process_input_lines(process)
        expected_lines = [
            "path in_fastq_array",
            "path in_fastq",
            "path in_file",
            "val in_int_array",
            "val in_int",
            "val in_str",
        ]
        self.assertEqual(actual_lines, expected_lines)

    def test_process_directives(self) -> None:
        wf = DirectivesTestWF()
        _, _, subtasks = translate(wf, dest_fmt='nextflow', to_console=False)
        process = [x[1] for x in subtasks if x[0] == 'modules/resources_test_tool.nf'][0]
        print(process)
        actual_order = _get_process_directive_lines(process)
        expected_order = [
            'container "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"',
            'publishDir "${params.outdir}/resources_test_tool"',
            'cpus "${params.resources_test_tool.cpus}"',
            'disk "${params.resources_test_tool.disk}"',
            'memory "${params.resources_test_tool.memory}"',
            'time "${params.resources_test_tool.time}"',
        ]
        for actual, expected in zip(actual_order, expected_order):
            self.assertEqual(actual, expected)

    def test_subworkflow_call(self) -> None:
        actual = _get_task_call_lines(self.mainstr, 'STP2')
        expected = [
            "in_fastq,",
            "in_file,",
            "ch_in_fastq_array,",
            "params.in_int,",
            "params.in_int_array,",
            "params.in_str",
        ]
        self.assertEqual(expected, actual)
        
        # from process outputs
        actual = _get_task_call_lines(self.mainstr, 'STP4')
        expected = [
            "STP1.out.outFastq,",
            "STP1.out.outFile,",
            "STP1.out.outFastqArray,",
            "STP1.out.outInt,",
            "STP1.out.outIntArray,",
            "STP1.out.outStr",
        ]
        self.assertEqual(expected, actual)

        # from subworkflow outputs
        actual = _get_task_call_lines(self.mainstr, 'STP6')
        expected = [
            "STP2.out.outFastq,",
            "STP2.out.outFile,",
            "STP2.out.outFastqArray,",
            "STP2.out.outInt,",
            "STP2.out.outIntArray,",
            "STP2.out.outStr",
        ]
        self.assertEqual(expected, actual)
    
    def test_subworkflow_inputs(self) -> None:
        subwf = [x[1] for x in self.subtasks if x[0] == 'subworkflows/multi_type_test_wf.nf'][0]
        actual_inputs = _lines_within_section(subwf, 'take:', 'main:')
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
    
