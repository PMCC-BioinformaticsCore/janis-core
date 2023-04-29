

from .basics import BasicIOTestWF
from .basics import WildcardSelectorOutputTestWF
from .basics import InputSelectorTestWF
from .basics import DirectivesTestWF
from .basics import ArrayIOTestWF
from .basics import ArrayIOExtrasTestWF

from .components import ComponentsMandatoryTestWF
from .components import ComponentsOptionalTestWF
from .components import ComponentsMandatoryArrayTestWF
from .components import ComponentsOptionalArrayTestWF

from .components import StepInputsWFInputTestWF
from .components import StepInputsStaticTestWF
from .components import StepInputsPartialStaticTestWF
from .components import StepInputsMinimalTestWF
from .components import StepConnectionsTestWF
from .components import ArrayStepInputsTestWF
from .components import ArrayStepConnectionsTestWF

from .scatter import BasicScatterTestWF
from .scatter import ChainedScatterTestWF
from .scatter import ScatterDotTestWF
from .scatter import ScatterCrossTestWF
from .scatter import ComprehensiveScatterTestWF

from .secondaries import SecondariesTestWF
from .secondaries import SecondariesTestTool
from .filenames import FilenameTestWF1
from .filenames import FilenameTestWF2

from .combos import ScatterSecondariesTestWF
from .outputs import OutputCollectionTestWF
from .unwrap import UnwrapTestWF
from .subworkflow import SubworkflowTestWF
from .subworkflow2 import Subworkflow2TestWF
from .subworkflow2 import Subworkflow3TestWF
from .data_sources import DataSourceTestWF
from .naming import NamingTestWF

from .file_pairs import AllFilePairsTestWF
from .file_pairs import FilePairsTestWF0
from .file_pairs import FilePairsTestWF1
from .file_pairs import FilePairsTestWF2
from .file_pairs import FilePairsTestWF3
from .file_pairs import FilePairsOptionalTestWF0
from .file_pairs import FilePairsOptionalTestWF1
from .file_pairs import FilePairsOptionalTestWF2
from .file_pairs import FilePairsOptionalTestWF3
from .file_pairs import FilePairsArrayTestWF
from .file_pairs import FilePairsArrayOptionalTestWF

from .files_directories_to_create import FilesDirectoriesToCreateTestWF
from .process_inputs import ProcessInputsTestWF
from .ordering import OrderingTestWF
from .plumbing_edge_cases import PlumbingEdgeCaseTestWF
from .plumbing_type_mismatch import PlumbingTypeMismatchTestWF
from .trace import EntityTraceTestWF
from .optional_types import OptionalInputTypesTestWF
from .mandatory_types import MandatoryInputTypesTestWF
from .duplicate_tasks import DuplicateTasksTestWF
from .minimal_task_inputs import MinimalTaskInputsTestWF1
from .minimal_task_inputs import MinimalTaskInputsTestWF2
from .minimal_task_inputs import MinimalTaskInputsTestWF3
from .minimal_task_inputs import MinimalTaskInputsTestWF4
from .minimal_task_inputs import MinimalTaskInputsTestWF5
from .minimal_task_inputs import MinimalTaskInputsTestWF6
from .all_input_types import AllInputTypesTestWF


from .assembly import w as AssemblyTestWF
from .string_formatter import StringFormatterTestWF

from .additional_features import (
    StepInputExpressionTestWF,
    ConditionStepTestWF,
    AliasSelectorTestWF,
    ArraysOfSecondaryFilesOutputsTestWF,
    ForEachTestWF,
    IndexOperatorTestWF,
)

from .codetools import (
    InputsPythonToolTestWF,
    OutputsPythonToolTestWF
)

