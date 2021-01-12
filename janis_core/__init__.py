"""
Welcome to the source code for Janis!

Janis is a framework creating specialised, simple workflow definitions that are
transpiled to Common Workflow Language or Workflow Definition Language.

Below you'll find the classes you can access by importing Janis:

>>  import janis as j

Some noteworthy classes are Workflow, CommandTool and the janis.translations module

Some terminology:
    - Edge:
        - An edge may exist between 2 nodes, it represents a dependency
        - Every edge will have a source_map which is indexed by the tag of the input on the finish node
        - This source_map value is singular or a list, as you can connect multiple sources to one input



_Proudly made on Planet Earth._

"""

# PEP396:  https://www.python.org/dev/peps/pep-0396/
from janis_core.__meta__ import __version__
from janis_core.toolbox.toolbox import JanisShed
import janis_core.toolbox.entrypoints as entrypoints

# Tools
from janis_core.tool.tool import Tool, ToolType, TOutput, TInput
from janis_core.workflow.workflow import (
    Workflow,
    WorkflowBuilder,
    WorkflowBase,
    DynamicWorkflow,
)
from janis_core.tool.commandtool import (
    CommandTool,
    CommandToolBuilder,
    ToolArgument,
    ToolInput,
    ToolOutput,
)
from janis_core.code.pythontool import PythonTool, CodeTool

# Types

from janis_core.types.data_types import DataType
from janis_core.types.common_data_types import (
    Boolean,
    String,
    Int,
    Float,
    Double,
    File,
    Directory,
    Array,
    Filename,
    Stdout,
    Stderr,
)
from janis_core.operators import (
    InputSelector,
    WildcardSelector,
    MemorySelector,
    CpuSelector,
    StringFormatter,
    standard,
    logical,
)

# Misc
from janis_core.tool.documentation import *
from janis_core.utils.logger import Logger, LogLevel
from janis_core.translations import SupportedTranslation
from janis_core.utils.scatter import ScatterDescription, ScatterMethod, ScatterMethods
from janis_core.hints import CaptureType, Engine, HINTS, Hint, HintEnum, HintArray
from janis_core.utils import get_value_for_hints_and_ordered_resource_tuple
from janis_core.utils.metadata import Metadata, WorkflowMetadata, ToolMetadata
from janis_core.utils.secondary import apply_secondary_file_format_to_filename
from janis_core.transformation import JanisTransformation, JanisTransformationGraph
