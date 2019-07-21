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

import sys
import pkg_resources

print("NAME: " + __name__)

# for entrypoint in pkg_resources.iter_entry_points(group="janis.extension"):
#     try:
#         m = entrypoint.load()
#         globals()[entrypoint.name] = entrypoint.load
#         sys.modules["janis." + entrypoint.name] = entrypoint.dist.location
#     except ImportError:
#         print("Import Error: " + entrypoint.name)
#         continue

# PEP396:  https://www.python.org/dev/peps/pep-0396/
from core.__meta__ import __version__

from core.hints import CaptureType, Engine, HINTS, Hint, HintEnum, HintArray
from core.tool.commandtool import CommandTool
from core.tool.tool import Tool, ToolArgument, ToolInput, ToolOutput
from core.translations import SupportedTranslations
from core.types import (
    InputSelector,
    WildcardSelector,
    MemorySelector,
    CpuSelector,
    StringFormatter,
)
from core.types.common_data_types import (
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
)
from core.types.data_types import DataType
from core.unix import *
from core.utils.logger import Logger, LogLevel
from core.utils.metadata import Metadata, WorkflowMetadata, ToolMetadata
from core.workflow.input import Input
from core.workflow.output import Output
from core.workflow.step import Step
from core.workflow.workflow import Workflow
