import re
from abc import ABC, abstractmethod
from typing import Optional, List, Dict, Any, Union

from janis_core.types import (
    ParseableType,
    Selector,
    Array,
    get_instantiated_type,
    DataType,
)
from janis_core.utils.logger import Logger
from janis_core.utils.metadata import Metadata
from janis_core.utils.validators import Validators

ToolType = str


class ToolTypes:
    Workflow: ToolType = "workflow"
    CommandTool: ToolType = "command-tool"
    ExpressionTool: ToolType = "expression-tool"


class TInput(object):
    def __init__(self, tag: str, intype: DataType, default=None):
        self.tag = tag
        self.intype = intype
        self.default = default

    def id(self):
        return self.tag


class TOutput(object):
    def __init__(self, tag, outtype):
        self.tag = tag
        self.outtype = outtype

    def id(self):
        return self.tag


class Tool(ABC, object):
    """
    One of Workflow, CommandLineTool, ExpressionTool* (* unimplemented)
    """

    def __init__(self, metadata_class=Metadata, **connections):
        """
        :param metadata_class:
        :param connections:
        """

        self.metadata: metadata_class = metadata_class()
        meta = self.bind_metadata()
        if meta:
            self.metadata = meta

        self.connections = connections

    @classmethod
    @abstractmethod
    def type(cls) -> ToolType:
        raise Exception(f"'{cls}' must implement type() method")

    @abstractmethod
    def id(self) -> str:
        raise Exception("Must implement id() method")

    @staticmethod
    def tool_module():
        return None

    @staticmethod
    def tool_provider():
        return None

    @abstractmethod
    def tool_inputs(self) -> List[TInput]:
        raise Exception("Must implement inputs() method")

    @abstractmethod
    def tool_outputs(self) -> List[TOutput]:
        raise Exception("Must implement outputs() method")

    def inputs_map(self) -> Dict[str, TInput]:
        return {inp.tag: inp for inp in self.tool_inputs()}

    def outputs_map(self) -> Dict[str, TOutput]:
        return {outp.tag: outp for outp in self.tool_outputs()}

    def friendly_name(self) -> Optional[str]:
        """
        Overriding this method is not required UNLESS you distribute your tool.
        Generating the docs will fail if your tool does not provide a name.

        :return: A friendly name of your tool
        """
        return None

    def all_input_keys(self) -> List[str]:
        return [t.id() for t in self.tool_inputs()]

    @abstractmethod
    def generate_inputs_override(
        self, additional_inputs=None, with_resource_overrides=False, hints=None
    ):
        pass

    def __call__(self, **connections):
        self.connections = connections
        return self

    @staticmethod
    @abstractmethod
    def version():
        return None

    def doc(self) -> Optional[str]:
        return None

    @abstractmethod
    def translate(
        self, translation: str, with_docker=True, with_resource_overrides=False
    ):
        raise Exception("Subclass must provide implementation for 'translate()' method")

    def bind_metadata(self):
        """
        A convenient place to add metadata about the tool. You are guaranteed that self.metadata will exist.
        It's possible to return a new instance of the ToolMetadata / WorkflowMetadata which will be rebound.
        This is usually called after the initialiser, though it may be called multiple times.
        :return:
        """
        return self.metadata
