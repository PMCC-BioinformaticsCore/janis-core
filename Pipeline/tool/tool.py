from abc import ABC, abstractmethod
from typing import Optional, List, Dict, Any
import re

from Pipeline.utils.logger import Logger
from Pipeline.types.data_types import DataType, NativeTypes
from Pipeline.translations.cwl.cwl import Cwl

ToolType = str


class ToolTypes:
    Workflow: ToolType = "workflow"
    CommandTool: ToolType = "command-tool"
    ExpressionTool: ToolType = "expression-tool"


class ToolArgument:
    expr_pattern = "\$\(.*\)"

    def __init__(self, value: str, prefix: Optional[str] = None, position: Optional[int] = 0,
                 separate_value_from_prefix=True, doc: Optional[str]=None):
        self.prefix: Optional[str] = prefix
        self.value = value
        self.position: Optional[int] = position
        self.is_expression = re.match(self.expr_pattern, self.value) is not None
        self.separate_value_from_prefix = separate_value_from_prefix
        self.doc = doc

        if self.prefix and not self.separate_value_from_prefix and not self.prefix.endswith("="):
            # I don't really know what this means.
            Logger.warn(f"Argument ({self.prefix} {self.value}) is not separating and did not end with ='")

    def cwl(self):
        from cwlgen import CommandLineBinding
        return CommandLineBinding(
            # load_contents=False,
            position=self.position,
            prefix=self.prefix,
            separate=self.separate_value_from_prefix,
            # item_separator=None,
            value_from=self.value,
            # shell_quote=True,
        )

    def wdl(self):
        return (self.prefix if self.prefix is not None else "") \
               + (" " if self.separate_value_from_prefix else "") \
               + (self.value if self.value is not None else "")


class ToolInput(ToolArgument):
    def __init__(self, tag: str, input_type: DataType, position: Optional[int] = None, prefix: Optional[str] = None,
                 separate_value_from_prefix: bool = True, default: Any = None, doc: Optional[str]=None):
        """
        :param tag: tag for input, what the yml will reference (eg: input1: path/to/file)
        :param input_type:
        """
        super().__init__("", prefix, position, separate_value_from_prefix)

        if default is not None:
            input_type.optional = True

        self.tag: str = tag
        self.input_type: DataType = input_type
        self.optional = self.input_type.optional
        self.default = default
        self.doc = doc

    def cwl(self):
        from cwlgen import CommandInputParameter, CommandLineBinding

        default = self.default if self.default is not None else self.input_type.default()

        return CommandInputParameter(
            param_id=self.tag,
            label=self.tag,
            secondary_files=self.input_type.secondary_files(),
            # streamable=False,
            doc=self.doc,
            default=default,
            input_binding=CommandLineBinding(
                # load_contents=self.load_contents,
                position=self.position,
                prefix=self.prefix,
                separate=self.separate_value_from_prefix
                # item_separator=self.item_separator,
                # value_from=self.value_from,
                # shell_quote=self.shell_quote
            ),
            param_type=self.input_type.cwl2_type()
        )


class ToolOutput:
    def __init__(self, tag: str, output_type: DataType, glob: Optional[str] = None, doc: Optional[str]=None):
        self.tag = tag
        self.output_type: DataType = output_type
        self.glob = glob
        self.doc = doc

    def cwl(self):
        from cwlgen import CommandOutputParameter, CommandOutputBinding
        return CommandOutputParameter(
            param_id=self.tag,
            label=self.tag,
            secondary_files=self.output_type.secondary_files(),
            # param_format=None,
            # streamable=False,
            doc=self.doc,
            output_binding=CommandOutputBinding(
                glob=self.glob,
                # load_contents=False,
                # output_eval=None
            ),
            param_type=self.output_type.cwl2_type()
        )


class Tool(ABC, object):
    """
    One of Workflow, CommandLineTool, ExpressionTool* (* unimplemented)
    """

    @classmethod
    def type(cls) -> ToolType:
        print(cls)
        raise Exception("Must implement type() method")

    @abstractmethod
    def id(self) -> str:
        raise Exception("Must implement id() method")

    @abstractmethod
    def inputs(self) -> List[ToolInput]:
        raise Exception("Must implement inputs() method")

    @abstractmethod
    def outputs(self) -> List[ToolOutput]:
        raise Exception("Must implement outputs() method")

    def inputs_map(self) -> Dict[str, ToolInput]:
        return {inp.tag: inp for inp in self.inputs()}

    def outputs_map(self) -> Dict[str, ToolOutput]:
        return {outp.tag: outp for outp in self.outputs()}

    @abstractmethod
    def cwl(self, with_docker=True) -> Dict[str, Any]:
        raise Exception("Must implement cwl() method")

    @staticmethod
    def doc() -> Optional[str]:
        return None

    def help(self):
        import inspect
        path = inspect.getfile(self.__class__)

        ins = sorted(self.inputs(), key=lambda i: i.position)

        def input_format(t: ToolInput):
            prefix_with_space = ""
            if t.prefix is not None:
                prefix_with_space = (t.prefix + ": ") if t.separate_value_from_prefix else t.prefix
            return f"\t\t{t.tag} ({prefix_with_space}{t.input_type.id()}" \
                f"{('=' + str(t.default)) if t.default is not None else ''}): {'' if t.doc is None else t.doc}"

        output_format = lambda t: f"\t\t{t.tag} ({t.output_type.id()}): {'' if t.doc is None else t.doc}"

        requiredInputs = "\n".join(input_format(x) for x in ins if not x.optional)
        optionalInputs = "\n".join(input_format(x) for x in ins if x.optional)
        outputs = "\n".join(output_format(o) for o in self.outputs())

        return f"""
Pipeline tool: {path} ({self.id()})
NAME
    {self.id()}

DESCRIPTION
    {self.doc() if self.doc is not None else "No documentation provided"}

INPUTS:
REQUIRED:
{requiredInputs}

OPTIONAL:
{optionalInputs}

OUTPUTS:
{outputs}
    """

