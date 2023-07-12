import operator
import sys
import os
from abc import ABC, abstractmethod
from enum import Enum
from typing import Optional, List, Dict, Set, Any
from uuid import uuid4

from janis_core.types.common_data_types import Array
from janis_core.tool.documentation import (
    InputDocumentation,
    OutputDocumentation,
    InputQualityType,
)
from janis_core.types import get_instantiated_type, DataType
from janis_core.utils import find_duplicates
from janis_core.utils.metadata import Metadata
from janis_core.utils.validators import Validators
from janis_core.tool.test_classes import (
    TTestCase,
    TTestExpectedOutput,
    TTestPreprocessor,
)
from nose.tools import nottest

class ToolType(Enum):
    Workflow = "workflow"
    CommandTool = "command-tool"
    CodeTool = "code-tool"

    def __str__(self):
        if self == ToolType.Workflow:
            return "Workflow"
        elif self == ToolType.CommandTool:
            return "CommandTool"
        elif self == ToolType.CodeTool:
            return "CodeTool"
        return "".join(a.title() for a in self.value.split("-"))


class TInput(object):
    def __init__(
        self, tag: str, intype: DataType, default=None, doc: InputDocumentation = None
    ):
        self.tag = tag
        self.intype = get_instantiated_type(intype)
        self.default = default
        self.doc = doc

    def __repr__(self):
        items = ["{self.id()}", self.intype.id()]
        if self.default is not None:
            items.append("default=" + str(self.default))
        return f"TInput({', '.join(items)})"

    def id(self) -> str:
        return self.tag


class TOutput(object):
    def __init__(self, tag, outtype, doc: OutputDocumentation = None):
        self.tag = tag
        self.outtype = get_instantiated_type(outtype)
        self.doc: Optional[OutputDocumentation] = doc

    def __repr__(self):
        return f'TOutput("{self.id()}", {self.outtype.id()})'

    def id(self) -> str:
        return self.tag


class Tool(ABC, object):
    """
    One of Workflow, CommandLineTool, ExpressionTool* (* unimplemented)
    """

    TEST_DATA_FOLDER = "test_data"

    def __init__(self, metadata_class=Metadata, **connections):
        """
        :param metadata_class:
        :param connections:
        """
        self.uuid: str = str(uuid4())
        self.metadata: metadata_class = metadata_class()
        meta = self.bind_metadata()
        if meta:
            self.metadata = meta

        self.connections: dict[str, Any] = connections

    def __repr__(self):
        return f"{str(self.type())}<{self.id()}>"

    @classmethod
    @abstractmethod
    def type(cls) -> ToolType:
        raise Exception(f"'{cls}' must implement type() method")

    @abstractmethod
    def containers(self) -> Dict[str, str]:
        pass

    @abstractmethod
    def id(self) -> str:
        raise Exception("Must implement id() method")

    def versioned_id(self) -> str:
        if self.version() is not None:
            return Validators.transform_identifier_to_be_valid(
                f"{self.id()}/{self.version()}", "_"
            )
        return self.id()

    def tool_module(self):
        return None

    def tool_provider(self):
        return None

    @abstractmethod
    def tool_inputs(self) -> List[TInput]:
        raise Exception("Must implement inputs() method")

    @abstractmethod
    def tool_outputs(self) -> List[TOutput]:
        raise Exception("Must implement outputs() method")

    def inputs_map(self) -> Dict[str, TInput]:
        ins = self.tool_inputs()
        indict = {inp.tag: inp for inp in ins}

        if len(ins) != len(indict):
            dups = find_duplicates([i.tag for i in ins])
            dupstext = ", ".join(dups)
            raise Exception(
                f"There are {len(dups)} duplicate values in {self.id()}'s inputs: {dupstext}"
            )

        return indict

    def outputs_map(self) -> Dict[str, TOutput]:
        outs = self.tool_outputs()
        outdict = {outp.tag: outp for outp in outs}

        if len(outs) != len(outdict):
            dups = find_duplicates([o.tag for o in outs])
            dupstext = ", ".join(dups)
            raise Exception(
                f"There are {len(dups)} duplicate values in {self.id()}'s outputs: {dupstext}"
            )

        return outdict

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
    def has_tool_with_no_container(self):
        pass

    @abstractmethod
    def generate_inputs_override(
        self,
        additional_inputs=None,
        with_resource_overrides=False,
        hints=None,
        include_defaults=True,
        values_to_ignore: Set[str] = None,
        quality_type: List[InputQualityType] = None,
    ):
        pass

    def __call__(self, **connections):
        self.connections = connections
        return self

    @abstractmethod
    def version(self):
        return None

    def doc(self) -> Optional[str]:
        return None

    def translate(
        self,
        dest_fmt: str,
        mode: Optional[str] = None,

        # file io
        to_disk: Optional[bool] = None,
        export_path: Optional[str] = None,
        should_zip: Optional[bool] = None,   
        to_console: Optional[bool] = None,
        tool_to_console: Optional[bool] = None,
        write_inputs_file: Optional[bool] = None,
        
        # inputs
        additional_inputs: Optional[dict[str, str]] = None,
        hints: Optional[dict[str, str]] = None,
        
        # containers
        with_container: Optional[bool] = None,
        allow_empty_container: Optional[bool] = None,
        container_override: Optional[str | dict[str, Any]] = None,
        
        # resouces
        with_resource_overrides: Optional[bool] = None,
        merge_resources: Optional[bool] = None,
        max_cores: Optional[int] = None,
        max_mem: Optional[int] = None,
        max_duration: Optional[int] = None,
        
        # misc
        render_comments: Optional[bool] = None,
        should_validate: Optional[bool] = None,
        ) -> Any:
        from janis_core.translations import translate
        return translate(
            entity=self,
            dest_fmt=dest_fmt,
            mode=mode,

            # file io
            to_disk=to_disk,
            export_path=export_path,
            should_zip=should_zip,
            to_console=to_console,
            tool_to_console=tool_to_console,
            write_inputs_file=write_inputs_file,
            
            # inputs
            additional_inputs=additional_inputs,
            hints=hints,
            
            # containers
            with_container=with_container,
            allow_empty_container=allow_empty_container,
            container_override=container_override,
            
            # resouces
            with_resource_overrides=with_resource_overrides,
            merge_resources=merge_resources,
            max_cores=max_cores,
            max_mem=max_mem,
            max_duration=max_duration,
            
            # misc
            render_comments=render_comments,
            should_validate=should_validate,
        )

    def bind_metadata(self):
        """
        A convenient place to add metadata about the tool. You are guaranteed that self.metadata will exist.
        It's possible to return a new instance of the ToolMetadata / WorkflowMetadata which will be rebound.
        This is usually called after the initialiser, though it may be called multiple times.
        :return:
        """
        return self.metadata

    def help(self):
        import inspect

        tb = " " * 4

        path = inspect.getfile(self.__class__)

        ins = self.tool_inputs()

        metadata = self.metadata

        def input_format(t: TInput):
            return (
                f"{2 * tb}{t.id()} ({t.intype.id()}{('=' + str(t.default)) if t.default is not None else ''})"
                f": {'' if t.doc is None else t.doc}"
            )

        output_format = (
            lambda t: f"{2 * tb}{t.id()} ({t.outtype.id()}): {'' if t.doc is None else t.doc}"
        )

        requiredInputs = "\n".join(
            input_format(x) for x in ins if not x.intype.optional
        )
        optionalInputs = "\n".join(input_format(x) for x in ins if x.intype.optional)
        outputs = "\n".join(output_format(o) for o in self.tool_outputs())

        return f"""
Pipeline tool: {path} ({self.id()})
NAME
    {self.id()} ({self.friendly_name()})

DOCUMENTATION URL
    {metadata.documentationUrl if metadata.documentationUrl else "No url provided"}

DESCRIPTION
    {metadata.documentation if metadata.documentation else "No documentation provided"}

INPUTS:
    REQUIRED:
{requiredInputs}

    OPTIONAL:
{optionalInputs}

OUTPUTS:
{outputs}
"""

    @nottest
    def tests(self) -> Optional[List[TTestCase]]:
        """
        A list of test cases for this tool
        """
        return None

    def minimal_test(self) -> List[TTestExpectedOutput]:
        """
        A minimal test simply checks if output files exist (if their sizes are bigger than 0).
        It should be used when we don't know what outputs we can expect from a tool. It should
        be called within a TTestCase. Be aware that we still need to know the inputs.

        :return: List of expected outputs
        :rtype: List[TTestExpectedOutput]
        """
        outcome = []
        for i in self.tool_outputs():
            preprocessor = (
                TTestPreprocessor.ListOfFilesExist
                if i.outtype.is_base_type(Array)
                else TTestPreprocessor.FileSize
            )
            comparison = operator.eq if i.outtype.is_base_type(Array) else operator.gt
            expected_value = True if i.outtype.is_base_type(Array) else 0
            secondary_files_suffixes = (
                i.outtype.fundamental_type().secondary_files()
                if i.outtype.is_base_type(Array)
                else i.outtype.secondary_files()
            )
            outcome += [
                TTestExpectedOutput(
                    tag=i.tag,
                    preprocessor=preprocessor,
                    operator=comparison,
                    expected_value=expected_value,
                )
            ]
            if secondary_files_suffixes is not None:
                for suffix in secondary_files_suffixes:
                    outcome += [
                        TTestExpectedOutput(
                            tag=i.tag,
                            suffix_secondary_file=suffix,
                            preprocessor=preprocessor,
                            operator=comparison,
                            expected_value=expected_value,
                        )
                    ]
        return outcome

    @classmethod
    @nottest
    def test_data_path(cls):
        module_path = os.path.dirname(sys.modules[cls.__module__].__file__)
        return os.path.join(module_path, cls.TEST_DATA_FOLDER)

    @classmethod
    @nottest
    def skip_test(cls) -> bool:
        """
        Sometimes, we may want to skip tests for some tools because they are not ready yet
        """
        return False
