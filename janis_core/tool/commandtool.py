import re
from abc import ABC, abstractmethod
from inspect import isclass
from typing import List, Dict, Optional, Any, Union, Callable, Set

from janis_core.tool.documentation import (
    InputDocumentation,
    DocumentationMeta,
    OutputDocumentation,
    InputQualityType,
)
from janis_core.utils.validators import Validators

from janis_core.types import ParseableType, get_instantiated_type, Stdout, Stderr

from janis_core.types.common_data_types import String, Filename
from janis_core.tool.tool import Tool, ToolTypes, TInput, TOutput
from janis_core.translationdeps.supportedtranslations import SupportedTranslation
from janis_core.utils.logger import Logger
from janis_core.operators.selectors import Selector
from janis_core.utils.metadata import ToolMetadata, Metadata


class ToolArgument:
    expr_pattern = "\$\(.*\)"

    def __init__(
        self,
        value: Any,
        prefix: Optional[str] = None,
        position: Optional[int] = 0,
        separate_value_from_prefix=None,
        doc: Optional[Union[str, DocumentationMeta]] = None,
        shell_quote: bool = None,
    ):
        """
        A ``ToolArgument`` is a CLI parameter that cannot be override (at runtime).
        The value can


        :param value:
        :type value: ``str`` | ``janis.InputSelector`` | ``janis.StringFormatter``
        :param position: The position of the input to be applied. (Default = 0, after the base_command).
        :param prefix: The prefix to be appended before the element. (By default, a space will also be applied, see ``separate_value_from_prefix`` for more information)
        :param separate_value_from_prefix: (Default: True) Add a space between the prefix and value when ``True``.
        :param doc: Documentation string for the argument, this is used to generate the tool documentation and provide
        :param shell_quote: Stops shell quotes from being applied in all circumstances, useful when joining multiple commands together.
        """

        self.prefix: Optional[str] = prefix
        self.value = value
        self.position: Optional[int] = position
        self.is_expression = (
            isinstance(self.value, Selector)
            or (re.match(self.expr_pattern, self.value) is not None)
            if self.value
            else None
        )
        self.separate_value_from_prefix = separate_value_from_prefix
        self.doc: DocumentationMeta = doc if isinstance(
            doc, InputDocumentation
        ) else InputDocumentation(doc)
        self.shell_quote = shell_quote

        if (
            self.prefix
            and self.separate_value_from_prefix is not None
            and not self.separate_value_from_prefix
            and not self.prefix.endswith("=")
        ):
            # I don't really know what this means.
            Logger.warn(
                f"Argument ({self.prefix}{self.value}) is not separating and did not end with ='"
            )


# This should really be a CommandToolInput
class ToolInput(ToolArgument):
    def __init__(
        self,
        tag: str,
        input_type: ParseableType,
        position: Optional[int] = None,
        prefix: Optional[str] = None,
        separate_value_from_prefix: bool = None,
        prefix_applies_to_all_elements: bool = None,
        presents_as: str = None,
        secondaries_present_as: Dict[str, str] = None,
        separator: str = None,
        shell_quote: bool = None,
        localise_file: bool = None,
        default: Any = None,
        doc: Optional[Union[str, InputDocumentation]] = None,
    ):
        """
        A ``ToolInput`` represents an input to a tool, with parameters that allow it to be bound on the command line.
        The ToolInput must have either a position or prefix set to be bound onto the command line.

        :param tag: The identifier of the input (unique to inputs and outputs of a tool)
        :param input_type: The data type that this input accepts
        :type input_type: ``janis.ParseableType``
        :param position: The position of the input to be applied. (Default = 0, after the base_command).
        :param prefix: The prefix to be appended before the element. (By default, a space will also be applied, see ``separate_value_from_prefix`` for more information)
        :param separate_value_from_prefix: (Default: True) Add a space between the prefix and value when ``True``.
        :param prefix_applies_to_all_elements: Applies the prefix to each element of the array (Array inputs only)
        :param shell_quote: Stops shell quotes from being applied in all circumstances, useful when joining multiple commands together.
        :param separator: The separator between each element of an array (defaults to ' ')
        :param localise_file: Ensures that the file(s) are localised into the execution directory.
        :param default: The default value to be applied if the input is not defined.
        :param doc: Documentation string for the ToolInput, this is used to generate the tool documentation and provide
        hints to the user.
        """
        super().__init__(
            value=None,
            prefix=prefix,
            position=position,
            separate_value_from_prefix=separate_value_from_prefix,
            doc=None,
            shell_quote=shell_quote,
        )

        self.doc: InputDocumentation = doc if isinstance(
            doc, DocumentationMeta
        ) else InputDocumentation(doc=doc)

        # if default is not None:
        #     input_type.optional = True

        if not Validators.validate_identifier(tag):
            raise Exception(
                f"The identifier '{tag}' was not validated because {Validators.reason_for_failure(tag)}"
            )

        self.tag: str = tag
        self.input_type: ParseableType = get_instantiated_type(input_type)
        self.default = default
        self.prefix_applies_to_all_elements = prefix_applies_to_all_elements
        self.separator = separator

        self.localise_file = localise_file
        self.presents_as = presents_as
        self.secondaries_present_as = secondaries_present_as
        # if isinstance(input_type, Array):
        #     if self.prefix_applies_to_all_elements is None and self.separator is None:
        # self.separator = " "

    def id(self):
        return self.tag


# This should really be a CommandToolOutput
class ToolOutput:
    def __init__(
        self,
        tag: str,
        output_type: ParseableType,
        glob: Optional[Union[Selector, str]] = None,
        presents_as: str = None,
        secondaries_present_as: Dict[str, str] = None,
        doc: Optional[Union[str, OutputDocumentation]] = None,
    ):
        """
        A ToolOutput instructs the the engine how to collect an output and how
        it may be referenced in a workflow.

        :param tag: The identifier of a output, must be unique in the inputs and outputs.
        :param output_type: The type of output that is being collected.
        :param glob: How to collect this output, can accept any :class:`janis.Selector`.
        :param doc: Documentation on what the output is, used to generate docs.
        """

        if not Validators.validate_identifier(tag):
            raise Exception(
                f"The identifier '{tag}' was invalid because {Validators.reason_for_failure(tag)}"
            )

        self.tag = tag
        self.output_type: ParseableType = get_instantiated_type(output_type)

        if not glob and not (
            isinstance(self.output_type, Stdout) or isinstance(self.output_type, Stderr)
        ):
            raise Exception(
                "ToolOutput expects a glob when the output type is not Stdout / Stderr"
            )

        self.glob = glob
        self.presents_as = presents_as
        self.secondaries_present_as = secondaries_present_as
        self.doc = (
            doc
            if isinstance(doc, OutputDocumentation)
            else OutputDocumentation(doc=doc)
        )

    def id(self):
        return self.tag


class CommandTool(Tool, ABC):
    """
    A CommandTool is an interface between Janis and a program to be executed.
    Simply put, a CommandTool has a name, a command, inputs, outputs and a container to run in.

    This class can be inherited to created a CommandTool, else a CommandToolBuilder may be used.
    """

    def __init__(self, **connections):
        super().__init__(metadata_class=ToolMetadata, **connections)

    # Tool base
    @abstractmethod
    def tool(self) -> str:
        """
        Unique identifier of the tool
        :return:
        """
        pass

    @abstractmethod
    def base_command(self) -> Optional[Union[str, List[str]]]:
        """
        The command of the tool to execute, usually the tool name or path and not related to any inputs.
        This field will always come before any inputs or arguments, though it's possible to omit this
        field and the program will use the first ordered argument / position.
        :return: Optional[Union[str, List[str]]]
        """
        pass

    @abstractmethod
    def inputs(self) -> List[ToolInput]:
        """
        A list of named tool inputs that will be used to create the command line. See :class:`janis.ToolInput`
        for options on how to configure this command line binding.
        :return: List[janis.ToolInput]
        """
        pass

    def arguments(self) -> Optional[List[ToolArgument]]:
        """
        A list of arguments that will be used to create the command line. Although they are not directly
        addressable as inputs, it's possible to use use a :class:`janis.InputSelector` or
        :class:`janis.StringFormatter` in the value field. See :class:`janis.ToolArgument` for
        options on how to configure a this command line binding.
        :return: List[janis.ToolArgument]
        """
        return None

    @abstractmethod
    def outputs(self) -> List[ToolOutput]:
        """
        A list of named outputs of the tool. Each :class:`janis.ToolOutput` has a ``glob`` field that
        can be used to select the outputs, see its documentation for more information.
        :return:
        """
        pass

    def env_vars(self) -> Optional[Dict[str, Union[str, Selector]]]:
        return None

    # Tool versions

    @abstractmethod
    def container(self) -> str:
        """
        A link to an OCI compliant container accessible by your engine. Previously, docker().
        :return: str
        """
        pass

    @abstractmethod
    def version(self) -> str:
        """
        Version of the tool. Janis supports multiple versions of tools with the same ``.tool()`` value.
        The recommended format is `SemVer <https://semver.org/>`_, though you should reflect the tool version.
        :return: str
        """
        pass

    ## Other studd

    def id(self):
        return self.tool()

    def __hash__(self):
        return hash(self.tool())

    def full_name(self):
        if self.version() is not None:
            return f"{self.tool()}/{self.version()}"
        return self.tool()

    def memory(self, hints: Dict[str, Any]) -> Optional[float]:
        """
        These values are used to generate a separate runtime.json / runtime.yaml input
        that can be passed to the execution engine to fill in for the specified hints.

        These are now (2019-04-10) to be kept out of the workflow, to leave the workflow
        truly portable.

        This memory must be in GB!
        :param hints: Dict[Key: value] of hints
        :return: Optional[int]
        """
        return None

    def cpus(self, hints: Dict[str, Any]) -> Optional[int]:
        """
        These values are used to generate a separate runtime.json / runtime.yaml input
        that can be passed to the execution engine to fill in for the specified hints.

        These are now (2019-04-10) to be kept out of the workflow, to leave the workflow
        truly portable.

        The CPU must be a whole number. If your tool contains threads
        :return:
        """
        return None

    @classmethod
    def type(cls):
        return ToolTypes.CommandTool

    def translate(
        self,
        translation: SupportedTranslation,
        to_console=True,
        to_disk=False,
        export_path=None,
        with_docker=True,
        with_resource_overrides=False,
        allow_empty_container=False,
    ):
        import janis_core.translations

        return janis_core.translations.translate_tool(
            self,
            translation,
            to_console=to_console,
            to_disk=to_disk,
            export_path=export_path,
            with_docker=with_docker,
            with_resource_overrides=with_resource_overrides,
            allow_empty_container=allow_empty_container,
        )

    def tool_inputs(self) -> List[TInput]:
        return [
            TInput(t.id(), t.input_type, default=t.default, doc=t.doc)
            for t in self.inputs()
        ]

    def tool_outputs(self) -> List[TOutput]:
        return [TOutput(t.id(), t.output_type, doc=t.doc) for t in self.outputs()]

    def all_input_keys(self):
        return super().all_input_keys() + [
            "runtime_memory",
            "runtime_cpu",
            "runtime_disks",
        ]

    def help(self):
        import inspect

        tb = " " * 4
        path = inspect.getfile(self.__class__)

        ins = sorted(
            self.inputs(), key=lambda i: i.position if i.position is not None else 0
        )
        args = ""
        if self.arguments():
            args = " " + " ".join(
                f"{(a.prefix if a.prefix is not None else '') + ' ' if (a.prefix is not None and a.separate_value_from_prefix) else ''}{a.value}"
                for a in self.arguments()
            )

        prefixes = " -" + "".join(
            i.prefix.replace("-", "").replace(" ", "")
            for i in ins
            if i.prefix is not None
        )

        metadata = self.metadata
        docker = self.container()

        base = (
            (
                self.base_command()
                if isinstance(self.base_command(), str)
                else " ".join(self.base_command())
            )
            if self.base_command()
            else ""
        )
        command = base + args + prefixes

        def input_format(t: ToolInput):
            prefix_with_space = ""
            if t.prefix is not None:
                prefix_with_space = (
                    (t.prefix + ": ")
                    if (t.separate_value_from_prefix is not False)
                    else t.prefix
                )
            return (
                f"{2 * tb}{t.tag} ({prefix_with_space}{t.input_type.id()}{('=' + str(t.default)) if t.default is not None else ''})"
                f": {'' if t.doc is None else t.doc}"
            )

        output_format = (
            lambda t: f"{2 * tb}{t.tag} ({t.output_type.id()}): {'' if t.doc is None else t.doc}"
        )

        requiredInputs = "\n".join(
            input_format(x) for x in ins if not x.input_type.optional
        )
        optionalInputs = "\n".join(
            input_format(x) for x in ins if x.input_type.optional
        )
        outputs = "\n".join(output_format(o) for o in self.outputs())

        return f"""
    Pipeline tool: {path} ({self.id()})
NAME
    {self.id()}

SYNOPSIS
    {command}

DOCKER
    {docker}

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

    def generate_inputs_override(
        self,
        additional_inputs=None,
        with_resource_overrides=False,
        hints=None,
        include_defaults=True,
        values_to_ignore: Set[str] = None,
        quality_type: List[InputQualityType] = None,
    ):
        """
        Generate the overrides to be used with Janis. Although it may work with
        other
        :return:
        """
        d, ad = {}, additional_inputs or {}
        for i in self.inputs():
            if (
                (
                    not i.input_type.optional
                    or i.id() in ad
                    or (include_defaults and i.default)
                )
                and not (values_to_ignore and i.id() in values_to_ignore)
                and (not (i.doc and quality_type) or i.doc.quality in quality_type)
            ):
                d[i.id()] = ad.get(i.id(), i.default)

        if with_resource_overrides:
            cpus = self.cpus(hints) or 1
            mem = self.memory(hints)
            d.update(
                {
                    "runtime_memory": mem,
                    "runtime_cpu": cpus,
                    "runtime_disks": "local-disk 60 SSD",
                }
            )

        return d

    def wrapped_in_wf(self):
        from copy import copy
        from janis_core.workflow.workflow import WorkflowBuilder

        wf = WorkflowBuilder(self.id() + "Wf")
        inpmap = {}
        for i in self.inputs():

            if isinstance(i.input_type, Filename):
                intp = String(optional=True)
            else:
                intp = copy(i.input_type)
                if i.default:
                    intp.optional = True

            inpmap[i.id()] = wf.input(i.id(), intp)

        stp = wf.step(self.tool().lower(), self(**inpmap))

        for o in self.outputs():
            wf.output(o.id(), source=stp[o.id()])

        return wf


class CommandToolBuilder(CommandTool):
    def tool(self) -> str:
        return self._tool

    def friendly_name(self):
        return self._friendly_name

    def base_command(self) -> Optional[Union[str, List[str]]]:
        return self._base_command

    def inputs(self) -> List[ToolInput]:
        return self._inputs

    def arguments(self):
        return self._arguments

    def outputs(self) -> List[ToolOutput]:
        return self._outputs

    def container(self) -> str:
        return self._container

    def version(self) -> str:
        return self._version

    def tool_provider(self):
        return self._tool_provider

    def tool_module(self):
        return self._tool_module

    def env_vars(self):
        return self._env_vars

    def cpus(self, hints: Dict[str, Any]):
        if self._cpu is None:
            return None
        if isinstance(self._cpu, int) or isinstance(self._cpu, float):
            return self._cpu

        if callable(self._cpu):
            return self._cpu(hints)

        raise Exception(
            f"Janis does not recognise {type(self._cpu)} as a valid CPU type"
        )

    def memory(self, hints: Dict[str, Any]):
        if self._memory is None:
            return None
        if isinstance(self._memory, int) or isinstance(self._memory, float):
            return self._memory

        if callable(self._memory):
            return self._memory(hints)

        raise Exception(
            f"Janis does not recognise {type(self._memory)} as a valid memory type"
        )

    def __init__(
        self,
        tool: str,
        base_command: Optional[Union[str, List[str]]],
        inputs: List[ToolInput],
        outputs: List[ToolOutput],
        container: str,
        version: str,
        friendly_name: Optional[str] = None,
        arguments: List[ToolArgument] = None,
        env_vars: Dict = None,
        tool_module: str = None,
        tool_provider: str = None,
        metadata: ToolMetadata = None,
        cpu: Union[int, Callable[[Dict[str, Any]], int]] = None,
        memory: Union[int, Callable[[Dict[str, Any]], int]] = None,
    ):
        """
        Builder for a CommandTool.

        :param tool: Unique identifier of the tool
        :param friendly_name: A user friendly name of your tool (must be implemented for generated docs)
        :param base_command: The command of the tool to execute, usually the tool name or path and not related to any inputs.
        :param inputs: A list of named tool inputs that will be used to create the command line.
        :param outputs: A list of named outputs of the tool; a ``ToolOutput`` declares how the output is captured.
        :param arguments: A list of arguments that will be used to create the command line.
        :param container: A link to an OCI compliant container accessible by the engine.
        :param version: Version of the tool.
        :param env_vars: A dictionary of environment variables that should be defined within the container.
        :param tool_module: Unix, bioinformatics, etc.
        :param tool_provider: The manafacturer of the tool, eg: Illumina, Samtools
        :param metadata: Metadata object describing the Janis tool interface
        :param cpu: An integer, or function that takes a dictionary of hints and returns an integer
        :param memory: An integer, or function that takes a dictionary of hints and returns an integer
        """

        super().__init__()

        self._tool = tool
        self._friendly_name = friendly_name
        self._base_command = base_command
        self._inputs = inputs
        self._outputs = outputs
        self._container = container
        self._version = version
        self._arguments = arguments
        self._env_vars = env_vars
        self._tool_module = tool_module
        self._tool_provider = tool_provider
        self._metadata = metadata
        self._cpu = cpu
        self._memory = memory
