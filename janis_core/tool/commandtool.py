import re
from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Any, Union

from janis_core.utils.validators import Validators

from janis_core.types import (
    Selector,
    Logger,
    ParseableType,
    get_instantiated_type,
    Stdout,
    Stderr,
)

from janis_core.types.common_data_types import String, Filename
from janis_core.tool.tool import Tool, ToolTypes, TInput, TOutput
from janis_core.enums.supportedtranslations import SupportedTranslation
from janis_core.utils.metadata import ToolMetadata, Metadata


class ToolArgument:
    expr_pattern = "\$\(.*\)"

    def __init__(
        self,
        value: Any,
        prefix: Optional[str] = None,
        position: Optional[int] = 0,
        separate_value_from_prefix=None,
        doc: Optional[str] = None,
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
        self.doc = doc
        self.shell_quote = shell_quote

        if (
            self.prefix
            and self.separate_value_from_prefix is not None
            and not self.separate_value_from_prefix
            and not self.prefix.endswith("=")
        ):
            # I don't really know what this means.
            Logger.warn(
                f"Argument ({self.prefix} {self.value}) is not separating and did not end with ='"
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
        separator: str = None,
        shell_quote: bool = None,
        localise_file: bool = None,
        default: Any = None,
        doc: Optional[str] = None,
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
            doc=doc,
            shell_quote=shell_quote,
        )

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
        doc: Optional[str] = None,
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
        self.doc = doc

    def id(self):
        return self.tag


class CommandTool(Tool, ABC):
    """
    A CommandTool is an interface between Janis and a program to be executed.
    Simply put, a CommandTool has a name, a command, inputs, outputs and a container to run in.
    """

    def __init__(self, **connections):
        super().__init__(metadata_class=ToolMetadata, **connections)

    # Tool base
    @staticmethod
    @abstractmethod
    def tool() -> str:
        """
        Unique identifier of the tool
        :return:
        """
        pass

    @staticmethod
    @abstractmethod
    def base_command() -> Optional[Union[str, List[str]]]:
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

    @staticmethod
    @abstractmethod
    def container() -> str:
        """
        A link to an OCI compliant container accessible by your engine. Previously, docker().
        :return: str
        """
        pass

    @staticmethod
    @abstractmethod
    def version() -> str:
        """
        Version of the tool. Janis supports multiple versions of tools with the same ``.tool()`` value.
        The recommended format is `SemVer <https://semver.org/>`_, though you should reflect the tool version.
        :return: str
        """
        pass

    ## Other studd

    def id(self):
        return self.tool()

    @classmethod
    def __hash__(cls):
        return cls.tool()

    @classmethod
    def full_name(cls):
        if cls.version() is not None:
            return f"{cls.tool()}/{cls.version()}"
        return cls.tool()

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
        )

    def tool_inputs(self) -> List[TInput]:
        return [TInput(t.id(), t.input_type, default=t.default) for t in self.inputs()]

    def tool_outputs(self) -> List[TOutput]:
        return [TOutput(t.id(), t.output_type) for t in self.outputs()]

    def all_input_keys(self):
        return super().all_input_keys() + [
            "runtime_memory",
            "runtime_cpu",
            "runtime_disks",
        ]

    def help(self):
        import inspect

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
                f"\t\t{t.tag} ({prefix_with_space}{t.input_type.id()}{('=' + str(t.default)) if t.default is not None else ''})"
                f": {'' if t.doc is None else t.doc}"
            )

        output_format = (
            lambda t: f"\t\t{t.tag} ({t.output_type.id()}): {'' if t.doc is None else t.doc}"
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
        self, additional_inputs=None, with_resource_overrides=False, hints=None
    ):
        """
        Generate the overrides to be used with Janis. Although it may work with
        other
        :return:
        """
        d, ad = {}, additional_inputs or {}
        for i in self.inputs():
            if not i.input_type.optional or i.default or i.id() in ad:
                d[i] = ad.get(i.id(), i.default)

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

        stp = wf.step(self.tool().lower(), self.__class__(**inpmap))

        for o in self.outputs():
            wf.output(o.id(), source=stp[o.id()])

        return wf
