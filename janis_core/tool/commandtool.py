from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Any, Union

from janis_core.tool.tool import Tool, ToolArgument, ToolInput, ToolTypes, ToolOutput
from janis_core.enums.supportedtranslations import SupportedTranslation
from janis_core.utils.metadata import ToolMetadata, Metadata


class CommandTool(Tool, ABC):
    """
    A CommandTool is an interface between Janis and a program to be executed.
    Simply put, a CommandTool has a name, a command, inputs, outputs and a container to run in.
    """

    def __init__(self):
        super().__init__()
        self._metadata = ToolMetadata()

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

    def metadata(self) -> ToolMetadata:
        return self._metadata

    @classmethod
    def type(cls):
        return ToolTypes.CommandTool

    def translate(
        self,
        translation: SupportedTranslation,
        to_console=True,
        to_disk=False,
        with_docker=True,
        with_resource_overrides=False,
    ):
        import janis_core.translations

        return janis_core.translations.translate_tool(
            self,
            translation,
            to_console=to_console,
            with_docker=with_docker,
            with_resource_overrides=with_resource_overrides,
        )

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

        metadata = self.metadata() if self.metadata() else Metadata()
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
