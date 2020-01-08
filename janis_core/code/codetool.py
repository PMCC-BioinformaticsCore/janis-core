from abc import ABC, abstractmethod
from typing import List

from janis_core.tool.tool import Tool, TOutput, TInput, ToolType, ToolTypes


class CodeTool(Tool, ABC):

    # User should inherit from these blocks

    @abstractmethod
    def inputs(self) -> List[TInput]:
        pass

    @abstractmethod
    def outputs(self) -> List[TOutput]:
        pass

    @staticmethod
    @abstractmethod
    def code_block(**kwargs):
        """
        This code block must be 100% self contained. All libraries and functions must be
        imported and declared from within this block.
        :param kwargs:
        :return:
        """
        pass

    # Janis developer should inherit these methods

    @abstractmethod
    def base_command(self):
        pass

    @abstractmethod
    def script_name(self):
        pass

    @abstractmethod
    def container(self):
        pass

    @abstractmethod
    def prepared_script(self):
        pass

    # Other internal methods

    @classmethod
    def type(cls) -> ToolType:
        return ToolTypes.CodeTool

    def tool_inputs(self) -> List[TInput]:
        return self.inputs()

    def tool_outputs(self) -> List[TOutput]:
        return self.outputs()

    def generate_inputs_override(
        self, additional_inputs=None, with_resource_overrides=False, hints=None
    ):
        return {}

    def translate(
        self,
        translation: str,
        to_console=True,
        to_disk=False,
        with_docker=True,
        with_resource_overrides=False,
    ):
        from janis_core import translations

        return translations.translate_code_tool(
            self,
            translation=translation,
            to_console=to_console,
            to_disk=to_disk,
            with_docker=with_docker,
            # export_path=export_path,
        )
