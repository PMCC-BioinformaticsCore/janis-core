from typing import Tuple, Dict

from janis_core.translations import TranslatorBase


class NextflowTranslator(TranslatorBase):
    @classmethod
    def translate_workflow(
        cls,
        workflow,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> Tuple[any, Dict[str, any]]:
        pass

    @classmethod
    def translate_tool_internal(
        cls,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        pass

    @classmethod
    def translate_code_tool_internal(
        cls,
        tool,
        with_docker=True,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        pass

    @classmethod
    def unwrap_expression(cls, expression):
        pass

    @staticmethod
    def stringify_translated_workflow(wf):
        pass

    @staticmethod
    def stringify_translated_tool(tool):
        pass

    @staticmethod
    def stringify_translated_inputs(inputs):
        pass

    @staticmethod
    def workflow_filename(workflow):
        pass

    @staticmethod
    def inputs_filename(workflow):
        pass

    @staticmethod
    def tool_filename(tool):
        pass

    @staticmethod
    def resources_filename(workflow):
        pass


@staticmethod
def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
    pass
