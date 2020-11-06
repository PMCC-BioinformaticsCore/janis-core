import os
from typing import Dict, Union
from pkg_resources import parse_version

from janis_core import ToolType, Tool, Workflow
from janis_core.utils.metadata import ToolMetadata

from janis_core.tool import test_helpers
from janis_core import CommandTool, CodeTool


class ToolEvaluator:
    @classmethod
    def evaluate(cls, tool: Tool) -> Union[str, bool]:
        """
        Evaluate a Janis tool whether they satisfy certain criteria for them to be publishable

        :param tool: Janis tool
        :type tool: Tool

        :return: error message or True if valid
        :rtype: Union[str, bool]
        """
        if tool.type() == ToolType.Workflow:
            return cls.evaluate_workflow(tool)
        elif tool.type() == ToolType.CommandTool:
            return cls.evaluate_command_tool(tool)
        elif tool.type() == ToolType.CodeTool:
            return cls.evaluate_code_tool(tool)
        raise Exception("Unrecognised tool type: " + str(tool.type()))

    @classmethod
    def evaluate_command_tool(cls, tool: CommandTool) -> Union[str, bool]:
        """
        Evaluate a Janis command line tool whether they satisfy certain criteria for them to be publishable

        :param tool: Janis command line tool
        :type tool: CommandTool

        :return: error message or True if valid
        :rtype: Union[str, bool]
        """
        evaluation = cls.evaluate_generic(tool)
        return cls._read_evaluation(evaluation)

    @classmethod
    def evaluate_code_tool(cls, tool: CodeTool) -> Union[str, bool]:
        """
        Evaluate a Janis code tool whether they satisfy certain criteria for them to be publishable

        :param tool: Janis code tool
        :type tool: CodeTool

        :return: error message or True if valid
        :rtype: Union[str, bool]
        """
        evaluation = cls.evaluate_generic(tool)
        return cls._read_evaluation(evaluation)

    @classmethod
    def evaluate_workflow(cls, wf: Workflow) -> Union[str, bool]:
        """
        Evaluate a Janis workflow whether they satisfy certain criteria for them to be publishable

        :param tool: Janis workflow
        :type tool: Workflow

        :return: error message or True if valid
        :rtype: Union[str, bool]
        """
        evaluation = cls.evaluate_generic(wf)
        return cls._read_evaluation(evaluation)

    @classmethod
    def evaluate_generic(cls, tool: Tool) -> Dict[str, str]:
        """
        Generic evaluations to be applied to all Janis tools

        :param tool: Janis tool
        :type tool: Tool

        :return: evaluation outcome keyed by evaluated category
        :rtype: Dict[str, str]
        """
        evaluation = {}

        evaluation["friendly_name"] = cls.evaluate_friendly_name(tool)
        evaluation["metadata"] = cls.evaluate_metadata(tool)
        # TODO: turn this on when we have implemented all unit tests
        # evaluation["unit_tests_exists"] = cls.evaluate_unit_test_exists(tool)
        evaluation["container"] = cls.evaluate_container(tool)
        evaluation["translation"] = cls.evaluate_translation(tool)

        return evaluation

    @staticmethod
    def evaluate_unit_test_exists(tool: Tool) -> Union[str, bool]:
        """
        Evaluate if test suite for this tool is provided

        :param tool: Janis tool
        :type tool: Tool

        :return:  error message or True if unit tests for this tool exists
        :rtype: Union[str, bool]
        """
        if tool.tests():
            return True

        return "Mising unit tests"

    @staticmethod
    def evaluate_friendly_name(tool: Tool) -> Union[str, bool]:
        """
        Evaluate if a friendly name for documentation is provided

        :param tool: Janis tool
        :type tool: Tool

        :return:  error message or True if a friendly name for this tool exists
        :rtype: Union[str, bool]
        """
        if not tool.friendly_name():
            return "Missing friendly name"

        return True

    @staticmethod
    def evaluate_metadata(tool: Tool) -> Union[str, bool]:
        """
        Evaluate if important metadata for documentation is provided

        :param tool: Janis tool
        :type tool: Tool

        :return:  error message or True if all required metadata for this tool exists
        :rtype: Union[str, bool]
        """
        METADATA_KEY_CONTRIBUTORS = "contributors"
        METADATA_KEY_CREATED_DATE = "created date"
        METADATA_KEY_INSTITUTION = "institution"

        if isinstance(tool.metadata, ToolMetadata):
            required = {
                METADATA_KEY_CONTRIBUTORS: tool.metadata.contributors,
                METADATA_KEY_CREATED_DATE: tool.metadata.dateCreated,
                METADATA_KEY_INSTITUTION: tool.metadata.institution,
            }

            missing = []
            for key, field in required.items():
                if field is None or not field:
                    missing.append(key)

            # special case, tool_provider() value overwrites contributors in the documentation
            if METADATA_KEY_INSTITUTION in missing:
                if tool.tool_provider():
                    missing.remove(METADATA_KEY_INSTITUTION)

            if missing:
                return f"Missing metadata: {', '.join(missing)}"
        # elif isinstance(self.metadata, ...):
        else:
            return "Incorrect metadata class"

        return True

    @staticmethod
    def evaluate_container(tool: Tool) -> Union[str, bool]:
        """
        Evaluate if the container specified for this tool exists in the remote registry

        :param tool: Janis tool
        :type tool: Tool

        :return:  error message or True if listed container for this tool exists in the remote registry
        :rtype: Union[str, bool]
        """
        # If there is no container, we don't need to check if the container exists in the registry
        if tool.containers() is None:
            return True

        test_helpers.verify_janis_assistant_installed()
        from janis_assistant.data.container import get_digests_from_containers

        containers = [v for k, v in tool.containers().items()]
        cache_location = os.path.join(os.getcwd(), "tests_output", "containers")
        digest = get_digests_from_containers(containers, cache_location=cache_location)

        if digest:
            return True
        else:
            return f"image {tool.container()} not found"

    @staticmethod
    def evaluate_translation(tool: Tool) -> Union[str, bool]:
        """
        Evaluate if we can successfully translate to wdl and cwl


        :param tool: Janis tool
        :type tool: Tool

        :return:  error message or True if we can successfully translate to wdl and cwl
        :rtype: Union[str, bool]
        """
        engines = test_helpers.get_available_engines()
        output_dir = os.path.join(os.getcwd(), "tests_output", tool.id())

        errors = []
        for engine in engines:
            try:
                translator = engines[engine]
                translator.translate(
                    tool,
                    export_path=output_dir,
                    should_validate=True,
                    to_console=False,
                    to_disk=True,
                )
            except Exception as e:
                errors.append(f"{translator.name}: validation failed {str(e)}")

        if errors:
            return ", ".join(errors)

        return True

    @staticmethod
    def _read_evaluation(evaluation: Dict[str, str]) -> Union[str, bool]:
        """
        Translate evaluation results into reportable format

        :param evaluation:
        :type evaluation:
        :return:
            - True if no error is found in each of the evaluation criteria
            - A string of error messages if at least one of the evaluation fails
        :rtype: Union[str, bool]
        """
        errors = []
        for field in evaluation:
            if evaluation[field] is not True:
                errors.append(evaluation[field])

        if not errors:
            return True

        return "; ".join(errors)
