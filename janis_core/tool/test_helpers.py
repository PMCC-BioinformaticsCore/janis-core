# import docker
import unittest
from typing import Dict, Union, List, Set
from tabulate import tabulate

import janis_core as jc
from janis_core.translationdeps.supportedtranslations import SupportedTranslation
from janis_core.translations.cwl import CwlTranslator
from janis_core.translations.wdl import WdlTranslator
from janis_core.utils.metadata import ToolMetadata


def get_all_tools(modules: List):
    shed = jc.JanisShed
    shed.hydrate(force=True, modules=modules)

    return shed.get_all_tools()


def print_test_report(failed: Dict[str, str], succeeded: Set):
    headers = ["Tool", "Status", "Description"]
    formatted_failed = [(tid, "FAILED", terror) for tid, terror in failed.items()]
    formatted_passed = [(tid, "PASSED", "") for tid in succeeded]
    print(tabulate([*formatted_failed, *formatted_passed], headers=headers))


def evaluate(tool) -> Dict[str, str]:
    evaluation = {}

    evaluation["friendly_name"] = evaluate_friendly_name(tool)
    evaluation["metadata"] = evaluate_metadata(tool)
    evaluation["unit_tests_exists"] = evaluate_unit_test_exists(tool)
    # evaluation['container'] = evaluate_container(tool)
    evaluation["translation"] = evaluate_translation(tool)

    return evaluation


def run_tool_unit_tests(tool: jc.Tool):
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tool.UnitTestClass)
    return unittest.TextTestRunner().run(suite)


def evaluate_unit_test_exists(tool: jc.Tool) -> Union[str, bool]:
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(tool.UnitTestClass)
    if suite.countTestCases() > 0:
        return True

    return "Mising unit tests"


def evaluate_friendly_name(tool: jc.Tool) -> Union[str, bool]:
    if tool.friendly_name() is None:
        return "Missing friendly name"

    return True


def evaluate_metadata(tool: jc.Tool) -> Union[str, bool]:
    if isinstance(tool.metadata, ToolMetadata):
        required = {
            "contributors": tool.metadata.contributors,
            "created date": tool.metadata.dateCreated,
            "institution": tool.metadata.institution,
        }

        missing = []
        for key, field in required.items():
            if field is None or not field:
                missing.append(key)

        if missing:
            return f"Missing metadata: {', '.join(missing)}"
    # elif isinstance(self.metadata, ...):
    else:
        return "Incorrect metadata class"

    return True


# def evaluate_container(tool: jc.Tool) -> Union[str, bool]:
#     """
#     Evaluate if the image specified for this tool exists in the remote registry
#     """
#     client = docker.from_env()
#     try:
#         client.images.get_registry_data(tool.container())
#     except docker.errors.NotFound as e:
#         return f"image {tool.container()} not found"
#     except Exception as e:
#         return f"image {tool.container()}: {str(e)}"
#
#     return True


def evaluate_translation(tool: jc.Tool):
    cwl_file_path = f"/tmp/janis/tests/{tool.id()}/cwl"
    wdl_file_path = f"/tmp/janis/tests/{tool.id()}/wdl"

    tool.translate(
        SupportedTranslation.CWL,
        to_console=False,
        to_disk=True,
        export_path=cwl_file_path,
    )
    tool.translate(
        SupportedTranslation.WDL,
        to_console=False,
        to_disk=True,
        export_path=wdl_file_path,
    )

    # TODO: translate and validate
    CwlTranslator.validate_command_for(cwl_file_path, "", "", "")
    WdlTranslator.validate_command_for(wdl_file_path, "", "", "")

    return True
