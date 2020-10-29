from typing import Dict, Union, List, Set
from tabulate import tabulate

import janis_core as jc

janis_assistant_version_required_min = "0.10.6"


def get_all_tools(modules: List):
    shed = jc.JanisShed
    shed.hydrate(force=True, modules=modules)

    return shed.get_all_tools()


def get_one_tool(tool_id: str, modules: List, version: str = None):
    shed = jc.JanisShed
    shed.hydrate(force=True, modules=modules)

    return shed.get_tool(tool=tool_id, version=version)


def print_test_report(failed: Dict[str, str], succeeded: Set, first_column_header: str = "Tool"):
    headers = [first_column_header, "Status", "Description"]
    formatted_failed = [(tid, "FAILED", terror) for tid, terror in failed.items()]
    formatted_passed = [(tid, "PASSED", "") for tid in succeeded]
    print(tabulate([*formatted_failed, *formatted_passed], headers=headers))

