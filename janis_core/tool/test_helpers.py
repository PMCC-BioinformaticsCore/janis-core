import unittest
from typing import Dict, Union, List, Set
from tabulate import tabulate

import janis_core as jc


def get_all_tools(modules: List):
    shed = jc.JanisShed
    shed.hydrate(force=True, modules=modules)

    return shed.get_all_tools()


def print_test_report(failed: Dict[str, str], succeeded: Set):
    headers = ["Tool", "Status", "Description"]
    formatted_failed = [(tid, "FAILED", terror) for tid, terror in failed.items()]
    formatted_passed = [(tid, "PASSED", "") for tid in succeeded]
    print(tabulate([*formatted_failed, *formatted_passed], headers=headers))
