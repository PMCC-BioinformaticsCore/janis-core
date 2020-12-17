import hashlib
from typing import Dict, Union, List, Set, Optional
from datetime import datetime
from tabulate import tabulate
from pkg_resources import parse_version
from nose.tools import nottest

import janis_core as jc
from janis_core.translations.translationbase import TranslatorBase

janis_assistant_version_required_min = "0.11.0"


def verify_janis_assistant_installed():
    """
    Check if the correct version of janis assistant is installed
    """
    min_version_required = janis_assistant_version_required_min

    try:
        import janis_assistant

        if parse_version(janis_assistant.__version__) < parse_version(
            min_version_required
        ):
            raise Exception(
                f"to run this test, janis_asisstant >= {min_version_required}"
                f" must be installed. Installed version is {janis_assistant.__version__}"
            )

    except Exception as e:
        raise e


def get_available_engines() -> Dict[str, TranslatorBase]:
    """
    Get a list of available engines to run the test suite against
    """
    verify_janis_assistant_installed()
    from janis_assistant.engines.enginetypes import EngineType
    from janis_assistant.engines import (
        get_ideal_specification_for_engine,
        get_engine_type,
    )

    engines = {}
    for engine in EngineType.engines():
        engine_type = get_engine_type(engine)
        supported_translation = get_ideal_specification_for_engine(engine_type)
        engines[engine] = supported_translation.get_translator()

    return engines


def get_all_tools(modules: List) -> List[List[jc.Tool]]:
    """
    Get all available tools in the list of modules

    :param modules: a list of python modules to find Janis tools from
    :type modules: List

    :return: A list of a list of versioned tools
    :rtype: List[List[jc.Tool]]
    """
    shed = jc.JanisShed
    shed.hydrate(force=True, modules=modules)

    return shed.get_all_tools()


def get_one_tool(
    tool_id: str, modules: Optional[List] = None, version: Optional[str] = None
) -> jc.Tool:
    """
    Get one tool given id and the modules where to search it from

    :param tool_id: Janis tool id
    :type tool_id: str
    :param modules: a list of python modules to find Janis tool from
    :type modules: List
    :param version: tool version
    :type version:  Optional[str]

    :return: Janis tool
    :rtype: Tool
    """
    shed = jc.JanisShed
    shed.hydrate(force=True, modules=modules)

    return shed.get_tool(tool=tool_id, version=version)


@nottest
def print_test_report(
    failed: Dict[str, str],
    succeeded: Set,
    skipped: Optional[Set] = set(),
    id_column_header: str = "Tool",
):
    """
    Print test report in tabular text format

    :param failed: a dictionary of error messages keyed by tool id or other unique identifier
    :type failed: Dict[str, str]
    :param succeeded: a unique set of tool ids or other unique identifiers that have passed tests
    :type succeeded: Set
    :param id_column_header: if not provided, we will call it "Tool"
    :type id_column_header: str

    """
    headers = [id_column_header, "Status", "Description"]
    formatted_failed = [(tid, "FAILED", terror) for tid, terror in failed.items()]
    formatted_skipped = [(tid, "SKIPPED", "") for tid in skipped]
    formatted_passed = [(tid, "PASSED", "") for tid in succeeded]
    print(
        tabulate(
            [*formatted_failed, *formatted_skipped, *formatted_passed], headers=headers
        )
    )


@nottest
def hash_filename(source: str, last_modified: str) -> str:
    """

    :return:
    :rtype:
    """
    # If we cannot find lat modified date, we will always download again
    components = source + (last_modified or str(datetime.now()))
    hash_md5 = hashlib.md5(str.encode(components))
    hashed_filename = hash_md5.hexdigest()

    return hashed_filename
