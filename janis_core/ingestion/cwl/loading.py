
import ruamel.yaml
from typing import Any, Optional
from janis_core import settings
from janis_core.messages import log_warning
from .preprocessing import convert_cwl_types_to_python
from .preprocessing import handle_inline_cltool_identifiers


DEFAULT_PARSER_VERSION = "v1.2"


def load_cwl_version(doc: str) -> str:
    """loads a cwl document & returns the version field"""
    
    # adjust file path
    if doc.startswith("file://"):
        doc = doc[6:]

    # load tool into memory
    with open(doc) as fp:
        tool_dict = ruamel.yaml.load(fp, Loader=ruamel.yaml.Loader) # type: ignore
    
    if "cwlVersion" not in tool_dict:
        if settings.ingest.cwl.REQUIRE_CWL_VERSION: 
            raise Exception(f"Couldn't find cwlVersion in tool {doc}")
        else:
            msg = f'no cwl version was specified in {doc}. fell back to cwl v1.2 for ingestion.'
            log_warning(uuid=None, msg=msg)
            return DEFAULT_PARSER_VERSION
    
    # return version
    return tool_dict["cwlVersion"]


def load_cwl_document(doc: str, version: Optional[str]=None) -> Any:
    """loads a cwl document & returns the in-memory cwlutils object"""
    
    if not version:
        version = load_cwl_version(doc)

    cwl_utils = load_cwl_utils_from_version(version)
    loaded_doc = cwl_utils.load_document(doc)  # type: ignore

    # convert yaml datatypes to python datatypes
    loaded_doc = convert_cwl_types_to_python(loaded_doc, cwl_utils)
    
    # convert random ids (occurs for inline clt definition) to meaningful ids
    if isinstance(loaded_doc, cwl_utils.Workflow):
        loaded_doc = handle_inline_cltool_identifiers(loaded_doc, cwl_utils)

    return loaded_doc



def convert_etool_to_cltool(etool: Any, version: str) -> Any:
    """uses cwlutils etool_to_cltool() to convert an ExpressionTool into a CommandLineTool"""
    
    etool_to_cltool = load_etool_to_cltool_from_version(version)
    cltool = etool_to_cltool(etool)
    return cltool

def load_cwl_utils_from_version(version: str) -> Any:
    if version == "v1.0":
        import cwl_utils.parser.cwl_v1_0 as cwlutils
    elif version == "v1.1":
        import cwl_utils.parser.cwl_v1_1 as cwlutils
    elif version == "v1.2":
        import cwl_utils.parser.cwl_v1_2 as cwlutils
    else:
        print(
            f"Didn't recognise CWL version {version}, loading default: {DEFAULT_PARSER_VERSION}"
        )
        cwlutils = load_cwl_utils_from_version(DEFAULT_PARSER_VERSION)
    return cwlutils

def load_etool_to_cltool_from_version(version: str) -> Any:
    if version == "v1.0":
        from cwl_utils.cwl_v1_0_expression_refactor import etool_to_cltool
    elif version == "v1.1":
        from cwl_utils.cwl_v1_0_expression_refactor import etool_to_cltool
    elif version == "v1.2":
        from cwl_utils.cwl_v1_2_expression_refactor import etool_to_cltool
    else:
        print(
            f"Didn't recognise CWL version {version}, loading default: {DEFAULT_PARSER_VERSION}"
        )
        etool_to_cltool = load_etool_to_cltool_from_version(DEFAULT_PARSER_VERSION)
    return etool_to_cltool