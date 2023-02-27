
import os
import ruamel.yaml
from typing import Any, Optional

DEFAULT_PARSER_VERSION = "v1.2"


def load_cwl_document(doc: str, base_uri: Optional[str]=None) -> Any:
    """loads a cwl document & returns the in-memory cwlutils object"""
    initial_wd = os.getcwd()
    
    if base_uri:
        if base_uri.startswith("file://"):
            base_uri = base_uri[6:]
        os.chdir(base_uri)
    
    version = get_cwl_version_from_doc(doc)
    cwlgen = load_cwlgen_from_version(version)
    loaded_doc = cwlgen.load_document(doc)  # type: ignore
    clsname = loaded_doc.__class__.__name__  # type: ignore
    if clsname == 'ExpressionTool':
        loaded_doc = convert_etool_to_cltool(loaded_doc, version)

    if base_uri:
        os.chdir(initial_wd)
    return loaded_doc

def convert_etool_to_cltool(etool: Any, version: str) -> Any:
    etool_to_cltool = load_etool_to_cltool_from_version(version)
    cltool = etool_to_cltool(etool)
    for out in cltool.outputs:  # type: ignore
        out_id = out.id.split(".")[-1]  # type: ignore
        out.outputEval = f"JANIS (potentially unimplemented): j.ReadJsonOperator(j.Stdout)[{out_id}]"
    return cltool

def get_cwl_version_from_doc(doc: str) -> str:
    # load tool into memory
    if doc.startswith("file://"):
        doc = doc[7:]
    with open(doc) as fp:
        tool_dict = ruamel.yaml.load(fp, Loader=ruamel.yaml.Loader) # type: ignore
    if "cwlVersion" not in tool_dict:
        raise Exception(f"Couldn't find cwlVersion in tool {doc}")
    return tool_dict["cwlVersion"]

def load_cwlgen_from_version(cwl_version: str) -> Any:
    if cwl_version == "v1.0":
        import cwl_utils.parser.cwl_v1_0 as cwlutils
    elif cwl_version == "v1.1":
        import cwl_utils.parser.cwl_v1_1 as cwlutils
    elif cwl_version == "v1.2":
        import cwl_utils.parser.cwl_v1_2 as cwlutils
    else:
        print(
            f"Didn't recognise CWL version {cwl_version}, loading default: {DEFAULT_PARSER_VERSION}"
        )
        cwlutils = load_cwlgen_from_version(DEFAULT_PARSER_VERSION)
    return cwlutils

def load_etool_to_cltool_from_version(cwl_version: str) -> Any:
    if cwl_version == "v1.0":
        from cwl_utils.cwl_v1_0_expression_refactor import etool_to_cltool
    elif cwl_version == "v1.1":
        from cwl_utils.cwl_v1_0_expression_refactor import etool_to_cltool
    elif cwl_version == "v1.2":
        from cwl_utils.cwl_v1_2_expression_refactor import etool_to_cltool
    else:
        print(
            f"Didn't recognise CWL version {cwl_version}, loading default: {DEFAULT_PARSER_VERSION}"
        )
        etool_to_cltool = load_etool_to_cltool_from_version(DEFAULT_PARSER_VERSION)
    return etool_to_cltool