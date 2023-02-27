

"""
to_console=True,
tool_to_console=False,
with_resource_overrides=False,
to_disk=False,
write_inputs_file=True,
export_path=ExportPathKeywords.default,
should_validate=False,
should_zip=True,   
merge_resources=False,
hints=None,
allow_null_if_not_optional=True,
additional_inputs: Dict = None,
max_cores=None,
max_mem=None,
max_duration=None,
with_container=True,
allow_empty_container=False,
container_override=None,
render_comments: bool = True
with_docker=True,
with_hints=False,
"""

from typing import Optional, Any
from janis_core.translation_deps.exportpath import ExportPathKeywords



EXPORT_PATH:                str = ExportPathKeywords.default # base output directory

STRICT_IDENTIFIERS:         bool = True  # whether to enforce rules about tool / workflow / tool input etc identifiers
ALLOW_EMPTY_CONTAINER:      bool = False # makes docker containers optional
MERGE_RESOURCES:            bool = False # merge resource requirements into inputs config
RENDER_COMMENTS:            bool = True  # whether to render info comments in the translation
SHOULD_VALIDATE:            bool = False # whether to validate translated files
SHOUD_ZIP:                  bool = False # whether to zip translated tool folder
TO_DISK:                    bool = False # whether to write translated files to disk
TO_CONSOLE:                 bool = True  # whether to write main translated file to console
TOOL_TO_CONSOLE:            bool = False # whether to write translated tool files to console 
                                         # (workflow translation only)
WITH_CONTAINER:             bool = True  # with_container=True, with_docker=True,
WITH_RESOURCE_OVERRIDES:    bool = False # whether to add computational resources to inputs dict. 
                                         # uses resources specified in ingested workflow, else default values.
WRITE_INPUTS_FILE:          bool = True  # whether to write the inputs file to disk

ADDITIONAL_INPUTS:          Optional[dict[str, Any]] = None # key, val map for additional tool / workflow inputs supplied by user
CONTAINER_OVERRIDES:        Optional[dict[str, Any]] = None # key, val map where keys are tool ids, vals are containers to use
HINTS:                      Optional[dict[str, Any]] = None # key, val map for cwl type hints
MAX_CORES:                  Optional[int] = None            # ceiling value for cores resource
MAX_DURATION:               Optional[int] = None            # ceiling value for duration resource
MAX_MEM:                    Optional[int] = None            # ceiling value for memory resource

