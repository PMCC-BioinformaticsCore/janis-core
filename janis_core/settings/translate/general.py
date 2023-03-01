

from typing import Optional, Any
from janis_core.translation_deps.exportpath import ExportPathKeywords



EXPORT_PATH:                str = ExportPathKeywords.default # base output directory

STRICT_IDENTIFIERS:         bool = True  # whether to enforce rules about tool / workflow / tool input etc identifiers
ALLOW_EMPTY_CONTAINER:      bool = True  # makes docker containers optional
MERGE_RESOURCES:            bool = False # merge resource requirements into inputs config
RENDER_COMMENTS:            bool = True  # whether to render info comments in the translation
SHOULD_VALIDATE:            bool = False # whether to validate translated files
SHOULD_ZIP:                 bool = False # whether to zip translated tool folder
TO_DISK:                    bool = False # whether to write translated files to disk
TO_CONSOLE:                 bool = True  # whether to write main translated file to console
TOOL_TO_CONSOLE:            bool = False # whether to write translated tool files to console 
                                         # (workflow translation only)
WITH_CONTAINER:             bool = True  # with_container=True, with_docker=True,
WITH_RESOURCE_OVERRIDES:    bool = False # whether to add computational resources to inputs dict. 
                                         # uses resources specified in ingested workflow, else default values.
WRITE_INPUTS_FILE:          bool = True  # whether to write the inputs file to disk

CONTAINER_OVERRIDE:         Optional[str | dict[str, Any]] = None # val, or key, val map where keys are tool ids, vals are containers to use
ADDITIONAL_INPUTS:          Optional[dict[str, Any]] = None # key, val map for additional tool / workflow inputs supplied by user
HINTS:                      Optional[dict[str, Any]] = None # key, val map for cwl type hints
MAX_CORES:                  Optional[int] = None            # ceiling value for cores resource
MAX_DURATION:               Optional[int] = None            # ceiling value for duration resource
MAX_MEM:                    Optional[int] = None            # ceiling value for memory resource

