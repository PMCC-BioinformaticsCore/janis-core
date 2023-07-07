

import json
import os
from typing import Any, Optional, Tuple
from janis_core import settings
from janis_core import Tool
from janis_core.ingestion.galaxy import runtime
from janis_core.ingestion.galaxy import internal_mapping
from janis_core.ingestion.galaxy.runtime.startup import tool_setup
from janis_core.ingestion.galaxy.gxtool.parsing import load_xmltool
from janis_core.ingestion.galaxy.gxtool.command import gen_command
from janis_core.ingestion.galaxy.containers import resolve_dependencies_as_container

from janis_core.ingestion.galaxy.internal_model.tool.generate import gen_tool
from janis_core.ingestion.galaxy.internal_model.tool import ITool as InternalTool
from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow import StepMetadata
from janis_core.ingestion.galaxy.utils import galaxy as galaxy_utils

from janis_core.ingestion.galaxy.gxworkflow.parsing.metadata import ingest_metadata
from janis_core.ingestion.galaxy.gxworkflow.parsing.inputs import ingest_workflow_inputs
from janis_core.ingestion.galaxy.gxworkflow.parsing.step import ingest_workflow_steps
from janis_core.ingestion.galaxy.gxworkflow.parsing.tool_step.prepost import ingest_workflow_steps_prepost
from janis_core.ingestion.galaxy.gxworkflow.parsing.tool_step.outputs import ingest_workflow_steps_outputs
from janis_core.ingestion.galaxy.gxwrappers import request_single_wrapper

from janis_core.ingestion.galaxy.gxworkflow.values import handle_step_connection_inputs
from janis_core.ingestion.galaxy.gxworkflow.values import handle_step_runtime_inputs
from janis_core.ingestion.galaxy.gxworkflow.values import handle_step_static_inputs
from janis_core.ingestion.galaxy.gxworkflow.values import handle_step_default_inputs

from janis_core.ingestion.galaxy.gxworkflow.updates import update_component_knowledge
from janis_core.ingestion.galaxy.gxworkflow.connections import handle_scattering
from janis_core.ingestion.galaxy.gxworkflow.values.scripts import handle_step_script_configfile_inputs

from janis_core.ingestion.galaxy.gxwrappers import Wrapper
from janis_core.ingestion.galaxy.gxwrappers import WrapperCache
from janis_core.ingestion.galaxy.gxwrappers.downloads.wrappers import get_builtin_tool_path
from janis_core.ingestion.galaxy.gxwrappers.downloads.wrappers import fetch_xml

from janis_core.ingestion.galaxy import datatypes
# from janis_core.ingestion.galaxy.runtime.startup import setup_data_folder
from janis_core.messages import info_ingesting_tool
from janis_core.messages import info_ingesting_workflow
from janis_core.ingestion.galaxy.janis_mapping import to_janis_tool
from janis_core.ingestion.galaxy.janis_mapping import to_janis_workflow


### MODULE ENTRY POINTS ###

def parse_galaxy(uri: str) -> Tool:
    if _is_galaxy_local_tool(uri):
        name = os.path.splitext(uri)[0]
        info_ingesting_tool('galaxy', name)
        tool = ingest_tool(uri)
        return to_janis_tool(tool)
    
    elif _is_galaxy_toolshed_tool(uri):
        wrapper = _request_wrapper_info(uri)
        wrapper_path = fetch_xml(
            owner= wrapper.owner,
            repo= wrapper.repo,
            revision= wrapper.revision,
            tool_id= wrapper.tool_id
        )
        info_ingesting_tool('galaxy', wrapper.tool_id)
        internal_tool = ingest_tool(wrapper_path)
        _set_wrapper_export_paths(wrapper)
        return to_janis_tool(internal_tool)
    
    elif _is_galaxy_workflow(uri):
        name = os.path.splitext(uri)[0]
        info_ingesting_workflow('galaxy', name)
        internal_workflow = ingest_workflow(uri)
        _set_wrapper_export_paths(internal_workflow)
        return to_janis_workflow(internal_workflow)
    
    else:
        raise ValueError("file uri for galaxy ingestion must be either:\n- a tool id (starting with toolshed.g2.bx.psu.edu/)\n- a path ending in '.xml' (local tool)\n - a path ending in '.ga' (local workflow)")
    
def ingest_tool(path: str, gxstep: Optional[dict[str, Any]]=None) -> InternalTool:
    """
    ingests a galaxy tool xml file into a Tool (internal representation).
    'galaxy' is the galaxy tool representation, and
    'internal' is the internal tool representation we will build. 
    """
    # setup_data_folder()
    datatypes.populate()
    runtime.tool.tool_path = path
    galaxy = load_xmltool(path)
    command = gen_command(galaxy, gxstep)
    container = resolve_dependencies_as_container(galaxy)
    internal = gen_tool(galaxy, command, container, gxstep)
    return internal

def ingest_workflow(path: str) -> Workflow:
    """
    Ingests a galaxy workflow file into an *internal* representation. 
    'galaxy' is the actual .ga galaxy workflow, and '*internal*' is a ingestion.galaxy internal representation.
    
    The galaxy workflow doesn't actually map 1-to-1 to the eventual janis_core workflow because of the way tools work in galaxy.
    For this reason, the *internal* representation is kind of a middle-man which we build up over multiple passes of the galaxy workflow. 
    Order seems weird but there is reason for this ordering. 
    Eg: we can't read in the workflow step inputs / outputs until we have parsed the step tool. 
    
    Overall process for galaxy ingest is: galaxy -> *internal* -> janis_core model.
    """
    # setup_data_folder()
    datatypes.populate()
    galaxy = _load_galaxy_workflow(path)
    internal = Workflow()

    # ingesting workflow entities to internal
    ingest_metadata(internal, galaxy)
    ingest_workflow_inputs(internal, galaxy)
    ingest_workflow_steps(internal, galaxy)         # creates steps, but only the metadata
    ingest_workflow_tools(internal, galaxy)         # has to happen after ingesting step metadata, but before step inputs / outputs
    ingest_workflow_steps_prepost(internal, galaxy)
    ingest_workflow_steps_outputs(internal, galaxy) 

    # assigning step input values
    handle_step_connection_inputs(internal, galaxy)   # all these have to happen after step outputs are parsed
    handle_step_runtime_inputs(internal, galaxy)
    handle_step_script_configfile_inputs(internal)
    handle_step_static_inputs(internal, galaxy)
    handle_step_default_inputs(internal)

    # post ingestion tasks
    update_component_knowledge(internal)
    handle_scattering(internal)
    return internal


### HELPER METHODS ###

# (this function should probably be elsewhere)
def ingest_workflow_tools(janis: Workflow, galaxy: dict[str, Any]) -> None:
    for gx_step in galaxy['steps'].values():
        if gx_step['type'] == 'tool':
            j_step = internal_mapping.step(gx_step['id'], janis, galaxy)
            args = _gen_ingest_settings_for_step(j_step.metadata)
            tool_setup(args)
            tool = ingest_tool(runtime.tool.tool_path, gx_step)
            j_step.set_tool(tool)

def _is_galaxy_local_tool(uri: str) -> bool:
    _, ext = os.path.splitext(uri)
    if ext == '.xml':
        return True
    return False
    
def _is_galaxy_toolshed_tool(uri: str) -> bool:
    if uri.startswith('toolshed.g2.bx.psu.edu/'):
        return True
    return False
    
def _is_galaxy_workflow(uri: str) -> bool:
    _, ext = os.path.splitext(uri)
    if ext == '.ga':
        return True
    return False

def _request_wrapper_info(uri: str) -> Wrapper:
    # scrape toolshed for wrappers and update cache
    shed, _, owner, repo, tool_id, tool_build = uri.split('/')
    wrapper = request_single_wrapper(
        tool_shed=shed,
        owner=owner,
        repo=repo,
        tool_id=tool_id,
        tool_build=tool_build
    )
    cache = WrapperCache()
    cache.add(wrapper)
    return wrapper

def _load_galaxy_workflow(path: str) -> dict[str, Any]:
    with open(path, 'r') as fp:
        return json.load(fp)

def _gen_ingest_settings_for_step(metadata: StepMetadata) -> dict[str, Any]:
    tool_id = metadata.wrapper.tool_id
    if metadata.wrapper.inbuilt:
        xml_path = get_builtin_tool_path(tool_id)
        assert(xml_path)
        return {
            'infile': xml_path,
            'remote': None,
            'outdir': None
            #'outdir': f'{paths.wrapper(tool_id, tool_id)}'
        }
    else:
        revision = metadata.wrapper.revision
        owner = metadata.wrapper.owner
        repo = metadata.wrapper.repo
        return {
            'infile': None,
            'remote': f'{owner},{repo},{tool_id},{revision}',
            'outdir': None
            #'outdir': f'{paths.wrapper(tool_id, revision)}'
        }

def _set_wrapper_export_paths(entity: Workflow | Wrapper) -> None:
    settings.general.SOURCE_FILES: list[Tuple[str, str]] = []  # type: ignore

    if isinstance(entity, Workflow):
        for step in entity.steps:
            _set_wrapper_export_path(step.metadata.wrapper)
    elif isinstance(entity, Wrapper): # type: ignore
        _set_wrapper_export_path(entity)
    else:
        raise RuntimeError

def _set_wrapper_export_path(wrapper: Wrapper) -> None:
    assert(isinstance(settings.general.SOURCE_FILES, list))
    dest_dirname = _get_dest_dir(wrapper)
    src_files = _get_wrapper_files_src(wrapper)
    for src in src_files:
        dest = os.path.join(dest_dirname, os.path.basename(src))
        settings.general.SOURCE_FILES.append((src, dest))

def _get_dest_dir(wrapper: Wrapper) -> str:
    if wrapper.revision != 'None':
        return f'{wrapper.tool_id}-{wrapper.revision}'
    return f'{wrapper.tool_id}'

def _get_wrapper_files_src(wrapper: Wrapper) -> list[str]:
    tool_xml = fetch_xml(
        owner= wrapper.owner,
        repo= wrapper.repo,
        revision= wrapper.revision,
        tool_id= wrapper.tool_id
    )
    dirname = os.path.dirname(tool_xml)
    macro_xmls = galaxy_utils.get_macros(dirname)
    xmls = [tool_xml] + macro_xmls
    return sorted(xmls)
