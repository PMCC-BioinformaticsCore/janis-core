

import json
import os
from typing import Any, Optional, Tuple
from janis_core import settings
from janis_core import Tool
from janis_core.ingestion.galaxy import runtime
from janis_core.ingestion.galaxy import mapping
from janis_core.ingestion.galaxy.startup import tool_setup
from janis_core.ingestion.galaxy.gx.gxtool import load_xmltool
from janis_core.ingestion.galaxy.gx.command import gen_command
from janis_core.ingestion.galaxy.containers import resolve_dependencies_as_container

from janis_core.ingestion.galaxy.model.tool.generate import gen_tool
from janis_core.ingestion.galaxy.model.tool import Tool as InternalTool
from janis_core.ingestion.galaxy.model.workflow import Workflow
from janis_core.ingestion.galaxy.model.workflow import WorkflowStep
from janis_core.ingestion.galaxy.model.workflow import StepMetadata
from janis_core.ingestion.galaxy.utils import galaxy as galaxy_utils

from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.metadata import ingest_metadata
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.inputs import ingest_workflow_inputs
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.step import ingest_workflow_steps
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.tool_step.prepost import ingest_workflow_steps_prepost
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.tool_step.outputs import ingest_workflow_steps_outputs
from janis_core.ingestion.galaxy.gx.wrappers.requests.versions import request_single_wrapper

from janis_core.ingestion.galaxy.gx.gxworkflow.values import handle_step_connection_inputs
from janis_core.ingestion.galaxy.gx.gxworkflow.values import handle_step_runtime_inputs
from janis_core.ingestion.galaxy.gx.gxworkflow.values import handle_step_static_inputs
from janis_core.ingestion.galaxy.gx.gxworkflow.values import handle_step_default_inputs

from janis_core.ingestion.galaxy.gx.gxworkflow.updates import update_component_knowledge
from janis_core.ingestion.galaxy.gx.gxworkflow.connections import handle_scattering
from janis_core.ingestion.galaxy.gx.gxworkflow.values.scripts import handle_step_script_configfile_inputs

from janis_core.ingestion.galaxy.gx.wrappers import Wrapper
from janis_core.ingestion.galaxy.gx.wrappers import WrapperCache
from janis_core.ingestion.galaxy.gx.wrappers.downloads.wrappers import get_builtin_tool_path
from janis_core.ingestion.galaxy.gx.wrappers.downloads.wrappers import fetch_wrapper

from janis_core.ingestion.galaxy import datatypes
from janis_core.ingestion.galaxy.startup import setup_data_folder
from janis_core.messages import info_ingesting_tool
from janis_core.messages import info_ingesting_workflow
from janis_core.ingestion.galaxy.janis_mapping import to_janis_tool
from janis_core.ingestion.galaxy.janis_mapping import to_janis_workflow


def parse_galaxy(uri: str) -> Tool:
    if _is_galaxy_local_tool(uri):
        name = os.path.splitext(uri)[0]
        info_ingesting_tool('galaxy', name)
        tool = ingest_tool(uri)
        return to_janis_tool(tool)
    
    elif _is_galaxy_toolshed_tool(uri):
        wrapper = _request_wrapper_info(uri)
        wrapper_path = fetch_wrapper(
            owner= wrapper.owner,
            repo= wrapper.repo,
            revision= wrapper.revision,
            tool_id= wrapper.tool_id
        )
        info_ingesting_tool('galaxy', wrapper.tool_id)
        tool = ingest_tool(wrapper_path)
        return to_janis_tool(tool)
    
    elif _is_galaxy_workflow(uri):
        name = os.path.splitext(uri)[0]
        info_ingesting_workflow('galaxy', name)
        workflow = ingest_workflow(uri)
        return to_janis_workflow(workflow)
    
    else:
        raise ValueError("file uri for galaxy ingestion must be either:\n- a tool id (starting with toolshed.g2.bx.psu.edu/)\n- a path ending in '.xml' (local tool)\n - a path ending in '.ga' (local workflow)")
    
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

def ingest_tool(path: str, gxstep: Optional[dict[str, Any]]=None) -> InternalTool:
    """
    ingests a galaxy tool xml file into a Tool (internal representation).
    'galaxy' is the galaxy tool representation, and
    'internal' is the internal tool representation we will build. 
    """
    setup_data_folder()
    datatypes.populate()
    runtime.tool.tool_path = path
    galaxy = load_xmltool(path)
    command = gen_command(galaxy)
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
    setup_data_folder()
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
    handle_source_files(internal)
    return internal

def ingest_workflow_tools(janis: Workflow, galaxy: dict[str, Any]) -> None:
    for g_step in galaxy['steps'].values():
        if g_step['type'] == 'tool':
            j_step = mapping.step(g_step['id'], janis, galaxy)
            tool = _parse_step_tool(j_step.metadata, g_step)
            j_step.set_tool(tool)

def _load_galaxy_workflow(path: str) -> dict[str, Any]:
    with open(path, 'r') as fp:
        return json.load(fp)

def _parse_step_tool(metadata: StepMetadata, gxstep: dict[str, Any]) -> Tool:
    args = _create_tool_settings_for_step(metadata)
    tool_setup(args)
    # gxstep['tool_state'] = load_tool_state(
    #     gxstep, 
    #     additional_filters=[
    #         'ReplaceNullWithVarname'
    #         'ReplaceConnectedWithVarname',
    #         'ReplaceRuntimeWithVarname',
    #     ]
    # )
    tool = ingest_tool(runtime.tool.tool_path, gxstep)
    return tool

def _create_tool_settings_for_step(metadata: StepMetadata) -> dict[str, Any]:
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

def handle_source_files(janis: Workflow) -> None:
    settings.general.SOURCE_FILES: list[Tuple[str, str]] = []  # type: ignore
    for step in janis.steps:
        dest_dir = get_dest_dir(step)
        src_files = get_wrapper_files_src(step)
        dest_files = get_wrapper_files_dest(src_files, dest_dir)
        for src, dest in zip(src_files, dest_files):
            settings.general.SOURCE_FILES.append((src, dest))

def get_dest_dir(step: WorkflowStep) -> str:
    if step.metadata.wrapper.revision != 'None':
        return f'{step.metadata.wrapper.tool_id}-{step.metadata.wrapper.revision}'
    return f'{step.metadata.wrapper.tool_id}'

def get_wrapper_files_src(step: WorkflowStep) -> list[str]:
    wrapper = step.metadata.wrapper
    wrapper_path = fetch_wrapper(
        owner= wrapper.owner,
        repo= wrapper.repo,
        revision= wrapper.revision,
        tool_id= wrapper.tool_id
    )
    wrapper_dir = wrapper_path.rsplit('/', 1)[0]
    macro_xmls = galaxy_utils.get_macros(wrapper_dir)
    xmls = [wrapper_path] + macro_xmls
    return sorted(xmls)

def get_wrapper_files_dest(src_files: list[str], dest_dir: str) -> list[str]:
    out: list[str] = []
    for src in src_files:
        xmlpath = src.rsplit('/', 1)[-1]
        out.append(os.path.join(dest_dir, xmlpath))
    return sorted(out)