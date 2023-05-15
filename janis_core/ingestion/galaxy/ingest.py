

import json
from typing import Any, Optional
from janis_core.ingestion.galaxy import runtime
from janis_core.ingestion.galaxy import mapping

from janis_core.ingestion.galaxy.startup import tool_setup
from janis_core.ingestion.galaxy.gx.gxtool import load_xmltool
from janis_core.ingestion.galaxy.gx.command import gen_command
from janis_core.ingestion.galaxy.containers import resolve_dependencies_as_container

from janis_core.ingestion.galaxy.model.tool.generate import gen_tool
from janis_core.ingestion.galaxy.model.tool import Tool
from janis_core.ingestion.galaxy.model.workflow import Workflow
from janis_core.ingestion.galaxy.model.workflow import StepMetadata

from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.metadata import ingest_metadata
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.inputs import ingest_workflow_inputs
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.step import ingest_workflow_steps
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.tool_step.prepost import ingest_workflow_steps_prepost
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.tool_step.outputs import ingest_workflow_steps_outputs
from janis_core.ingestion.galaxy.gx.gxworkflow.parsing.tool_state import load_tool_state

from janis_core.ingestion.galaxy.gx.gxworkflow.values import handle_step_connection_inputs
from janis_core.ingestion.galaxy.gx.gxworkflow.values import handle_step_runtime_inputs
from janis_core.ingestion.galaxy.gx.gxworkflow.values import handle_step_static_inputs
from janis_core.ingestion.galaxy.gx.gxworkflow.values import handle_step_default_inputs

from janis_core.ingestion.galaxy.gx.gxworkflow.updates import update_component_knowledge
from janis_core.ingestion.galaxy.gx.gxworkflow.connections import handle_scattering
from janis_core.ingestion.galaxy.gx.gxworkflow.values.scripts import handle_step_script_inputs

from janis_core.ingestion.galaxy.gx.wrappers.downloads.wrappers import get_builtin_tool_path

from janis_core.ingestion.galaxy import datatypes
from janis_core.ingestion.galaxy.startup import setup_data_folder

# TODO future 
# from janis_core.ingestion.galaxy.gx.xmltool.tests import write_tests


def ingest_tool(path: str, gxstep: Optional[dict[str, Any]]=None) -> Tool:
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
    ingests a galaxy workflow file into a Workflow (internal representation).
    'galaxy' is the galaxy workflow representation, and
    'internal' is the internal workflow representation we will build. 
    order seems weird but trust me there is reason for this ordering.
    """
    setup_data_folder()
    datatypes.populate()
    galaxy = _load_galaxy_workflow(path)
    internal = Workflow()

    # ingesting workflow entities to internal
    ingest_metadata(internal, galaxy)
    ingest_workflow_inputs(internal, galaxy)
    ingest_workflow_steps(internal, galaxy)
    ingest_workflow_tools(internal, galaxy)
    ingest_workflow_steps_prepost(internal, galaxy)
    ingest_workflow_steps_outputs(internal, galaxy) 

    # assigning step input values
    handle_step_connection_inputs(internal, galaxy)
    handle_step_runtime_inputs(internal, galaxy)
    handle_step_script_inputs(internal)
    handle_step_static_inputs(internal, galaxy)
    handle_step_default_inputs(internal)

    # post ingestion tasks
    update_component_knowledge(internal)
    handle_scattering(internal)
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
    gxstep['tool_state'] = load_tool_state(gxstep)
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
