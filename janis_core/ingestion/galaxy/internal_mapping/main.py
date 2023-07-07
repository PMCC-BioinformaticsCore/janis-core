

from __future__ import annotations
from typing import TYPE_CHECKING, Any, Optional

if TYPE_CHECKING:
    from janis_core.ingestion.galaxy.internal_model.tool.tool import ITool
    from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
    from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowInput
    from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowStep
    from janis_core.ingestion.galaxy.internal_model.workflow import StepOutput

    from janis_core.ingestion.galaxy.gxtool.command.components import InputComponent
    from janis_core.ingestion.galaxy.gxtool.command.components import OutputComponent


# split this file into mapping.workflow, mapping.tool?

def tool_input(galaxy_param_name: str, tool: ITool) -> Optional[InputComponent]:
    """given a galaxy param name, get the corresponding Tool InputComponent if exists"""
    for inp in tool.inputs:
        if inp.gxparam and inp.gxparam.name == galaxy_param_name:
            return inp

def tool_output(galaxy_param_name: str, tool: ITool) -> Optional[OutputComponent]:
    """given a galaxy param name, get the corresponding Tool OutputComponent if exists"""
    for out in tool.outputs:
        if out.gxparam and out.gxparam.name == galaxy_param_name:
            return out

def emitter(query_id: str, query_output: str, internal: Workflow, galaxy: dict[str, Any]) -> StepOutput | WorkflowInput:
    """given a galaxy step id and output name, get the corresponding internal StepOutput or WorkflowInput"""
    g_step = _get_galaxy_step(query_id, galaxy)
    if g_step['type'] == 'tool':
        i_step = step(g_step['id'], internal, galaxy)
        return _get_internal_step_output(query_output, i_step, g_step)
    else:
        return _get_internal_workflow_input(g_step['id'], internal, galaxy)

def step(query_id: str, internal: Workflow, galaxy: dict[str, Any]) -> WorkflowStep:
    """lookup the corresponding internal step for a galaxy step with id=query_id"""
    g_step = _get_galaxy_step(query_id, galaxy)
    if 'janis_uuid' not in g_step:
        raise RuntimeError('Galaxy step not linked to internal step')
    internal_step_uuid = g_step['janis_uuid']
    return _get_internal_step(internal_step_uuid, internal)

def _get_galaxy_step(query_id: str, galaxy: dict[str, Any]) -> dict[str, Any]:
    for step in galaxy['steps'].values():
        if step['id'] == query_id:
            return step
    raise RuntimeError(f'No galaxy step with id {query_id}')

def _get_internal_step(query_uuid: str, internal: Workflow) -> WorkflowStep:
    for step in internal.steps:
        if step.uuid == query_uuid:
            return step
    raise RuntimeError(f'No internal step with uuid {query_uuid}')

def _get_internal_step_output(query_output: str, j_step: WorkflowStep, g_step: dict[str, Any]) -> StepOutput:
    """sorry"""
    j_outs = j_step.outputs.list()
    g_outs = g_step['outputs']
    for out in g_outs:
        if out['name'] == query_output:
            for output in j_outs:
                if output.uuid == out['janis_uuid']:
                    return output
    raise RuntimeError(f'No internal step output for galaxy step output {query_output}')

def _get_internal_workflow_input(query_id: str, internal: Workflow, galaxy: dict[str, Any]) -> WorkflowInput:
    g_step = _get_galaxy_step(query_id, galaxy)
    for inp in internal.inputs:
        if inp.uuid == g_step['janis_uuid']:
            return inp
    raise RuntimeError(f'No internal workflow input with uuid {query_id}')



