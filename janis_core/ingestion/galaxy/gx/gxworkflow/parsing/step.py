

from typing import Any

from janis_core.ingestion.galaxy.model.workflow import Workflow
from janis_core.ingestion.galaxy.model.workflow import WorkflowStep

from .tool_step.metadata import parse_step_metadata


def ingest_workflow_steps(janis: Workflow, galaxy: dict[str, Any]) -> None:
    for step in galaxy['steps'].values():
        if step['type'] == 'tool':
            workflow_step = parse_tool_step(step)
            step['janis_uuid'] = workflow_step.uuid # crucial!
            janis.add_step(workflow_step)

def parse_tool_step(step: dict[str, Any]) -> WorkflowStep:
    metadata = parse_step_metadata(step)
    return WorkflowStep(metadata=metadata)