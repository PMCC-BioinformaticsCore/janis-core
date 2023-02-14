

from typing import Any

from janis_core.ingestion.galaxy.model.workflow import StepOutput
from janis_core.ingestion.galaxy.model.workflow import WorkflowStep
from janis_core.ingestion.galaxy.model.workflow import Workflow

from janis_core.ingestion.galaxy import mapping
from janis_core.ingestion.galaxy import settings


def ingest_workflow_steps_outputs(janis: Workflow, galaxy: dict[str, Any]) -> None:
    ingestor = OutputParser(janis, galaxy)
    ingestor.parse()


class OutputParser:
    def __init__(self, janis: Workflow, galaxy: dict[str, Any]):
        self.janis = janis
        self.galaxy = galaxy
    
    def parse(self) -> None:
        for step in self.galaxy['steps'].values():
            if step['type'] == 'tool':
                self.parse_step_outputs(step)

    def parse_step_outputs(self, galaxy_step: dict[str, Any]) -> None:
        janis_step = mapping.step(galaxy_step['id'], self.janis, self.galaxy)
        settings.tool.set(from_wrapper=janis_step.metadata.wrapper)
        for galaxy_out in galaxy_step['outputs']:
            janis_out = init_tool_step_output(janis_step, galaxy_step, galaxy_out)
            galaxy_out['janis_uuid'] = janis_out.uuid # entity linking
            janis_step.outputs.add(janis_out)


def init_tool_step_output(janis_step: WorkflowStep, step: dict[str, Any], output: dict[str, Any]) -> StepOutput:
    tool_output = mapping.tool_output(output['name'], janis_step.tool)
    assert(tool_output)
    return StepOutput(
        step_uuid=janis_step.uuid,
        is_wflow_out=is_workflow_output(step, output),
        tool_output=tool_output,
    )

def is_workflow_output(step: dict[str, Any], output: dict[str, Any]) -> bool:
    if step['workflow_outputs']:
        for details in step['workflow_outputs']:
            if details['output_name'] == output['name']:
                return True
    elif not step['workflow_outputs']:
        return True
    return False