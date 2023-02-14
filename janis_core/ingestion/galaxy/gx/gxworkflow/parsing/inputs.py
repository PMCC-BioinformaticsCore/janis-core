

import json
from typing import Any

from janis_core.ingestion.galaxy.model.workflow import WorkflowInput
from janis_core.ingestion.galaxy.model.workflow import Workflow

from janis_core.ingestion.galaxy import datatypes


def ingest_workflow_inputs(janis: Workflow, galaxy: dict[str, Any]) -> None:
    for step in galaxy['steps'].values():
        if step['type'] in ['data_input', 'data_collection_input']:
            workflow_input = parse_input_step(step)
            step['janis_uuid'] = workflow_input.uuid # crucial!
            janis.add_input(workflow_input)

def parse_input_step(step: dict[str, Any]) -> WorkflowInput:
    return WorkflowInput(
        _name=format_input_step_name(step),
        array=format_input_step_array(step),
        is_runtime=False,
        datatype=datatypes.get(step, 'GalaxyInputStep'),
        optional=format_input_step_optionality(step)
    )

def format_input_step_name(step: dict[str, Any]) -> str:
    if step['label'] and step['label'] != '':
        return step['label']
    elif step['inputs']:
        return step['inputs'][0]['name']
    else:
        raise RuntimeError()
    
def format_input_step_array(step: dict[str, Any]) -> bool:
    if step['type'] and step['type'] == 'data_collection_input':
        return True
    return False

def format_input_step_optionality(step: dict[str, Any]) -> bool:
    tool_state = json.loads(step['tool_state'])
    if 'optional' in tool_state:
        return bool(tool_state['optional'])
    return False
