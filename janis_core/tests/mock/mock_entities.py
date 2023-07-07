

import json

from janis_core.ingestion.galaxy.gxworkflow.parsing.tool_step.metadata import parse_step_metadata
from janis_core.ingestion.galaxy.gxworkflow.parsing.inputs import parse_input_step

from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowInputInputValue
from janis_core.ingestion.galaxy.internal_model.workflow import StaticInputValue
from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowStep
from janis_core.ingestion.galaxy.internal_model.workflow import StepOutput

from .mock_tool import MOCK_TOOL_ABRICATE


gx_workflow_path = 'janis_core/tests/data/galaxy/wf_abricate.ga'  
step1_wrapper_details = {
    "owner": "iuc", 
    "repo": "abricate", 
    "revision": "c2ef298da409", 
    "tool_id": "abricate", 
    "tool_version": "1.0.1", 
    "tool_build": "1.0.1", 
    "date_created": "2014-01-27 14:29:14", 
    "requirements": [{"type": "conda", "name": "abricate", "version": "1.0.1"}], 
    "url": "https://toolshed.g2.bx.psu.edu/repos/iuc/abricate/archive/c2ef298da409.tar.gz"
}

with open(gx_workflow_path, 'r') as fp:
    wf = json.load(fp)

MOCK_WORKFLOW_INPUT1 = parse_input_step(wf['steps']['0'])
MOCK_STEP1 = WorkflowStep(parse_step_metadata(wf['steps']['1'])) # avoids wrapper cache lookup
MOCK_STEP1.set_tool(MOCK_TOOL_ABRICATE)

MOCK_STEP_INPUT1 = WorkflowInputInputValue(
    component=MOCK_TOOL_ABRICATE.inputs[0], # MOCK_POSITIONAL1
    input_uuid=MOCK_WORKFLOW_INPUT1.uuid,
    is_runtime=False,
)

MOCK_STEP_INPUT2 = StaticInputValue(
    component=MOCK_TOOL_ABRICATE.inputs[1], # MOCK_FLAG1
    str_value='False',
    is_default=False
)

MOCK_STEP_INPUT3 = StaticInputValue(
    component=MOCK_TOOL_ABRICATE.inputs[2], # MOCK_OPTION1
    str_value='[80.0, 100.0]',
    is_default=False
)
MOCK_STEP_INPUT3.scatter = True

MOCK_STEP_INPUT4 = StaticInputValue(
    component=MOCK_TOOL_ABRICATE.inputs[3], # MOCK_OPTION2
    str_value='card',
    is_default=False
)

MOCK_STEP1.inputs.add(MOCK_STEP_INPUT1)
MOCK_STEP1.inputs.add(MOCK_STEP_INPUT2)
MOCK_STEP1.inputs.add(MOCK_STEP_INPUT3)
MOCK_STEP1.inputs.add(MOCK_STEP_INPUT4)

MOCK_STEP_OUTPUT1 = StepOutput(
    step_uuid=MOCK_STEP1.uuid,
    is_wflow_out=True,
    tool_output=MOCK_TOOL_ABRICATE.outputs[0]
)
MOCK_STEP1.outputs.add(MOCK_STEP_OUTPUT1)
