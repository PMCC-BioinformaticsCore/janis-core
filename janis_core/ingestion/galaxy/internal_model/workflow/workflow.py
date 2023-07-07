

from typing import Optional
from uuid import uuid4
from janis_core.ingestion.galaxy import tags

from .step.outputs import StepOutput
from .step.step import WorkflowStep
from .metadata import WorkflowMetadata
from .input import WorkflowInput



class Workflow:
    """
    a Workflow() is the final representation of the galaxy workflow
    being parsed. Includes workflow metadata, inputs, tool steps, and outputs. 
    step_ids are used to lookup key information. 
    each input, step, and output is stored with step_id as key. 
    that way everything can be referenced properly. 
    """
    def __init__(self):
        self.uuid: str = str(uuid4())
        self.inputs: list[WorkflowInput] = []
        self._steps: list[WorkflowStep] = []
        self._metadata: Optional[WorkflowMetadata] = None

    @property
    def name(self) -> str:
        return self.metadata.name
    
    @property
    def tag(self) -> str:
        return tags.get(self.uuid)
    
    @property
    def steps(self) -> list[WorkflowStep]:
        the_list = self._steps
        the_list = sorted(the_list, key=lambda x: x.metadata.step_id)
        return the_list

    @property
    def metadata(self) -> WorkflowMetadata:
        if self._metadata:
            return self._metadata
        raise RuntimeError('no metadata set')

    @property
    def outputs(self) -> list[StepOutput]:
        workflow_outputs: list[StepOutput] = []
        for step in self.steps:
            for out in step.outputs.list():
                if out.is_wflow_out:
                    workflow_outputs.append(out)
        return workflow_outputs

    def set_metadata(self, metadata: WorkflowMetadata) -> None:
        self._metadata = metadata
        tags.new_group('workflow', self.uuid)
        tags.register(self)
    
    def add_input(self, w_inp: WorkflowInput) -> None:
        tags.switch_group(self.uuid)
        tags.register(w_inp)
        self.inputs.append(w_inp)

    def add_step(self, step: WorkflowStep) -> None:
        tags.switch_group(self.uuid)
        tags.register(step)
        self._steps.append(step)

    def get_input(self, query_uuid: str) -> WorkflowInput:
        for winp in self.inputs:
            if winp.uuid == query_uuid:
                return winp
        raise RuntimeError('could not find input with uuid')
    
    def get_input_children(self, query_uuid: str) -> list[WorkflowStep]:
        out: list[WorkflowStep] = []
        for step in self.steps:
            for winp_connection in step.inputs.workflow_inputs:
                if winp_connection.input_uuid == query_uuid:
                    out.append(step)
        return out
    
    def get_step(self, query_uuid: str) -> WorkflowStep:
        for step in self.steps:
            if step.uuid == query_uuid:
                return step
        raise RuntimeError('could not find step with uuid')
    
    def get_step_children(self, query_uuid: str) -> list[WorkflowStep]:
        out: list[WorkflowStep] = []
        for step in self.steps:
            for connection in step.inputs.connections:
                if connection.step_uuid == query_uuid:
                    out.append(step)
        return out
    
    def get_step_output(self, query_uuid: str) -> StepOutput:
        for step in self.steps:
            for s_out in step.outputs.list():
                if s_out.tool_output.uuid == query_uuid:
                    return s_out
        raise RuntimeError('could not find step output with uuid')
    

