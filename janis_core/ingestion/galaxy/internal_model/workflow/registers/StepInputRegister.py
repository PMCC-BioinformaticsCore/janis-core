

from typing import Optional

from ..step.inputs import InputValue, WorkflowInputInputValue
from ..step.inputs import ConnectionInputValue



class StepInputRegister:
    def __init__(self):
        self.inputs: list[InputValue] = []

    @property
    def all(self) -> list[InputValue]:
        return self.inputs
    
    @property
    def connections(self) -> list[ConnectionInputValue]:
        return [x for x in self.inputs if isinstance(x, ConnectionInputValue)]
    
    @property
    def workflow_inputs(self) -> list[WorkflowInputInputValue]:
        return [x for x in self.inputs if isinstance(x, WorkflowInputInputValue)]
    
    @property
    def linked(self) -> list[InputValue]:
        return [x for x in self.inputs if x.component]
    
    @property
    def unlinked(self) -> list[InputValue]:
        return [x for x in self.inputs if not x.component]

    def add(self, invalue: InputValue) -> None:
        self.inputs.append(invalue)
    
    def get(self, query_uuid: str) -> Optional[InputValue]:
        for invalue in self.inputs:
            if invalue.component and invalue.component.uuid == query_uuid:
                return invalue
