

from __future__ import annotations
from typing import Optional

from janis_core.ingestion.galaxy.internal_model.tool import ITool
from uuid import uuid4
from janis_core.ingestion.galaxy import tags

from ..registers.StepInputRegister import StepInputRegister
from ..registers.StepOutputRegister import StepOutputRegister
from .metadata import StepMetadata


class WorkflowStep:
    """represents a galaxy tool step"""

    def __init__(self, metadata: StepMetadata):
        self.metadata = metadata
        self.uuid: str = str(uuid4())
        self.inputs: StepInputRegister = StepInputRegister()
        self.outputs: StepOutputRegister = StepOutputRegister()
        self.preprocessing: Optional[str] = None
        self.postprocessing: Optional[str] = None
        self._tool: Optional[ITool] = None

    @property
    def name(self) -> str:
        return self.metadata.wrapper.tool_id
    
    @property
    def tag(self) -> str:
        return tags.get(self.uuid)
    
    @property
    def docstring(self) -> Optional[str]:
        return self.metadata.label

    @property
    def tool(self) -> ITool:
        if self._tool:
            return self._tool
        raise RuntimeError()

    def set_tool(self, tool: ITool) -> None:
        self._tool = tool

    def __repr__(self) -> str:
        return f'(WorkflowStep) step{self.metadata.step_id} - {self.name}'
