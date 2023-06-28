
from dataclasses import dataclass
from typing import Any, Optional

from janis_core.ingestion.galaxy.gx.wrappers import Wrapper


@dataclass
class StepMetadata:
    wrapper: Wrapper
    step_id: int
    step_name: str
    tool_state: dict[str, Any]
    workflow_outputs: list[dict[str, Any]]
    _label: Optional[str] = None

    @property
    def label(self) -> str:
        if self._label:
            return self._label
        return self.step_name

