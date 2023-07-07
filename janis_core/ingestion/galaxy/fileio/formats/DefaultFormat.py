


from janis_core.ingestion.galaxy.internal_model.tool import ITool
from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowStep

from .Format import Format


class DefaultFormat(Format):
    """renders output pages in same manner as appears in original workflow"""

    def workflow(self, workflow: Workflow) -> str:
        raise NotImplementedError()
    
    def subworkflow(self, workflow: Workflow) -> str:
        raise NotImplementedError()
    
    def step(self, step: WorkflowStep) -> str:
        raise NotImplementedError()

    def tool(self, tool: ITool) -> str:
        raise NotImplementedError()

    def script(self) -> str:
        raise NotImplementedError()
    
    def inputs(self, workflow: Workflow) -> str:
        raise NotImplementedError()
    
    def config(self, workflow: Workflow) -> str:
        raise NotImplementedError()


