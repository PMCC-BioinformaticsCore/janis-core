


from janis_core.ingestion.galaxy.model.tool import Tool
from janis_core.ingestion.galaxy.model.workflow import Workflow
from janis_core.ingestion.galaxy.model.workflow import WorkflowStep
from .Format import Format


class SubworkflowFormat(Format):
    """
    renders output pages in different manner as appears in original workflow.
    steps using tools which have pre/post process are broken into 3 step subworkflows. 
    """

    def workflow(self, workflow: Workflow) -> str:
        raise NotImplementedError()
    
    def subworkflow(self, workflow: Workflow) -> str:
        raise NotImplementedError()
    
    def step(self, step: WorkflowStep) -> str:
        raise NotImplementedError()

    def tool(self, tool: Tool) -> str:
        raise NotImplementedError()

    def script(self) -> str:
        raise NotImplementedError()
    
    def inputs(self, workflow: Workflow) -> str:
        raise NotImplementedError()
    
    def config(self, workflow: Workflow) -> str:
        raise NotImplementedError()


