

from abc import ABC, abstractmethod

from janis_core.ingestion.galaxy.internal_model.tool import ITool
from janis_core.ingestion.galaxy.internal_model.workflow import Workflow
from janis_core.ingestion.galaxy.internal_model.workflow import WorkflowStep


class Format(ABC):
    """
    controls the format of pages user will see.
    example 1: controls whether the workflow steps will be 
    broken into subworkflows or not. 
    example 2?: controls which workflows will have an inputs dict file created? 
    """

    @abstractmethod
    def workflow(self, workflow: Workflow) -> str:
        """returns the workflow page"""
        ...
    
    @abstractmethod
    def subworkflow(self, workflow: Workflow) -> str:
        """returns a subworkflow"""
        ...
    
    @abstractmethod
    def step(self, step: WorkflowStep) -> str:
        """writes a step"""
        ...

    @abstractmethod
    def tool(self, tool: ITool) -> str:
        """writes a tool"""
        ...

    @abstractmethod
    def script(self) -> str:
        """writes a script used with a tool"""
        ...
    
    @abstractmethod
    def inputs(self, workflow: Workflow) -> str:
        """writes the workflow inputs file"""
        ...
    
    @abstractmethod
    def config(self, workflow: Workflow) -> str:
        """writes the workflow config file"""
        ...


