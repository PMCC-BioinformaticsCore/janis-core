



from typing import Protocol
from janis_core.ingestion.galaxy.datatypes import JanisDatatype


class TaggedEntity(Protocol):
    """
    A workflow entity which has a tag.
    Includes Tools, WorkflowSteps, Workflow, WorkflowInput, ToolInputs, ToolOutputs.
    """
    uuid: str
    
    @property
    def name(self) -> str:
        """name of this entity. Derived from held information"""
        ...
    
    @property
    def tag(self) -> str:
        """provides the current tag for this entity"""
        ...



class TypedEntity(Protocol):
    """
    A workflow entity which has a datatype.
    Includes WorkflowInput, WorkflowOutput, ToolInput, ToolOutput
    """

    @property
    def datatype(self) -> JanisDatatype:
        """
        name of this entity. 
        should be derived from information gathered about the entity
        """
        ...
