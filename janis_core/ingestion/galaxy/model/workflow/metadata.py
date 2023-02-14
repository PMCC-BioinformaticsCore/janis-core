


from dataclasses import dataclass


@dataclass
class WorkflowMetadata:
    name: str
    uuid: str  
    annotation: str
    version: str
    tags: list[str]