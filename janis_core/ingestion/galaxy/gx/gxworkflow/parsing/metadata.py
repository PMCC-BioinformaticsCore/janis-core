


from typing import Any

from janis_core.ingestion.galaxy.model.workflow import Workflow
from janis_core.ingestion.galaxy.model.workflow import WorkflowMetadata


def ingest_metadata(janis: Workflow, galaxy: dict[str, Any]) -> None:
    metadata = parse_metadata(galaxy)
    janis.set_metadata(metadata)

def parse_metadata(galaxy: dict[str, Any]) -> WorkflowMetadata:
    """scrape basic workflow metadata"""
    return WorkflowMetadata(
        name=galaxy['name'],
        uuid=galaxy['uuid'],
        annotation=galaxy['annotation'],
        version=galaxy['version'],
        tags=galaxy['tags']
    )