

from janis_core.ingestion.galaxy import aliases # leave at top
from janis_core.ingestion.galaxy import settings
from typing import Optional

from janis_core.ingestion.galaxy.startup import general_setup
from janis_core.ingestion.galaxy.cli import CLIparser
from janis_core.ingestion.galaxy import paths

from janis_core.ingestion.galaxy.startup import tool_setup
from janis_core.ingestion.galaxy.startup import workflow_setup
from janis_core.ingestion.galaxy.ingest import ingest_tool
from janis_core.ingestion.galaxy.ingest import ingest_workflow

from janis_core.ingestion.galaxy.fileio import write_workflow
from janis_core.ingestion.galaxy.fileio import write_tool


import sys


"""
gxtool2janis program entry point
parses cli settings then hands execution to other files based on command
"""


def main():
    args = CLIparser(sys.argv).args
    general_setup(args)

    if args['command'] == 'tool':
        tool_mode(args)
    elif args['command'] == 'workflow':
        workflow_mode(args)


def tool_mode(args: dict[str, Optional[str]]) -> None:
    tool_setup(args)
    tool = ingest_tool(settings.tool.tool_path)
    path = paths.tool(tool.metadata.id)
    write_tool(tool, path=path) 

def workflow_mode(args: dict[str, Optional[str]]) -> None:
    workflow_setup(args)
    wf = ingest_workflow(settings.workflow.workflow_path)
    write_workflow(wf)


if __name__ == '__main__':
    main()
