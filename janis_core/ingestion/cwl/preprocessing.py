

from typing import Any


def handle_inline_cltool_identifiers(cwl_workflow: Any) -> Any:
    """
    alters tool ids in a cwl_workflow. 
    workflows with inline CommandLineTool definitions (rather than those in a dedicated tool.cwl file)
    will have random ids starting with _: for those inline tool definitions. 
    in these cases, this function overrides the random id, and replaces it with an id based on the step id.
    """
    for step in cwl_workflow.steps:
        tool = step.run
        if tool.id.startswith('_:'):
            step_name = step.id.rsplit('#', 1)[1]
            workdir = step.id.rsplit('#', 1)[0].rsplit('/', 1)[0]
            new_id = f'{workdir}/{step_name}_tool.cwl'
            tool.id = new_id
    return cwl_workflow