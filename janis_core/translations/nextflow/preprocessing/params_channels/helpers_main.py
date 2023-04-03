

from janis_core.workflow.workflow import Workflow

from .helpers_common import get_file_wf_inputs
from .helpers_common import get_scatter_wf_inputs
from .helpers_common import get_all_workflow_inputs

from ...scope import Scope


def get_param_inputs_to_register_main(wf: Workflow, scope: Scope) -> set[str]:
    return get_all_workflow_inputs(wf)

def get_channel_inputs_to_register_main(wf: Workflow, scope: Scope) -> set[str]:
    """
    Get the wf inputs for which we will create a nf channel in the main wf.
    inputs which are:
    - files
    - filenames (in some cases)
    - scattered on 
    will become nextflow channels. all other inputs will become global scope params. 
    """
    all_inputs = get_all_workflow_inputs(wf)
    file_inputs = all_inputs & get_file_wf_inputs(wf)
    scatter_inputs = all_inputs & get_scatter_wf_inputs(wf)

    final_inputs = file_inputs | scatter_inputs
    return final_inputs

