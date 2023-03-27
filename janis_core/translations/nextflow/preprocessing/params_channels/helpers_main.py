

from janis_core.workflow.workflow import Workflow

from .helpers_common import get_file_wf_inputs
from .helpers_common import get_filename_wf_inputs
from .helpers_common import get_scatter_wf_inputs
from .helpers_common import get_true_workflow_inputs

from ...scope import Scope


def get_param_inputs_to_register_main(wf: Workflow, scope: Scope) -> set[str]:
    final_inputs = get_true_workflow_inputs(wf)
    return final_inputs


def get_channel_inputs_to_register_main(wf: Workflow, scope: Scope) -> set[str]:
    """
    Get the wf inputs for which we will create a nf channel in the main wf.
    inputs which are:
    - files
    - filenames (in some cases)
    - scattered on 
    will become nextflow channels. all other inputs will become global scope params. 
    """
    file_inputs = get_file_wf_inputs(wf)
    filename_inputs = get_filename_wf_inputs(wf)
    scatter_inputs = get_scatter_wf_inputs(wf)
    # subworkflow_inputs = _get_subworkflow_inputs(wf)
    # null_value_to_subworkflow_inputs = _get_null_value_to_subworkflow_inputs(wf)

    final_inputs = file_inputs | filename_inputs | scatter_inputs
    return final_inputs

