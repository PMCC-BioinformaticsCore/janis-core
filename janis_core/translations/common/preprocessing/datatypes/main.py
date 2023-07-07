

"""
exists to standardise datatypes which have secondary files. 
currently if a source is secondary type, and dest is secondary type, and their secondary files don't match,
an error is thrown. 

These are the cases we fix in this file:
- source is defined secondary type, dest is GenericFileWithSecondaries(secondaries=[])
- source is GenericFileWithSecondaries(secondaries=[]), dest is defined secondary type
- source is GenericFileWithSecondaries, dest is GenericFileWithSecondaries, but their secondary files don't match. 
"""

from janis_core import WorkflowBuilder

def balance_mismatch_secondary_types(main_wf: WorkflowBuilder) -> WorkflowBuilder:
    source_dest_map = gather_mismatch_connections(main_wf)
    source_dest_map = balance_mismatch_connections(source_dest_map)
    wf = apply_balanced_connections(main_wf, source_dest_map)
    return wf

class Connection:
    pass

def gather_mismatch_connections(main_wf: WorkflowBuilder) -> list[Connection]:
    pass