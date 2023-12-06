


# from typing import Optional

# from janis_core import Workflow, CodeTool, CommandToolBuilder
# from ...model.process import NFProcess
# from ...model.workflow import NFWorkflow
# from ...model.files import NFFile

# from .process import generate_file_process
# from .workflow import generate_file_workflow

# def generate_files(main_wf: Workflow, nf_processes: dict[str, NFProcess], nf_workflows: dict[str, NFWorkflow]) -> dict[str, NFFile]:
#     """generates nextflow files for processes and workflows"""
#     nf_files: dict[str, NFFile] = {}
    
#     for tool_id, process in nf_processes.items():
#         tool = _get_tool(tool_id, main_wf)
#         nffile = generate_file_process(process, tool)
#         nf_files[tool_id] = nffile
    
#     for wf_id, workflow in nf_workflows.items():
#         is_subworkflow = True if wf_id != main_wf.id() else False
#         wf = _get_workflow(wf_id, main_wf)
#         nffile = generate_file_workflow(workflow, nf_processes, nf_workflows, wf, is_subworkflow)
#         nf_files[wf_id] = nffile

#     return nf_files





