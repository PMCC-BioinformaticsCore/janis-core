

from __future__ import annotations
from typing import TYPE_CHECKING
from janis_core.ingestion.galaxy.fileio import safe_init_folder

if TYPE_CHECKING:
    from janis_core.ingestion.galaxy.model.tool import Tool
    from janis_core.ingestion.galaxy.model.workflow import Workflow
    from janis_core.ingestion.galaxy.model.workflow import WorkflowStep

import shutil
from janis_core.ingestion.galaxy.runtime import paths

from janis_core.ingestion.galaxy.utils import galaxy as galaxy_utils
from janis_core.ingestion.galaxy.gx.wrappers import fetch_xml

from .text.workflow.InputsText import InputsText
from .text.workflow.WorkflowText import WorkflowText
from .text.tool.ConfigfileText import ConfigfileText
from .text.tool.UnstranslatedText import UntranslatedText
from .text.tool.ToolText import ToolText




def write_tool(tool: Tool, path: str) -> None:
    text = ToolText(tool)
    page = text.render()
    with open(path, 'w') as fp:
        fp.write(page)

def write_workflow(janis: Workflow) -> None:
    write_untranslated(janis)
    write_scripts(janis)
    write_wrappers(janis)
    write_main_workflow(janis)
    write_inputs(janis)
    write_tools(janis)
    #write_sub_workflows(janis)
    #write_config(janis)

def write_tools(janis: Workflow) -> None:
    for step in janis.steps:
        tool_id = step.metadata.wrapper.tool_id
        write_tool(step.tool, paths.tool(tool_id))

def write_untranslated(janis: Workflow) -> None:
    for step in janis.steps:
        if step.preprocessing or step.postprocessing:
            tool_id = step.metadata.wrapper.tool_id
            path = paths.untranslated(tool_id)
            text = UntranslatedText(step)
            page = text.render()
            with open(path, 'w') as fp:
                fp.write(page)

def write_scripts(janis: Workflow) -> None:
    for step in janis.steps:
        if step.tool.configfiles:
            tool_id = step.metadata.wrapper.tool_id
            for configfile in step.tool.configfiles:
                path = paths.configfile(tool_id, configfile.name)
                text = ConfigfileText(configfile)
                page = text.render()
                with open(path, 'w') as fp:
                    fp.write(page)

def write_wrappers(janis: Workflow) -> None:
    for step in janis.steps:
        src_files = get_wrapper_files_src(step)
        dest = get_dest_dir(step)
        safe_init_folder(dest)
        for src in src_files:
            shutil.copy2(src, dest)

def get_wrapper_files_src(step: WorkflowStep) -> list[str]:
    wrapper = step.metadata.wrapper
    wrapper_path = fetch_xml(
        owner= wrapper.owner,
        repo= wrapper.repo,
        revision= wrapper.revision,
        tool_id= wrapper.tool_id
    )
    wrapper_dir = wrapper_path.rsplit('/', 1)[0]
    macro_xmls = galaxy_utils.get_macros(wrapper_dir)
    return [wrapper_path] + macro_xmls

def get_dest_dir(step: WorkflowStep) -> str:
    tool_id = step.metadata.wrapper.tool_id
    revision = step.metadata.wrapper.revision
    return paths.wrapper(tool_id, revision)

def write_main_workflow(janis: Workflow) -> None:
    path = paths.workflow()
    text = WorkflowText(janis)
    page = text.render()
    with open(path, 'w') as fp:
        fp.write(page)

def write_inputs(janis: Workflow) -> None:
    path = paths.inputs(file_format='yaml')
    text = InputsText(janis, file_format='yaml')
    page = text.render()
    with open(path, 'w') as fp:
        fp.write(page)

def write_sub_workflows(janis: Workflow) -> None:
    raise NotImplementedError()

def write_config(janis: Workflow) -> None:
    raise NotImplementedError()




