


from typing import Any, Optional
from janis_core.ingestion.galaxy.gxtool.command.cmdstr.DynamicCommandStatement import DynamicCommandStatement
from janis_core.ingestion.galaxy.gxtool.command.cmdstr.generate import gen_command_string
from janis_core.ingestion.galaxy.gxtool.command.cmdstr.CommandString import CommandStringSource
from janis_core.ingestion.galaxy.internal_model.workflow.workflow import Workflow

from janis_core.ingestion.galaxy.gxtool.parsing import load_xmltool
from janis_core.ingestion.galaxy.gxworkflow import load_tool_state
from janis_core.ingestion.galaxy.gxtool.command import load_templated_command_str

from janis_core.ingestion.galaxy import internal_mapping
from janis_core.ingestion.galaxy import runtime


def ingest_workflow_steps_prepost(janis: Workflow, galaxy: dict[str, Any]) -> None:
    for g_step in galaxy['steps'].values():
        if g_step['type'] == 'tool':
            ingest_prepost(g_step, janis, galaxy)

def ingest_prepost(g_step: dict[str, Any], janis: Workflow, galaxy: dict[str, Any]):
    # get janis step & update settings
    j_step = internal_mapping.step(g_step['id'], janis, galaxy)
    runtime.tool.set(from_wrapper=j_step.metadata.wrapper)

    xmltool = load_xmltool(runtime.tool.tool_path)
    tool_state = load_tool_state(
        g_step, 
        additional_filters=[
            'ReplaceNullWithVarname'
            'ReplaceConnectedWithVarname',
            'ReplaceRuntimeWithVarname',
        ]
    )
    command = load_templated_command_str(inputs_dict=tool_state)
    cmdstr = gen_command_string(source=CommandStringSource.XML, text=command, xmltool=xmltool)
    j_step.preprocessing = extract_cmdline(cmdstr.preprocessing)
    j_step.postprocessing = extract_cmdline(cmdstr.postprocessing)

def extract_cmdline(statements: list[DynamicCommandStatement]) -> Optional[str]:
    if not statements:
        return None
    cmdlines = [x.cmdline for x in statements]
    return ' &&'.join(cmdlines)

