#!/usr/bin/env python3

import os

from typing import Any, Optional
import janis_core as j


DEFAULT_PARSER_VERSION = "v1.2"
from cwl_utils.parser.cwl_v1_0 import CommandLineTool as CommandLineTool_1_0
from cwl_utils.parser.cwl_v1_1 import CommandLineTool as CommandLineTool_1_1
from cwl_utils.parser.cwl_v1_2 import CommandLineTool as CommandLineTool_1_2
CommandLineTool = CommandLineTool_1_0 | CommandLineTool_1_1 | CommandLineTool_1_2


from janis_core.workflow.workflow import StepNode, OutputNode
from janis_core.utils.errors import UnsupportedError
from janis_core.messages import log_warning

from .identifiers import get_id_filename
from .identifiers import get_id_entity
from .preprocessing import handle_inline_cltool_identifiers
from .preprocessing import convert_cwl_types_to_python
from .loading import load_cwl_version
from .loading import load_cwlgen_from_version
from .loading import load_cwl_document
from .loading import convert_etool_to_cltool

from .graph import add_step_edges_to_graph

from .parsing.tool import CLTToolParser
from .parsing.workflow import WorkflowInputParser
from .parsing.workflow import WorkflowOutputParser
from .parsing.workflow import WorkflowStepInputsParser
from .parsing.workflow import WorkflowStepScatterParser

from janis_core import settings



parsed_cache = {}

def parse(doc: str, base_uri: Optional[str]=None) -> j.Tool:
    # main entry point to ingest a cwl file    
    initial_wd = os.getcwd()
    if base_uri:
        _swap_directory(base_uri)

    parser = CWlParser(doc, base_uri)
    cwl_entity = load_cwl_document(parser.doc, parser.version)
    janis_entity = parser.ingest(cwl_entity)

    if base_uri:
        _revert_directory(initial_wd)

    return janis_entity

def _swap_directory(directory: str) -> None:
    if directory.startswith("file://"):
        directory = directory[6:]
    os.chdir(directory)

def _revert_directory(directory: str) -> None:
    os.chdir(directory)



class CWlParser:

    """main class to parse a cwl_utils Workflow | CommandLineTool | ExpressionTool"""

    def __init__(self, doc: str, base_uri: Optional[str]=None) -> None:
        self.doc = doc
        self.base_uri = base_uri
        self.version = load_cwl_version(doc)
        self.cwl_utils = load_cwlgen_from_version(self.version)

    def ingest(self, cwl_entity: Any) -> j.Tool:
        if isinstance(cwl_entity, self.cwl_utils.Workflow):
            return self.ingest_workflow(cwl_entity)
        if isinstance(cwl_entity, self.cwl_utils.CommandLineTool):
            return self.ingest_command_line_tool(cwl_entity)
        if isinstance(cwl_entity, self.cwl_utils.ExpressionTool):
            return self.ingest_expression_tool(cwl_entity)
        else:
            raise UnsupportedError(
                f"Janis doesn't support ingesting from {type(cwl_entity).__name__}"
            )
        
    def ingest_workflow(self, workflow: Any):
        workflow = convert_cwl_types_to_python(workflow, self.cwl_utils)
        workflow = handle_inline_cltool_identifiers(workflow)
        identifier = get_id_filename(workflow.id)

        wf = j.WorkflowBuilder(
            identifier=identifier,
            friendly_name=workflow.label,
            doc=workflow.doc,
        )

        for inp in workflow.inputs:
            self.ingest_workflow_input(wf, inp)

        # first step ingest pass
        for step in workflow.steps:
            self.ingest_workflow_step(wf, step)
        
        # second step ingest pass
        for step in workflow.steps:
            self.ingest_workflow_step_scatter(wf, step)
        
        # third step ingest pass
        for step in workflow.steps:
            self.ingest_workflow_step_inputs(wf, step)

        for out in workflow.outputs:
            self.ingest_workflow_output(wf, out)

        return wf
    
    def ingest_workflow_input(self, wf: j.Workflow, inp: Any) -> j.InputNodeSelector:
        parser = WorkflowInputParser(self.cwl_utils)
        return parser.parse(inp, wf)

    def ingest_workflow_output(self, wf: j.Workflow, out: Any) -> OutputNode:
        parser = WorkflowOutputParser(self.cwl_utils)
        return parser.parse(out, wf)

    def ingest_workflow_step(self, wf: j.Workflow, cwlstp: Any) -> StepNode:

        if isinstance(cwlstp.run, (self.cwl_utils.CommandLineTool, self.cwl_utils.Workflow)):
            tool = self.ingest(cwlstp.run)
        else:
            tool = parse(cwlstp.run, os.path.dirname(self.doc))

        # if _foreach is not None:
        #     wf.has_scatter = True

        step_identifier = get_id_entity(cwlstp.id)
        return wf.step(
            identifier=step_identifier,
            tool=tool,
            doc=cwlstp.doc,
            ignore_missing=True
        ) 
        
    def ingest_workflow_step_scatter(self, wf: j.Workflow, cwlstp: Any) -> None:
        step_identifier = get_id_entity(cwlstp.id)
        jstep = wf.step_nodes[step_identifier]

        parser = WorkflowStepScatterParser(self.cwl_utils)
        scatter = parser.parse(cwlstp, jstep)

        if scatter is not None:
            jstep.scatter = scatter
            wf.has_scatter = True

    def ingest_workflow_step_inputs(self, wf: j.Workflow, cwlstp: Any) -> None:
        parser = WorkflowStepInputsParser(self.cwl_utils)
        inputs_dict = parser.parse(cwlstp, wf)

        step_identifier = get_id_entity(cwlstp.id)
        jstep = wf.step_nodes[step_identifier]
        jstep.tool.connections = inputs_dict
        add_step_edges_to_graph(jstep, wf)

    def ingest_command_line_tool(self, clt: Any):
        if settings.ingest.SAFE_MODE:
            # some more basic way of getting the stuff we want? 
            pass
        else:
            parser = CLTToolParser(self.cwl_utils)
            return parser.parse(clt)

            
    def ingest_expression_tool(self, expr_tool: Any) -> j.CommandTool:
        # j.Logger.warn(
        #     f"Expression tools aren't well converted to Janis as they rely on unimplemented functionality: {clt.id}"
        # )

        clt = convert_etool_to_cltool(expr_tool, self.version)
        for out in clt.outputs:
            out_id = get_id_entity(out.id)
            out.janis_collection_override = j.ReadJsonOperator(j.Stdout)[out_id]

        tool = self.ingest_command_line_tool(clt)
        msg = 'Translation of CWL ExpressionTools is currently an experimental feature of janis translate'
        log_warning(tool.uuid, msg)
        return tool

