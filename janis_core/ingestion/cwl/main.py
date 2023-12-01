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
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory

from .identifiers import get_id_filename
from .identifiers import get_id_entity
from .loading import load_cwl_version
from .loading import load_cwl_utils_from_version
from .loading import load_cwl_document
from .loading import convert_etool_to_cltool

from .graph import add_step_edges_to_graph

from .parsing.tool import CLTParser
from .parsing.workflow import WorkflowInputParser
from .parsing.workflow import WorkflowOutputParser
from .parsing.workflow import WorkflowStepInputsParser
from .parsing.workflow import WorkflowStepAttributesParser



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
        self.cwl_utils = load_cwl_utils_from_version(self.version)

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
        # TODO: currently no error handling here.
        identifier = get_id_filename(workflow.id)

        wf = j.WorkflowBuilder(
            identifier=identifier,
            friendly_name=workflow.label,
            doc=workflow.doc,
        )

        for inp in workflow.inputs:
            self.ingest_workflow_input(wf, inp)

        # first step ingest pass
        for cwl_step in workflow.steps:
            j_step = self.ingest_workflow_step(wf, cwl_step)
        
        # second step ingest pass
        for cwl_step in workflow.steps:
            self.ingest_workflow_step_attributes(wf, j_step, cwl_step)
        
        # third step ingest pass
        for cwl_step in workflow.steps:
            self.ingest_workflow_step_inputs(wf, j_step, cwl_step)
        
        for out in workflow.outputs:
            self.ingest_workflow_output(wf, out)

        return wf
    
    def ingest_workflow_input(self, wf: j.Workflow, inp: Any) -> j.InputNodeSelector:
        parser = WorkflowInputParser(cwl_utils=self.cwl_utils, entity=inp, wf=wf, entity_uuid=wf.uuid)
        return parser.parse()

    def ingest_workflow_output(self, wf: j.Workflow, out: Any) -> OutputNode:
        parser = WorkflowOutputParser(cwl_utils=self.cwl_utils, entity=out, wf=wf, entity_uuid=wf.uuid)
        return parser.parse()

    def ingest_workflow_step(self, wf: j.Workflow, cwlstp: Any) -> StepNode:
        step_identifier = get_id_entity(cwlstp.id)
        
        if isinstance(cwlstp.run, (self.cwl_utils.CommandLineTool, self.cwl_utils.Workflow)):
            tool = self.ingest(cwlstp.run)
        else:
            tool = parse(cwlstp.run, os.path.dirname(self.doc))

        # if _foreach is not None:
        #     wf.has_scatter = True

        return wf.step(
            identifier=step_identifier,
            tool=tool,
            doc=cwlstp.doc,
            ignore_missing=True
        ) 
        
    def ingest_workflow_step_attributes(self, wf: j.Workflow, j_step: Any, cwlstp: Any) -> None:
        parser = WorkflowStepAttributesParser(cwl_utils=self.cwl_utils, entity=cwlstp, wf=wf, entity_uuid=j_step.uuid)
        parser.parse()

        if parser.scatter is not None:
            jstep.scatter = parser.scatter
            wf.has_scatter = True
        
        if parser.when is not None:
            jstep.when = parser.when

    def ingest_workflow_step_inputs(self, wf: j.Workflow, j_step: Any, cwlstp: Any) -> None:
        parser = WorkflowStepInputsParser(cwl_utils=self.cwl_utils, entity=cwlstp, wf=wf, entity_uuid=j_step.uuid)
        inputs_dict = parser.parse()
        jstep.tool.connections = inputs_dict
        add_step_edges_to_graph(jstep, wf)
        self.ingest_workflow_step_inputs_pickvalue(jstep, cwlstp)

    def ingest_workflow_step_inputs_pickvalue(self, jstep: StepNode, cwlstp: Any) -> None:
        operator_map = {
            'first_non_null': j.FirstOperator,
            'the_only_non_null': j.FirstOperator,
            'all_non_null': j.FilterNullOperator,
        }

        # new: multiple input sources with selection method
        for inp in cwlstp.in_:
            inp_identifier = get_id_entity(inp.id)
            if inp_identifier not in jstep.sources:
                continue 
            
            # must have more than one source
            sti = jstep.sources[inp_identifier]
            if len(sti.source_map) <= 1:
                continue 

            # must have pickValue set on cwl step input
            if not hasattr(inp, 'pickValue') or inp.pickValue is None:
                continue 

            operator = operator_map[inp.pickValue]
            sti.operator = operator(sti.source_map)

    def ingest_command_line_tool(self, clt: Any, is_expression_tool: bool=False):
        parser = CLTParser(cwl_utils=self.cwl_utils, clt=clt, entity=clt, is_expression_tool=is_expression_tool)
        return parser.parse()
            
    def ingest_expression_tool(self, etool: Any) -> j.CommandTool:
        # j.Logger.warn(
        #     f"Expression tools aren't well converted to Janis as they rely on unimplemented functionality: {clt.id}"
        # )
        # cast to CommandLineTool then parse as CommandLineTool
        clt = self.parse_etool_to_cltool(etool)
        tool = self.ingest_command_line_tool(clt, is_expression_tool=True)
        msg = 'Translation of CWL ExpressionTools is currently an experimental feature of janis translate'
        log_message(tool.uuid, msg, ErrorCategory.EXPERIMENTAL)
        return tool
    
    def parse_etool_to_cltool(self, etool: Any) -> Any:
        clt = convert_etool_to_cltool(etool, self.version)
        for out in clt.outputs:
            out_id = get_id_entity(out.id)
            out.janis_collection_override = j.ReadJsonOperator('cwl.output.json')[f"'{out_id}'"]
        
        # change 'expression.js' in base_command to something more meaningful 
        # base_command: ['nodejs', 'expression.js'] -> ['nodejs', '{toolname}.js']
        # update the InitialWorkDirRequirement for expression.js to be {toolname}.js
        toolname = get_id_filename(clt.id)
        clt.baseCommand = ['nodejs', f'{toolname}.js']
        for req in clt.requirements:
            if isinstance(req, self.cwl_utils.InitialWorkDirRequirement):
                req.listing[0].entryname = f'{toolname}.js'
        
        return clt

