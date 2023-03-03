#!/usr/bin/env python3

import copy

from typing import Any, Optional
import janis_core as j


DEFAULT_PARSER_VERSION = "v1.2"
from cwl_utils.parser.cwl_v1_0 import CommandLineTool as CommandLineTool_1_0
from cwl_utils.parser.cwl_v1_1 import CommandLineTool as CommandLineTool_1_1
from cwl_utils.parser.cwl_v1_2 import CommandLineTool as CommandLineTool_1_2
CommandLineTool = CommandLineTool_1_0 | CommandLineTool_1_1 | CommandLineTool_1_2

from janis_core import StepOutputSelector 
from janis_core.utils import first_value

from janis_core.workflow.workflow import InputNode, StepNode, OutputNode
from janis_core.workflow.workflow import verify_or_try_get_source
from janis_core.types.data_types import is_python_primitive
from janis_core.types import Filename
from janis_core.tool.documentation import (
    InputDocumentation,
    InputQualityType,
)

from janis_core.utils.errors import UnsupportedError

from .identifiers import get_id_filename
from .identifiers import get_id_path
from .identifiers import get_id_entity
from .identifiers import remove_output_name_from_output_source
from .preprocessing import handle_inline_cltool_identifiers
from .loading import load_cwl_version
from .loading import load_cwlgen_from_version
from .loading import load_cwl_document
from .loading import convert_etool_to_cltool
from .types import ingest_cwl_type
from .types import cast_cwl_type_to_python
from .parsing import parse_basic_expression

import os



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
        identifier = get_id_entity(inp.id)
        return wf.input(
            identifier=identifier,
            datatype=ingest_cwl_type(inp.type, self.cwl_utils, secondary_files=inp.secondaryFiles),
            default=inp.default,
            doc=inp.doc,
        )

    def ingest_workflow_output(self, wf: j.Workflow, out: Any) -> OutputNode:
        out_identifier = get_id_entity(out.id)
        if isinstance(out.outputSource, list):
            cwl_source = [remove_output_name_from_output_source(x) for x in out.outputSource]
        else:
            cwl_source = remove_output_name_from_output_source(out.outputSource)
        
        sources = self.parse_sources(wf, cwl_source)
        if len(sources) == 1:
            source = sources[0]
        else:
            source = sources
        return wf.output(
            identifier=out_identifier,
            datatype=ingest_cwl_type(out.type, self.cwl_utils, secondary_files=out.secondaryFiles),
            source=source,
        )

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
        
    def ingest_workflow_step_scatter(self, wf: j.Workflow, cwlstp: Any) -> None:
        step_identifier = get_id_entity(cwlstp.id)
        jstep = wf.step_nodes[step_identifier]
        
        scatter = None
        if cwlstp.scatter:
            scatter_fields_raw = cwlstp.scatter
            if not isinstance(scatter_fields_raw, list):
                scatter_fields_raw = [scatter_fields_raw]

            scatter_fields = []
            for field in scatter_fields_raw:
                [*other_fields, input_to_scatter] = field.split("/")
                scatter_fields.append(input_to_scatter)

            scatter_method = cwlstp.scatterMethod
            scatter = j.ScatterDescription(
                fields=scatter_fields, method=self.ingest_scatter_method(scatter_method)
            )
        
        if scatter is not None:
            jstep.scatter = scatter
            wf.has_scatter = True

    def ingest_scatter_method(self, scatter_method: Optional[str]) -> Optional[j.ScatterMethod]:
        if scatter_method is None or scatter_method == "":
            return None
        elif scatter_method == "dotproduct":
            return j.ScatterMethod.dot
        elif scatter_method == "nested_crossproduct":
            j.Logger.warn(
                "Requesting nested_crossproduct, but Janis only supports flat_crossproduct. Will fallback to flat_crossproduct"
            )
            return j.ScatterMethod.cross
        elif scatter_method == "flat_crossproduct":
            return j.ScatterMethod.cross

        raise UnsupportedError(f"Unsupported scatter method '{scatter_method}'")

    def ingest_workflow_step_inputs(self, wf: j.Workflow, cwlstp: Any) -> None:
        step_identifier = get_id_entity(cwlstp.id)
        jstep = wf.step_nodes[step_identifier]

        valid_step_inputs = self.get_valid_step_inputs(cwlstp, jstep)

        connections = {}
        for inp in valid_step_inputs:
            inp_identifier = get_id_entity(inp.id)
            source = self.get_input_source(wf, inp)
            connections[inp_identifier] = source
        
        jstep.tool.connections = connections
        self.add_step_edges(jstep, wf)

    def get_valid_step_inputs(self, cwlstp: Any, jstep: StepNode) -> list[Any]:
        return [x for x in cwlstp.in_ if self.is_valid_step_input(x, jstep)]
    
    def is_valid_step_input(self, inp: Any, jstep: StepNode) -> bool:
        # ERROR HANDLING
        inp_identifier = get_id_entity(inp.id)
        if inp_identifier in jstep.tool.inputs_map():
            print(f'WARNING: [STEP: {jstep.id()}] [TOOL: {jstep.tool.id()}] - Tool has no input named "{inp_identifier}". Ignoring.')
            return True
        return False
    
    def get_input_source(self, wf: j.Workflow, inp: Any) -> Any:
        source = None
        if inp.source is not None:
            sources = self.parse_sources(wf, inp.source)
            if len(sources) == 1:
                source = sources[0]
            else:
                source = sources
        elif inp.valueFrom is not None:
            source = parse_basic_expression(inp.valueFrom)
        elif inp.default:
            source = cast_cwl_type_to_python(inp.default)

        if source is None:
            print(f"Source is None from object: {inp.save()}")
        
        return source

    def add_step_edges(self, jstep: StepNode, wf: j.Workflow) -> None:
        connections = jstep.tool.connections
        tinputs = jstep.tool.inputs_map()
        jstep.sources = {}
        
        added_edges = []
        for (k, v) in connections.items():
            
            # static values provided when creating step.
            # janis wants to create a workflow input for these.
            isfilename = isinstance(v, Filename)
            if is_python_primitive(v) or isfilename:
                inp_identifier = f"{jstep.id()}_{k}"
                referencedtype = copy.copy(tinputs[k].intype) if not isfilename else v

                referencedtype.optional = True

                indoc = tinputs[k].doc
                indoc.quality = InputQualityType.configuration

                v = wf.input(
                    inp_identifier,
                    referencedtype,
                    default=v.generated_filename() if isfilename else v,
                    doc=indoc,
                )

            if v is None:
                inp_identifier = f"{jstep.id()}_{k}"
                doc = copy.copy(InputDocumentation.try_parse_from(tinputs[k].doc))
                doc.quality = InputQualityType.configuration
                v = wf.input(inp_identifier, tinputs[k].intype, default=v, doc=doc)

            verifiedsource = verify_or_try_get_source(v)
            if isinstance(verifiedsource, list):
                for vv in verifiedsource:
                    added_edges.append(jstep._add_edge(k, vv))
            else:
                added_edges.append(jstep._add_edge(k, verifiedsource))

        for e in added_edges:
            si = e.finish.sources[e.ftag] if e.ftag else first_value(e.finish.sources)
            wf.has_multiple_inputs = wf.has_multiple_inputs or si.multiple_inputs
    
    def ingest_command_line_tool(self, clt: Any):
        docker_requirement = None  # : Optional[self.cwlgen.DockerRequirement]
        files_to_create = {}
        memory, cpus, time = None, None, None
        
        for req in clt.requirements or []:

            if isinstance(req, self.cwl_utils.DockerRequirement):
                docker_requirement = req

            elif isinstance(req, self.cwl_utils.InitialWorkDirRequirement):
                for dirent in req.listing:
                    if isinstance(dirent, str):
                        continue
                    files_to_create[
                        parse_basic_expression(dirent.entryname)
                    ] = dirent.entry

            elif isinstance(req, self.cwl_utils.ResourceRequirement):
                # maybe convert mebibytes to megabytes?
                memory = parse_basic_expression(req.ramMin or req.ramMax)
                cpus = parse_basic_expression(req.coresMin)

            elif hasattr(req, 'timelimit') and isinstance(req, self.cwl_utils.ToolTimeLimit):
                time = req.timelimit
            
            else:
                pass

        container = None
        if docker_requirement:
            container = docker_requirement.dockerPull

        identifier = get_id_filename(clt.id)

        jclt = j.CommandToolBuilder(
            tool=identifier,
            base_command=clt.baseCommand,
            inputs=[self.ingest_command_tool_input(inp) for inp in clt.inputs],
            outputs=[self.ingest_command_tool_output(out) for out in clt.outputs],
            arguments=[
                self.ingest_command_tool_argument(arg) for arg in (clt.arguments or [])
            ],
            version="v0.1.0",
            container=container or "ubuntu:latest",
            doc=clt.doc,
            friendly_name=clt.label,
            files_to_create=files_to_create or None,
            memory=memory,
            cpus=cpus,
            time=time,
        )
        return jclt

    def ingest_command_tool_argument(self, arg: Any):

        if isinstance(arg, str):
            return j.ToolArgument(parse_basic_expression(arg))
        else:
            return j.ToolArgument(
                value=parse_basic_expression(arg.valueFrom),
                position=arg.position,
                prefix=arg.prefix,
                separate_value_from_prefix=arg.separate,
                shell_quote=arg.shellQuote,
            )

    def ingest_command_tool_input(self, inp: Any):
        inpBinding = inp.inputBinding

        if inpBinding and inpBinding.valueFrom:
            j.Logger.warn(
                f"Won't translate the expression for input {inp.id}: {inpBinding.valueFrom}"
            )

        inp_type = ingest_cwl_type(inp.type, self.cwl_utils, secondary_files=inp.secondaryFiles)
        identifier = get_id_entity(inp.id)

        return j.ToolInput(
            tag=identifier,
            input_type=inp_type,
            position=inpBinding.position if inpBinding else None,
            prefix=inpBinding.prefix if inpBinding else None,
            separate_value_from_prefix=inpBinding.separate if inpBinding else None,
            separator=inpBinding.itemSeparator if inpBinding else None,
            shell_quote=inpBinding.shellQuote if inpBinding else None,
            default=inp.default,
        )
    
    def ingest_command_tool_output(self, out: Any, is_expression_tool: bool=False) -> j.ToolOutput:  
        # out: self.cwlgen.CommandOutputParameter
        
        # tag
        identifier = get_id_entity(out.id)
        
        # datatype
        output_type = ingest_cwl_type(out.type, self.cwl_utils, secondary_files=out.secondaryFiles)

        # selector
        selector = None
        if hasattr(out, 'janis_collection_override'):
            selector = out.janis_collection_override

        elif out.outputBinding:
            if out.outputBinding.glob:
                selector = j.WildcardSelector(
                    parse_basic_expression(out.outputBinding.glob)
                )
            elif out.outputBinding.outputEval:
                selector = parse_basic_expression(out.outputBinding.outputEval)

        return j.ToolOutput(identifier, output_type, selector)
    
    def ingest_expression_tool(self, expr_tool: Any) -> j.CommandToolBuilder:
        clt = convert_etool_to_cltool(expr_tool, self.version)
        j.Logger.warn(
            f"Expression tools aren't well converted to Janis as they rely on unimplemented functionality: {clt.id}"
        )
        for out in clt.outputs:
            out_id = get_id_entity(out.id)
            out.janis_collection_override = j.ReadJsonOperator(j.Stdout)[out_id]

        return self.ingest_command_line_tool(clt)

    def parse_sources(self, wf: j.Workflow, sources: str | list[str]) -> list[InputNode | StepOutputSelector]:
        """
        each source is a workflow input, step output, or complex expression.
        input nodes will all be parsed into janis wf at this stage.
        we can check if the source is an input on the janis wf, then if not, must be a step output.
        """
        out: list[InputNode | StepOutputSelector] = []
        
        if isinstance(sources, str):
            sources = [sources]

        for src in sources:
            # get the wfinp / step output identifier
            identifier = get_id_entity(src)

            # is complex expression?
            if identifier.startswith("$("):
                raise UnsupportedError(
                    f"This script can't parse expressions in the step input {step_input}"
                )
            
            # is step output?
            # if referencing step output, that step will have already been parsed into the janis wf
            if get_id_path(src) and get_id_path(src) in wf.step_nodes:
                stp_id = get_id_path(src)
                out_id = get_id_entity(src)
                stp = wf.step_nodes[stp_id]
                selector = stp.get_item(out_id)
                out.append(selector)
            
            # is wf input?
            elif identifier in wf.input_nodes:
                resolved_src = wf[identifier]
                out.append(resolved_src) 
            
            else:
                raise NotImplementedError
            
        return out

    
