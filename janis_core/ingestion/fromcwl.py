#!/usr/bin/env python3

import re
import os
from typing import Optional, Union, List

from janis_core.utils.validators import Validators

import janis_core as j


DEFAULT_PARSER_VERSION = "v1.2"


class CWlParser:

    parsed_cache = {}
    file_datatype_cache = {}

    def __init__(self, cwl_version: str, base_uri: str = None):
        self.cwl_version = cwl_version
        self.base_uri = base_uri
        self.cwlgen, self.cwlgen_etool_to_cltool = self.load_cwlgen_from_version(
            cwl_version=cwl_version
        )

    def get_data_type_from_secondaries(cls, secondaries: List[str], optional: bool):
        def calcluate_hash_of_set(s):
            return hash("|".join(sorted(set(s))))

        if not cls.file_datatype_cache:

            FastaGzType = None
            try:
                from janis_bioinformatics.data_types import FastaGz

                FastaGzType = FastaGz
            except ImportError:
                pass
            dts = j.JanisShed.get_all_datatypes()
            cls.file_datatype_cache = {
                calcluate_hash_of_set(dt.secondary_files()): dt
                for dt in dts
                if issubclass(dt, j.File)
                and dt().secondary_files()
                and (FastaGzType is not None and not issubclass(dt, FastaGzType))
            }

        sec_hash = calcluate_hash_of_set(secondaries)
        if sec_hash in cls.file_datatype_cache:
            return cls.file_datatype_cache[sec_hash](optional=optional)

        return j.GenericFileWithSecondaries(secondaries=secondaries)

    @staticmethod
    def from_doc(doc: str, base_uri=None):
        abs_path = os.path.relpath(doc)

        if abs_path in CWlParser.parsed_cache:
            return CWlParser.parsed_cache[abs_path]

        initial = None
        if base_uri:
            if base_uri.startswith("file://"):
                base_uri = base_uri[6:]
            initial = os.getcwd()
            os.chdir(base_uri)
        cwl_version = CWlParser.load_cwl_version_from_doc(doc)
        parser = CWlParser(cwl_version=cwl_version, base_uri=os.path.dirname(doc))

        tool = parser.from_document(doc)
        CWlParser.parsed_cache[abs_path] = tool
        if initial:
            os.chdir(initial)
        return tool

    def from_document(self, doc):

        loaded_doc = self.cwlgen.load_document(doc)
        return self.from_loaded_doc(loaded_doc)

    def from_loaded_doc(self, loaded_doc) -> j.Tool:

        if isinstance(loaded_doc, self.cwlgen.CommandLineTool):
            return self.ingest_command_line_tool(loaded_doc)

        elif isinstance(loaded_doc, self.cwlgen.Workflow):
            return self.ingest_workflow(loaded_doc)

        elif isinstance(loaded_doc, self.cwlgen.ExpressionTool):
            return self.ingest_expression_tool(loaded_doc)
        else:
            raise Exception(
                f"Janis doesn't support ingesting from {type(loaded_doc).__name__}"
            )

    def ingest_cwl_type(self, cwl_type, secondary_files):
        inp_type = self.from_cwl_inner_type(cwl_type)
        if secondary_files:
            array_optional_layers = []
            while isinstance(inp_type, j.Array):
                array_optional_layers.append(inp_type.optional)
                inp_type = inp_type.subtype()

            inp_type = self.get_data_type_from_secondaries(
                secondaries=self.process_secondary_files(secondary_files),
                optional=inp_type.optional,
            )
            for is_optional in array_optional_layers[::-1]:
                inp_type = j.Array(inp_type, optional=is_optional)

        return inp_type

    def from_cwl_inner_type(self, cwl_type) -> j.DataType:
        if isinstance(cwl_type, str):
            optional = "?" in cwl_type
            cwl_type = cwl_type.replace("?", "")
            array_count = 0
            while cwl_type.endswith("[]"):
                array_count += 1
                cwl_type = cwl_type[:-2]

            if cwl_type == "File":
                inner = j.File
            elif cwl_type == "Directory":
                inner = j.Directory
            elif cwl_type == "string":
                inner = j.String
            elif cwl_type == "int":
                inner = j.Int
            elif cwl_type == "float":
                inner = j.Float
            elif cwl_type == "boolean":
                inner = j.Boolean
            elif cwl_type == "stdout":
                inner = j.Stdout
            elif cwl_type == "stderr":
                inner = j.Stderr
            elif cwl_type == "Any":
                inner = j.String
            elif cwl_type == "long":
                inner = j.Int
            else:
                raise Exception(f"Can't detect type {cwl_type}")
            return inner(optional=optional)

        elif isinstance(cwl_type, list):
            optional = None
            types = []
            for c in cwl_type:
                if c == "null":
                    optional = True
                else:
                    types.append(self.ingest_cwl_type(c, []))

            if len(types) == 1:
                if optional is not None:
                    types[0].optional = optional
                return types[0]
            else:
                from janis_core.types.common_data_types import UnionType

                if optional is not None:
                    for inner in types:
                        inner.optional = optional

                return UnionType(*types)

        elif isinstance(cwl_type, self.cwlgen.CommandInputArraySchema):
            return j.Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwlgen.InputArraySchema):
            return j.Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwlgen.CommandOutputArraySchema):
            return j.Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwlgen.OutputArraySchema):
            return j.Array(self.from_cwl_inner_type(cwl_type.items))
        elif isinstance(cwl_type, self.cwlgen.InputEnumSchema):
            return j.String()

        else:
            raise Exception(f"Can't parse type {type(cwl_type).__name__}")

    @classmethod
    def get_tool_tag_from_identifier(cls, identifier):
        i = cls.get_tag_from_identifier(identifier)

        while not Validators.validate_identifier(i):
            i = str(
                input(
                    f"The tag for tool: '{i}' (fullID: {identifier}) was invalid, please choose another: "
                )
            )

        return i

    @classmethod
    def get_tag_from_identifier(cls, identifier: any):
        identifier = cls.get_source_from_identifier(identifier)
        if "/" in identifier:
            identifier = str(identifier.split("/")[-1])
        while "-" in identifier:
            identifier = identifier.replace("-", "_")

        if identifier == "input":
            return "inp"
        if identifier == "output":
            return "outp"

        return identifier

    @classmethod
    def get_source_from_identifier(cls, identifier):
        if not isinstance(identifier, str):
            identifier = str(identifier)
        if "#" in identifier:
            identifier = str(identifier.split("#")[-1])

        while "-" in identifier:
            identifier = identifier.replace("-", "_")

        if identifier == "input":
            return "inp"
        if identifier == "output":
            return "outp"

        return identifier

    def process_secondary_files(self, secondary_files: List):
        if not hasattr(self.cwlgen, "SecondaryFileSchema"):
            return secondary_files

        return [
            s.pattern if isinstance(s, self.cwlgen.SecondaryFileSchema) else s
            for s in secondary_files
        ]

    function_token_matcher = re.compile(
        "^\$\{\s*?return\s+?(.+?);*?\s*?\}$"
    )  # ${ return "arriba" }
    single_token_matcher = re.compile("^\$\((.+)\)$")  # "$(some value here)"
    inline_expression_matcher = re.compile(
        "\$\((.+?)\)"
    )  # valueFrom: "Hello, my name is $(name)
    input_selector_matcher = re.compile(
        "^inputs\.([A-z0-9_]+)$"
    )  # valueFrom: $(inputs.myName) -> InputSelector("myName")
    string_matcher = re.compile('^".+?"$')  # literal strings

    @staticmethod
    def parse_number_from_string(num):
        """
        Parse a string that is expected to contain a number.
        :param num: str. the number in string.
        :return: float or int. Parsed num.
        """
        if not isinstance(num, str):  # optional - check type
            raise TypeError("num should be a str. Got {}.".format(type(num)))
        if re.compile("^\s*\d+\s*$").search(num):
            return int(num)
        if re.compile("^\s*(\d*\.\d+)|(\d+\.\d*)\s*$").search(num):
            return float(num)
        raise ValueError("num is not a number. Got {}.".format(num))  # optional

    def parse_basic_expression(self, expr):
        if expr is None:
            return None
        if not isinstance(expr, str):
            return expr
        match = self.single_token_matcher.match(expr)
        if match:
            return self.convert_javascript_token(match.groups()[0])

        bigger_match = self.function_token_matcher.match(expr)
        if bigger_match:
            return self.convert_javascript_token(bigger_match.groups()[0])

        tokens = set(self.inline_expression_matcher.findall(expr))

        string_format = f"{expr}"
        token_replacers = {}

        for token, idx in zip(tokens, range(len(tokens))):
            key = f"JANIS_CWL_TOKEN_{idx+1}"
            string_format = string_format.replace(f"$({token})", f"{{{key}}}")
            token_replacers[key] = self.convert_javascript_token(token)

        if len(token_replacers) == 0:
            return string_format

        return j.StringFormatter(string_format, **token_replacers)

    def convert_javascript_token(self, token: str):
        input_selector_match = self.input_selector_matcher.match(token)
        if input_selector_match:
            return j.InputSelector(input_selector_match.groups()[0])

        if token.endswith(".size"):
            return j.FileSizeOperator(self.convert_javascript_token(token[:-5]))
        if token.endswith(".basename"):
            return j.BasenameOperator(self.convert_javascript_token(token[:-9]))
        if token.endswith(".path"):
            # Ignore it because Janis will automatically put this back in where relevant
            return self.convert_javascript_token(token[:-5])
        if token.endswith(".contents"):
            return j.ReadContents(self.convert_javascript_token(token[:-9]))

        is_string = self.string_matcher.match(token)
        if is_string:
            return token[1:-1]

        try:
            return self.parse_number_from_string(token)
        except ValueError:
            pass

        j.Logger.warn(
            f"Couldn't translate javascript token, will use the placeholder '<expr>{token}</expr>'"
        )
        return f"<expr>{token}</expr>"

    def ingest_command_tool_argument(self, arg):

        if isinstance(arg, str):
            return j.ToolArgument(self.parse_basic_expression(arg))
        else:
            return j.ToolArgument(
                value=self.parse_basic_expression(arg.valueFrom),
                position=arg.position,
                prefix=arg.prefix,
                separate_value_from_prefix=arg.separate,
                shell_quote=arg.shellQuote,
            )

    def ingest_command_tool_input(self, inp):
        inpBinding = inp.inputBinding

        if inpBinding and inpBinding.valueFrom:
            j.Logger.warn(
                f"Won't translate the expression for input {inp.id}: {inpBinding.valueFrom}"
            )

        inp_type = self.ingest_cwl_type(inp.type, secondary_files=inp.secondaryFiles)

        return j.ToolInput(
            tag=self.get_tag_from_identifier(inp.id),
            input_type=inp_type,
            position=inpBinding.position if inpBinding else None,
            prefix=inpBinding.prefix if inpBinding else None,
            separate_value_from_prefix=inpBinding.separate if inpBinding else None,
            separator=inpBinding.itemSeparator if inpBinding else None,
            shell_quote=inpBinding.shellQuote if inpBinding else None,
            default=inp.default,
        )

    def ingest_expression_tool_input(self, inp):
        inp_type = self.ingest_cwl_type(inp.type, secondary_files=inp.secondaryFiles)
        return j.ToolInput(
            tag=self.get_tag_from_identifier(inp.id), input_type=inp_type
        )

    def ingest_expression_tool_output(self, out):
        out_type = self.ingest_cwl_type(out.type, secondary_files=out.secondaryFiles)

        return j.ToolOutput(
            self.get_tag_from_identifier(out.id),
            output_type=out_type,
            skip_output_quality_check=True,
        )

    def ingest_command_tool_output(
        self, out
    ):  # out: self.cwlgen.CommandOutputParameter
        outBinding = out.outputBinding

        selector = None
        if outBinding:
            if outBinding.glob:
                selector = j.WildcardSelector(
                    self.parse_basic_expression(outBinding.glob)
                )
            elif outBinding.outputEval:
                selector = self.parse_basic_expression(outBinding.outputEval)
        elif out.outputEval:
            selector = self.parse_basic_expression(out.outputEval)

        return j.ToolOutput(
            tag=self.get_tag_from_identifier(out.id),
            output_type=self.ingest_cwl_type(
                out.type, secondary_files=out.secondaryFiles
            ),
            selector=selector,
        )

    def parse_workflow_source(
        self, wf: j.Workflow, step_input, potential_prefix: Optional[str] = None
    ):
        if step_input is None:
            return None
        if isinstance(step_input, list):
            return [
                self.parse_workflow_source(wf, si, potential_prefix=potential_prefix)
                for si in step_input
            ]

        if not isinstance(step_input, str):
            raise Exception(f"Can't parse step_input {step_input}")

        parsed_step_input = self.get_source_from_identifier(step_input)
        if parsed_step_input.startswith(wf.id() + "/"):
            parsed_step_input = parsed_step_input[len(wf.id()) + 1 :]
        if potential_prefix and parsed_step_input.startswith(potential_prefix + "/"):
            parsed_step_input = parsed_step_input[len(potential_prefix) + 1 :]

        if parsed_step_input.startswith("$("):
            raise Exception(
                f"This script can't parse expressions in the step input {step_input}"
            )

        [*ignore, source_str, tag_str] = (
            parsed_step_input.split("/")
            if "/" in parsed_step_input
            else (parsed_step_input, None)
        )

        tag_str = self.get_tag_from_identifier(tag_str)
        if source_str not in wf.nodes:
            raise KeyError(f"Couldn't find input / step {source_str} in nodes")
        source = wf[source_str]
        from janis_core.workflow.workflow import StepNode

        if tag_str and isinstance(source, StepNode):
            source = source.get_item(tag_str)
        return source

    def ingest_workflow_input(self, wf: j.Workflow, inp):

        return wf.input(
            identifier=self.get_tag_from_identifier(inp.id),
            datatype=self.ingest_cwl_type(inp.type, secondary_files=inp.secondaryFiles),
            default=inp.default,
            doc=inp.doc,
        )

    def ingest_workflow_output(self, wf: j.Workflow, out):
        import cwl_utils.parser_v1_2 as cwlgen

        out: cwlgen.WorkflowOutputParameter = out
        identifier = self.get_tag_from_identifier(out.id)
        out_source = self.parse_workflow_source(
            wf, out.outputSource, potential_prefix=identifier
        )
        return wf.output(
            identifier=identifier,
            datatype=self.ingest_cwl_type(out.type, secondary_files=out.secondaryFiles),
            source=out_source,
        )

    def ingest_workflow_step(self, wf: j.Workflow, stp):
        import cwl_utils.parser_v1_2 as cwlgen

        stp: cwlgen.WorkflowStep = stp
        step_identifier = self.get_tag_from_identifier(stp.id)

        if isinstance(stp.run, (self.cwlgen.CommandLineTool, self.cwlgen.Workflow)):
            tool = self.from_loaded_doc(stp.run)
        else:
            tool = CWlParser.from_doc(stp.run, base_uri=self.base_uri)

        inputs = {}
        for inp in stp.in_:
            inp: cwlgen.WorkflowStepInput = inp
            inp_identifier = self.get_tag_from_identifier(inp.id)

            source = None
            if inp.source is not None:
                source = self.parse_workflow_source(
                    wf, inp.source, potential_prefix=step_identifier
                )
            elif inp.valueFrom is not None:
                source = self.parse_basic_expression(inp.valueFrom)
            elif inp.default:
                source = inp.default

            if source is None:
                print(f"Source is None from object: {inp.save()}")
            inputs[inp_identifier] = source

        scatter = None
        if stp.scatter:
            scatter_fields_raw = stp.scatter
            if not isinstance(scatter_fields_raw, list):
                scatter_fields_raw = [scatter_fields_raw]

            scatter_fields = []
            for field in scatter_fields_raw:
                [*other_fields, input_to_scatter] = field.split("/")
                scatter_fields.append(input_to_scatter)

            scatter_method = stp.scatterMethod
            scatter = j.ScatterDescription(
                fields=scatter_fields, method=self.ingest_scatter_method(scatter_method)
            )

        return wf.step(
            identifier=step_identifier,
            tool=tool(**inputs),
            scatter=scatter,
            when=None,
            doc=stp.doc,
        )

    def ingest_scatter_method(self, scatter_method) -> j.ScatterMethod:
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

        raise Exception(f"Unrecognised scatter method '{scatter_method}'")

    def ingest_command_line_tool(self, clt):

        docker_requirement = None  # : Optional[self.cwlgen.DockerRequirement]
        files_to_create = {}
        memory, cpus, time = None, None, None
        for req in clt.requirements or []:
            if isinstance(req, self.cwlgen.DockerRequirement):
                docker_requirement = req

            elif isinstance(req, self.cwlgen.InitialWorkDirRequirement):
                for dirent in req.listing:
                    if isinstance(dirent, str):
                        continue
                    files_to_create[
                        self.parse_basic_expression(dirent.entryname)
                    ] = dirent.entry

            elif isinstance(req, self.cwlgen.ResourceRequirement):
                # maybe convert mebibytes to megabytes?
                memory = self.parse_basic_expression(req.ramMin or req.ramMax)
                cpus = self.parse_basic_expression(req.coresMin)
            elif (
                hasattr(self.cwlgen, "ToolTimeLimit")
                and isinstance(req, self.cwlgen.ToolTimeLimit)
                and req.timelimit > 0
            ):
                time = req.timelimit

        container = None
        if docker_requirement:
            container = docker_requirement.dockerPull

        tool_id = self.get_tool_tag_from_identifier(clt.id)
        jclt = j.CommandToolBuilder(
            tool=tool_id,
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

    def ingest_workflow(self, workflow):

        wf = j.WorkflowBuilder(
            identifier=self.get_tag_from_identifier(workflow.id),
            friendly_name=workflow.label,
            doc=workflow.doc,
        )

        for inp in workflow.inputs:
            self.ingest_workflow_input(wf, inp)

        steps_to_process = workflow.steps
        iters = len(steps_to_process) ** 2
        while len(steps_to_process) > 0:
            stp = steps_to_process.pop(0)
            try:
                self.ingest_workflow_step(wf, stp)
            except KeyError as e:
                steps_to_process.append(stp)
                if iters < 0:
                    print("")
            iters -= 1
            if iters < 0:
                print("Oh god...")
            if iters < -len(steps_to_process) ** 2:
                raise Exception("Probably stuck in some weird loop, we're ")
        # for stp in workflow.steps:
        #     self.ingest_workflow_step(wf, stp)

        for out in workflow.outputs:
            self.ingest_workflow_output(wf, out)

        return wf

    @classmethod
    def load_cwlgen_from_version(cls, cwl_version: str):

        if cwl_version == "v1.0":
            import cwl_utils.parser_v1_0 as cwlutils
            from cwl_utils.cwl_v1_0_expression_refactor import etool_to_cltool
        elif cwl_version == "v1.1":
            import cwl_utils.parser_v1_1 as cwlutils
            from cwl_utils.cwl_v1_0_expression_refactor import etool_to_cltool
        elif cwl_version == "v1.2":
            import cwl_utils.parser_v1_2 as cwlutils
            from cwl_utils.cwl_v1_2_expression_refactor import etool_to_cltool
        else:
            print(
                f"Didn't recognise CWL version {cwl_version}, loading default: {DEFAULT_PARSER_VERSION}"
            )
            cwlutils, etool_to_cltool = cls.load_cwlgen_from_version(
                DEFAULT_PARSER_VERSION
            )

        return cwlutils, etool_to_cltool

    @classmethod
    def load_cwl_version_from_doc(cls, doc: str) -> str:
        import ruamel.yaml

        # load tool into memory
        if doc.startswith("file://"):
            doc = doc[7:]
        with open(doc) as fp:
            tool_dict = ruamel.yaml.load(fp, Loader=ruamel.yaml.Loader)

        if "cwlVersion" not in tool_dict:
            raise Exception(f"Couldn't find cwlVersion in tool {doc}")

        return tool_dict["cwlVersion"]

    def ingest_expression_tool(self, expr_tool):
        # https://github.com/common-workflow-language/cwl-utils/pull/5
        clt = self.cwlgen_etool_to_cltool(expr_tool)
        j.Logger.warn(
            f"Expression tools aren't well converted to Janis as they rely on unimplemented functionality: {clt.id}"
        )
        for out in clt.outputs:
            out_id = out.id.split(".")[-1]
            out.outputEval = f"JANIS (potentially unimplemented): j.ReadJsonOperator(j.Stdout)[out_id]"

        return self.ingest_command_line_tool(clt)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        raise Exception("Expected 1 argument, the name of a CWL tool.")
    toolname = sys.argv[1]

    tool = CWlParser.from_doc(toolname)

    tool.translate("janis")
