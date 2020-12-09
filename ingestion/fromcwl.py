#!/usr/bin/env python3

import re
from typing import Optional, Union, List

import janis_core as j


DEFAULT_PARSER_VERSION = "v1.2"


class GenericFileWithSecondaries(j.File):
    def __init__(self, optional=False, secondaries: List[str] = None):
        super().__init__(optional=optional)
        self.secondaries = secondaries

    def secondary_files(self) -> Optional[List[str]]:
        return self.secondaries


class CWlParser:
    def __init__(self, cwl_version: str):
        self.cwlgen = self.load_cwlgen_from_version(cwl_version=cwl_version)

    @staticmethod
    def from_doc(doc):
        cwl_version = CWlParser.load_cwl_version_from_doc(doc)
        parser = CWlParser(cwl_version=cwl_version)
        return parser.from_document(doc)

    def from_document(self, doc):

        loaded_doc = self.cwlgen.load_document(doc)
        return self.from_loaded_doc(loaded_doc)

    def from_loaded_doc(self, loaded_doc) -> j.Tool:

        if isinstance(loaded_doc, self.cwlgen.CommandLineTool):
            return self.ingest_command_line_tool(loaded_doc)

        elif isinstance(loaded_doc, self.cwlgen.Workflow):
            raise Exception("Janis can't ingest from workflow yet")

        else:
            raise Exception(
                f"Janis doesn't support ingesting from {type(loaded_doc).__name__}"
            )

    def from_cwl_type(self, cwl_type):
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
                    types.append(self.from_cwl_type(c))

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
            return j.Array(self.from_cwl_type(cwl_type.items))

        else:
            raise Exception(f"Can't parse type {type(cwl_type).__name__}")

    @classmethod
    def get_tag_from_identifier(cls, identifier: any):
        if not isinstance(identifier, str):
            identifier = str(identifier)
        if "#" in identifier:
            identifier = str(identifier.split("#")[-1])
        if "/" in identifier:
            identifier = str(identifier.split("/")[-1])

        return identifier

    def process_secondary_files(self, secondary_files: List):
        if not hasattr(self.cwlgen, "SecondaryFileSchema"):
            return secondary_files

        return [
            s.pattern if isinstance(s, self.cwlgen.SecondaryFileSchema) else s
            for s in secondary_files
        ]

    single_token_matcher = re.compile("^\$\((.+)\)$")
    inline_expression_matcher = re.compile("\$\((.+)\)")
    input_selector_matcher = re.compile("^inputs\.(A-z0-9)$")

    def parse_basic_expression(self, expr):

        match = self.single_token_matcher.match(expr)
        if match:
            return self.convert_javascript_token(match.groups()[0])

        tokens = set(self.inline_expression_matcher.findall(expr))

        string_format = f"{expr}"
        token_replacers = {}

        for token, idx in zip(tokens, range(len(tokens))):
            key = f"JANIS_CWL_TOKEN_{idx+1}"
            string_format = string_format.replace(token, key)
            token_replacers[key] = self.convert_javascript_token(token)

        return j.StringFormatter(string_format, **token_replacers)

    def convert_javascript_token(self, token: str):
        input_selector_match = self.input_selector_matcher.match(token)
        if input_selector_match:
            return j.InputSelector(input_selector_match.group()[0])

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

        inp_type = self.from_cwl_type(inp.type)
        if inp.secondaryFiles:
            array_optional_layers = []
            while isinstance(inp_type, j.Array):
                array_optional_layers.append(inp_type.optional)
                inp_type = inp_type.subtype()

            inp_type = GenericFileWithSecondaries(
                secondaries=self.process_secondary_files(inp.secondaryFiles)
            )
            for is_optional in array_optional_layers[::-1]:
                inp_type = j.Array(inp_type, optional=is_optional)

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

        return j.ToolOutput(
            tag=self.get_tag_from_identifier(out.id),
            output_type=self.from_cwl_type(out.type),
            selector=selector,
        )

    def ingest_command_line_tool(self, clt):

        docker_requirement = None  # : Optional[self.cwlgen.DockerRequirement]
        for req in clt.requirements:
            if isinstance(req, self.cwlgen.DockerRequirement):
                docker_requirement = req

        container = None
        if docker_requirement:
            container = docker_requirement.dockerPull

        jclt = j.CommandToolBuilder(
            tool=self.get_tag_from_identifier(clt.id),
            base_command=clt.baseCommand,
            inputs=[self.ingest_command_tool_input(inp) for inp in clt.inputs],
            outputs=[self.ingest_command_tool_output(out) for out in clt.outputs],
            arguments=[self.ingest_command_tool_argument(arg) for arg in clt.arguments],
            version="v0.1.0",
            container=container or "ubuntu:latest",
        )
        return jclt

    @classmethod
    def load_cwlgen_from_version(cls, cwl_version: str):
        global cwlgen

        if cwl_version == "v1.0":
            import cwl_utils.parser_v1_0 as cwlutils
        elif cwl_version == "v1.1":
            import cwl_utils.parser_v1_1 as cwlutils
        elif cwl_version == "v1.2":
            import cwl_utils.parser_v1_2 as cwlutils
        else:
            print(
                f"Didn't recognise CWL version {cwl_version}, loading default: {DEFAULT_PARSER_VERSION}"
            )
            cwlutils = cls.load_cwlgen_from_version(DEFAULT_PARSER_VERSION)

        return cwlutils

    @classmethod
    def load_cwl_version_from_doc(cls, doc: str) -> str:
        import ruamel.yaml

        # load tool into memory
        with open(doc) as fp:
            tool_dict = ruamel.yaml.load(fp, Loader=ruamel.yaml.Loader)

        if "cwlVersion" not in tool_dict:
            raise Exception(f"Couldn't find cwlVersion in tool {doc}")

        return tool_dict["cwlVersion"]


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        raise Exception("Expected 1 argument, the name of a CWL tool.")
    toolname = sys.argv[1]

    tool = CWlParser.from_doc(toolname)

    tool.translate("wdl")
