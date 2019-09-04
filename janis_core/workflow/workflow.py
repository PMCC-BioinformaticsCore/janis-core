import os
import copy
from inspect import isclass
from typing import List, Type, Union, Optional, Dict, Tuple, Any

from janis_core import translations
from janis_core.utils import first_value
from janis_core.utils.metadata import WorkflowMetadata
from janis_core.utils.logger import Logger
from janis_core.graph.node import Node, NodeTypes
from janis_core.graph.stepinput import StepInput
from janis_core.translations import ExportPathKeywords
from janis_core.types import DataType, ParseableType, get_instantiated_type
from janis_core.tool.tool import Tool, ToolType, ToolTypes, ToolInput, ToolOutput
from janis_core.tool.commandtool import CommandTool
from janis_core.types.data_types import is_python_primitive
from janis_core.utils.validators import Validators

ConnectionSource = Union[Node, Tuple[Node, str]]


def verify_or_try_get_source(source: Union[ConnectionSource, List[ConnectionSource]]):
    if isinstance(source, list):
        return [verify_or_try_get_source(s) for s in source]
    node, tag = None, None
    if isinstance(source, tuple):
        node, tag = source
    else:
        node = source

    outs = node.outputs()
    if tag is None:
        if len(outs) > 1:
            raise Exception(
                f"Too many outputs of {node.id()} to guess the correct output"
            )
        tag = list(outs.keys())[0]

    if tag not in outs:
        tags = ", ".join([f"out.{o}" for o in outs.keys()])
        raise Exception(
            f"Couldn't find tag '{tag}' in outputs of '{node.id()}', "
            f"expected one of {tags}"
        )

    return node, tag


class InputNode(Node):
    def __init__(
        self, wf, identifier: str, datatype: DataType, default: any, doc: str = None
    ):
        super().__init__(wf, NodeTypes.INPUT, identifier)
        self.datatype = datatype
        self.default = default
        self.doc = doc

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {None: ToolOutput(self.identifier, self.datatype)}

    def inputs(self):
        return None


class StepNode(Node):
    def __init__(self, wf, identifier, tool: Tool, doc: str = None):
        super().__init__(wf, NodeTypes.STEP, identifier)
        self.tool = tool
        self.doc = doc

    def inputs(self):
        return self.tool.inputs_map()

    def outputs(self):
        return self.tool.outputs_map()

    def _add_edge(self, tag: str, source: ConnectionSource):
        node, outtag = verify_or_try_get_source(source)

        if tag not in self.sources:
            self.sources[tag] = StepInput(self, tag)

        return self.sources[tag].add_source(node, outtag)

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]

        return self.get_item(item)

    def __getitem__(self, item):
        return self.get_item(item)

    def get_item(self, item):
        ins = self.inputs()
        if item in ins:
            return self, item
        outs = self.outputs()
        if item in outs:
            return self, item

        tags = ", ".join(
            [f"in.{i}" for i in ins.keys()] + [f"out.{o}" for o in outs.keys()]
        )

        raise KeyError(
            f"Step '{self.id()}' with tool '{self.tool.id()}' has no identifier '{item}' ({tags})"
        )

    # "always set attributes keys"
    always_set = {"tool", "sources", "node_type", "identifier"}

    def __setitem__(self, key, value):
        print("Attempting to set: " + key)
        if key not in StepNode.always_set and key in self.inputs():
            if isinstance(value, list):
                return [self._add_edge(key, v) for v in value]
            return self._add_edge(key, value)

        self.__dict__[key] = value


class OutputNode(Node):
    def __init__(
        self,
        wf,
        identifier: str,
        datatype: DataType,
        source: ConnectionSource,
        doc: str = None,
    ):
        super().__init__(wf, NodeTypes.OUTPUT, identifier)
        self.datatype = datatype

        if source[0].node_type != NodeTypes.STEP:
            raise Exception(
                f"Unsupported connection type: {NodeTypes.OUTPUT} → {source[0].node_type}"
            )

        stype = source[0].outputs()[source[1]].output_type
        if not datatype.can_receive_from(stype):
            Logger.critical(
                f"Mismatch of types when joining '{source[0].id()}.{source[1]}' to '{identifier}' "
                f"({stype.id()} -/→ {datatype.id()})"
            )

        self.source = verify_or_try_get_source(source)
        self.doc = doc

    def inputs(self) -> Dict[str, ToolInput]:
        # Program will just grab first value anyway
        return {None: ToolInput(self.identifier, self.datatype)}

    def outputs(self):
        return None


class Workflow(Tool):
    def __init__(self, identifier: str, name: str = None):
        super().__init__()

        self.metadata = WorkflowMetadata()

        Logger.log(f"Creating workflow with identifier: '{identifier}'")

        if not Validators.validate_identifier(identifier):
            raise Exception(
                f"The identifier '{identifier}' was not validated by '{Validators.identifier_regex}' "
                f"(must start with letters, and then only contain letters, numbers and an underscore)"
            )

        self._identifier = identifier
        self._name = name

        # The following variables allow us to quickly check data about the graph
        self.nodes: Dict[str, Node] = {}

        self.input_nodes: Dict[str, InputNode] = {}
        self.step_nodes: Dict[str, StepNode] = {}
        self.output_nodes: Dict[str, OutputNode] = {}

        # Flags for different requirements that a workflow might need
        self.has_scatter = False
        self.has_subworkflow = False
        self.has_multiple_inputs = False

    def verify_identifier(self, identifier: str, component: str):

        if identifier in self.__dict__:
            raise Exception(
                f"'{identifier}' is a protected keyword for a janis workflow"
            )

        if identifier in self.nodes:
            existing = self.nodes[identifier]
            raise Exception(
                f"There already exists a node (and component) with id '{identifier}'. The added "
                f"component ('{component}') clashes with '{repr(existing)}')."
            )

        if not Validators.validate_identifier(identifier):
            raise Exception(
                f"The identifier '{identifier}' was not validated by '{Validators.identifier_regex}' "
                f"(must start with letters, and then only contain letters, numbers and an underscore)"
            )

    def input(
        self,
        identifier: str,
        datatype: ParseableType,
        default: any = None,
        doc: str = None,
    ):
        """
        Create an input node on a workflow
        :return:
        """

        self.verify_identifier(identifier, repr(datatype))

        datatype = get_instantiated_type(datatype)
        if default:
            datatype.optional = True

        inp = InputNode(
            self, identifier=identifier, datatype=datatype, default=default, doc=doc
        )
        self.nodes[identifier] = inp
        self.input_nodes[identifier] = inp
        return inp

    def output(
        self,
        identifier: str,
        datatype: Optional[ParseableType] = None,
        source: Union[StepNode, ConnectionSource] = None,
    ):
        self.verify_identifier(identifier, repr(datatype))

        if source is None:
            raise Exception("Output source must not be 'None'")

        node, tag = verify_or_try_get_source(source)
        if not datatype:
            datatype = node.outputs()[tag].output_type

        otp = OutputNode(
            self,
            identifier=identifier,
            datatype=get_instantiated_type(datatype),
            source=(node, tag),
        )
        self.nodes[identifier] = otp
        self.output_nodes[identifier] = otp
        return otp

    def step(
        self,
        identifier: str,
        tool: Union[Tool, Type[Tool]],
        ignore_missing=False,
        **connections,
    ):
        if isclass(tool) and issubclass(tool, Tool):
            tool = tool()

        self.verify_identifier(identifier, tool.id())

        tool.workflow = self
        inputs = tool.inputs_map()

        provided_keys = set(connections.keys())
        all_keys = set(inputs.keys())
        required_keys = set(i for i, v in inputs.items() if not v.input_type.optional)

        if not provided_keys.issubset(all_keys):
            unrecparams = ", ".join(provided_keys - all_keys)

            tags = ", ".join([f"in.{i}" for i in all_keys])

            raise Exception(
                f"Unrecognised parameters {unrecparams} when creating '{identifier}' ({tool.id()}). "
                f"Expected types: {tags}"
            )

        if not ignore_missing and not required_keys.issubset(provided_keys):
            missing = ", ".join(f"'{i}'" for i in (required_keys - provided_keys))
            raise Exception(
                f"Missing the parameters {missing} when creating '{identifier}' ({tool.id()})"
            )

        stp = StepNode(self, identifier=identifier, tool=tool)

        added_edges = []
        for (k, v) in connections.items():

            if is_python_primitive(v):
                inp_identifier = f"{identifier}_{k}"

                referencedtype = copy.copy(inputs[k].input_type)
                parsed_type = get_instantiated_type(v)

                if parsed_type and not referencedtype.can_receive_from(parsed_type):
                    raise TypeError(
                        f"The type {parsed_type.id()} inferred from the value '{v}' is not "
                        f"compatible with the '{identifier}.{k}' type: {referencedtype.id()}"
                    )

                referencedtype.optional = True

                v = self.input(inp_identifier, referencedtype, default=v)
            if v is None:
                inp_identifier = f"{identifier}_{k}"
                v = self.input(inp_identifier, inputs[k].input_type, default=v)

            verifiedsource = verify_or_try_get_source(v)
            if isinstance(verifiedsource, list):
                for vv in verifiedsource:
                    added_edges.append(stp._add_edge(k, vv))
            else:
                added_edges.append(stp._add_edge(k, verifiedsource))

        for e in added_edges:
            self.has_scatter = self.has_scatter or e.scatter
            si = e.finish.sources[e.ftag] if e.ftag else first_value(e.finish.sources)
            self.has_multiple_inputs = self.has_multiple_inputs or si.multiple_inputs

        self.nodes[identifier] = stp
        self.step_nodes[identifier] = stp
        return stp

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__.get(item)

        if item in self.nodes:
            return self.nodes[item]

        raise AttributeError(
            f"AttributeError: '{type(self).__name__}' object has no attribute '{item}'"
        )

    def __getitem__(self, item):

        if item in self.nodes:
            return self.nodes[item]

        raise KeyError(f"KeyError: '{type(self).__name__}' object has no node '{item}'")

    @classmethod
    def type(cls) -> ToolType:
        return ToolTypes.Workflow

    def id(self) -> str:
        return self._identifier

    def friendly_name(self):
        return self._name

    def inputs(self) -> List[ToolInput]:
        """
        List of ToolInputs of the workflow, we can toss out most of the metadata
        about positioning, prefixes, etc that the ToolInput class uses
        """
        return [ToolInput(i.id(), i.datatype) for i in self.input_nodes.values()]

    def outputs(self) -> List[ToolOutput]:
        """
        Similar to inputs, return a list of ToolOutputs of the workflow
        """
        return [ToolOutput(i.id(), i.datatype) for i in self.output_nodes.values()]

    @staticmethod
    def version():
        return "v0.1.0"

    def translate(
        self,
        translation: translations.SupportedTranslation,
        to_console=True,
        to_disk=False,
        write_inputs_file=True,
        with_docker=True,
        with_hints=False,
        with_resource_overrides=False,
        validate=False,
        should_zip=True,
        export_path=ExportPathKeywords.default,
        merge_resources=False,
        hints=None,
        allow_null_if_not_optional=True,
        additional_inputs: Dict = None,
        max_cores=None,
        max_mem=None,
    ):
        return translations.translate_workflow(
            self,
            translation=translation,
            to_console=to_console,
            to_disk=to_disk,
            with_docker=with_docker,
            with_resource_overrides=with_resource_overrides,
            should_zip=should_zip,
            export_path=export_path,
            write_inputs_file=write_inputs_file,
            should_validate=validate,
            merge_resources=merge_resources,
            hints=hints,
            allow_null_if_not_optional=allow_null_if_not_optional,
            additional_inputs=additional_inputs,
            max_cores=max_cores,
            max_mem=max_mem,
        )

    def generate_inputs_override(self):
        """
        Generate the overrides to be used with Janis. Although it may work with
        other
        :return:
        """
        return {
            i.id(): i.default
            for i in self.input_nodes.values()
            if i.default or not i.datatype.optional
        }

    def generate_resources_file(
        self,
        translation: translations.SupportedTranslation,
        hints: Dict[str, Any] = None,
        to_console=True,
    ):
        tr = translations.build_resources_input(self, translation, hints)
        if to_console:
            print(tr)
        return tr

    def get_tools(self) -> Dict[str, CommandTool]:
        tools: Dict[str, CommandTool] = {}
        for t in self._steps:
            tl = t.step.tool()
            if isinstance(tl, Workflow):
                tools.update(tl.get_tools())
            elif t.id() not in tools:
                tools[tl.id()] = tl
        return tools

    def report(self, to_console=True, tabulate_tablefmt=None):
        import tabulate

        tools = self.get_tools()
        keys = sorted(tools.keys(), key=lambda a: a[0].lower())

        header = ["tool", "version", "container"]
        data = []
        for t in keys:
            tool = tools[t]
            data.append(
                [
                    f"{tool.friendly_name()} ({tool.id()})",
                    tool.version(),
                    tool.container(),
                ]
            )

        retval = tabulate.tabulate(data, headers=header, tablefmt=tabulate_tablefmt)
        if to_console:
            print(retval)

        return retval

    def generate_resources_table(
        self,
        hints: Dict[str, Any],
        to_console=True,
        to_disk=False,
        output_type: str = "tsv",
    ):
        delim = "\t" if output_type == "tsv" else ","

        tools = self.get_tools()
        header = ["name", "cpu", "memory (GB)"]
        data = []

        for t in sorted(tools.keys()):
            tool = tools[t]
            data.append([tool.id(), tool.cpus(hints), tool.memory(hints)])

        data.sort(key=lambda a: a[0].lower())

        data.insert(0, header)

        if to_console:
            import tabulate

            print(tabulate.tabulate(data, headers="firstrow"))
        if to_disk:
            import csv

            d = ExportPathKeywords.resolve(
                ExportPathKeywords.default_no_spec, None, self.id()
            )
            path = d + f"resources.{output_type}"

            if not os.path.isdir(d):
                os.makedirs(d)

            Logger.info(f"Writing resources {output_type} to '{path}'")
            with open(path, "w+") as mf:
                writer = csv.writer(mf, delimiter=delim)
                for row in data:
                    writer.writerow(row)

        return data


if __name__ == "__main__":
    # from janis_unix import Echo

    wf = Workflow("simple")

    from janis_bioinformatics.tools.cutadapt import CutAdapt_1_18 as Cutadapt
    from janis_bioinformatics.tools.common import BwaMem_SamToolsView
    from janis_bioinformatics.tools.gatk4 import Gatk4SortSam_4_0
    from janis_bioinformatics.data_types import FastaWithDict, Fastq

    # Inputs
    wf.input("sampleName", str)
    wf.input("reference", FastaWithDict)
    wf.input("fastq", Fastq)

    # Steps
    wf.step("cutadapt", Cutadapt, fastq=wf.fastq)
    wf.step(
        "bwamem",
        BwaMem_SamToolsView,
        reads=wf.cutadapt.out,
        sampleName=wf.sampleName,
        reference=wf.reference,
    )
    wf.step(
        "sortsam",
        Gatk4SortSam_4_0,
        bam=wf.bwamem.out,
        sortOrder="coordinate",
        createIndex=True,
        validationStringency="SILENT",
        maxRecordsInRam=5000000,
        tmpDir=".",
    )

    # outputs
    wf.output("out", source=wf.sortsam)

    wf.translate("wdl")
