import os
import copy
from abc import abstractmethod
from inspect import isclass
from typing import List, Type, Union, Optional, Dict, Tuple, Any

from janis_core import translations
from janis_core.utils import first_value
from janis_core.utils.metadata import WorkflowMetadata
from janis_core.utils.logger import Logger
from janis_core.graph.node import Node, NodeTypes
from janis_core.graph.steptaginput import StepTagInput
from janis_core.translations import ExportPathKeywords
from janis_core.types import (
    DataType,
    ParseableType,
    get_instantiated_type,
    Array,
    InputSelector,
)
from janis_core.tool.tool import Tool, ToolType, ToolTypes, ToolInput, ToolOutput
from janis_core.tool.commandtool import CommandTool
from janis_core.types.data_types import is_python_primitive
from janis_core.utils.scatter import ScatterDescription
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
        self,
        wf,
        identifier: str,
        datatype: DataType,
        default: any,
        value: any,
        doc: str = None,
    ):
        super().__init__(wf, NodeTypes.INPUT, identifier)
        self.datatype = datatype
        self.default = default
        self.doc = doc
        self.value = value

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {None: ToolOutput(self.identifier, self.datatype)}

    def inputs(self):
        return None


class StepNode(Node):
    def __init__(
        self,
        wf,
        identifier,
        tool: Tool,
        doc: str = None,
        scatter: ScatterDescription = None,
    ):
        super().__init__(wf, NodeTypes.STEP, identifier)
        self.tool = tool
        self.doc = doc
        self.scatter = scatter

    def inputs(self):
        return self.tool.inputs_map()

    def outputs(self):
        return self.tool.outputs_map()

    def _add_edge(self, tag: str, source: ConnectionSource):
        node, outtag = verify_or_try_get_source(source)

        if tag not in self.sources:
            self.sources[tag] = StepTagInput(self, tag)

        # If tag is in scatter.fields, then we can
        scatter = self.scatter and tag in self.scatter.fields

        return self.sources[tag].add_source(node, outtag, should_scatter=scatter)

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
        output_tag: Union[str, InputSelector, List[Union[str, InputSelector]]] = None,
        output_prefix: Union[str, InputSelector] = None,
    ):
        super().__init__(wf, NodeTypes.OUTPUT, identifier)
        self.datatype = datatype

        if source[0].node_type != NodeTypes.STEP:
            raise Exception(
                f"Unsupported connection type: {NodeTypes.OUTPUT} → {source[0].node_type}"
            )

        stype = source[0].outputs()[source[1]].output_type
        snode = source[0]
        if isinstance(snode, StepNode) and snode.scatter:
            stype = Array(stype)

        if not datatype.can_receive_from(stype):
            Logger.critical(
                f"Mismatch of types when joining to output node '{source[0].id()}.{source[1]}' to '{identifier}' "
                f"({stype.id()} -/→ {datatype.id()})"
            )

        self.source = verify_or_try_get_source(source)
        self.doc = doc
        self.output_tag = output_tag
        self.output_prefix = output_prefix

    def inputs(self) -> Dict[str, ToolInput]:
        # Program will just grab first value anyway
        return {None: ToolInput(self.identifier, self.datatype)}

    def outputs(self):
        return None


class Workflow(Tool):
    def __init__(self, **connections):
        super().__init__(metadata_class=WorkflowMetadata)

        self.connections = connections

        Logger.log(f"Creating workflow with identifier: '{self.id()}'")

        if not Validators.validate_identifier(self.id()):
            raise Exception(
                f"The identifier '{self.id()}' was invalid because {Validators.reason_for_failure(self.id())}"
            )

        # The following variables allow us to quickly check data about the graph
        self.nodes: Dict[str, Node] = {}

        self.input_nodes: Dict[str, InputNode] = {}
        self.step_nodes: Dict[str, StepNode] = {}
        self.output_nodes: Dict[str, OutputNode] = {}

        # Flags for different requirements that a workflow might need
        self.has_scatter = False
        self.has_subworkflow = False
        self.has_multiple_inputs = False

        # Now that we've initialised everything, we can "construct" the workflows for that subclass this class
        # else, for the WorkflowBuilder it will do nothing and they'll add workflows later
        self.constructor()

    @abstractmethod
    def friendly_name(self):
        pass

    @abstractmethod
    def constructor(self):
        """
        A place to construct your workflows. This is called directly after initialisation.
        :return:
        """
        pass

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
                f"The identifier '{identifier}' was invalid because {Validators.reason_for_failure(identifier)}"
            )

    def input(
        self,
        identifier: str,
        datatype: ParseableType,
        default: any = None,
        value: any = None,
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
            self,
            identifier=identifier,
            datatype=datatype,
            default=default,
            doc=doc,
            value=value,
        )
        self.nodes[identifier] = inp
        self.input_nodes[identifier] = inp
        return inp

    def output(
        self,
        identifier: str,
        datatype: Optional[ParseableType] = None,
        source: Union[StepNode, ConnectionSource] = None,
        output_tag: Union[
            str,
            InputSelector,
            ConnectionSource,
            List[Union[str, InputSelector, ConnectionSource]],
        ] = None,
        output_prefix: Union[str, InputSelector, ConnectionSource] = None,
    ):
        """
        Create an output on a workflow

        :param identifier: The identifier for the output
        :param datatype: Optional data type of the output to check. This will be automatically inferred if not provided.
        :param source: The source of the output, must be an output to a step node
        :param output_tag: A janis annotation for grouping outputs by this value.  If a list is passed, it represents
        a structure of nested directories, the first element being the root directory.
        At most, one InputSelector can resolve to an array, and this behaviour is only defined if the output
        scattered source, and the number of elements is equal.
        :param output_prefix: Decides the prefix that an output will have, or acts as a map if the InputSelector
        resolves to an array with equal length to the number of shards (scatters). Any other behaviour is defined and
        may result in an unexpected termination.
        :return:
        """
        self.verify_identifier(identifier, repr(datatype))

        if source is None:
            raise Exception("Output source must not be 'None'")

        node, tag = verify_or_try_get_source(source)
        if not datatype:
            datatype = node.outputs()[tag].output_type

            if isinstance(node, StepNode) and node.scatter:
                datatype = Array(datatype)

        if output_prefix:
            output_prefix = self.verify_output_source_type(
                identifier, output_prefix, "output_prefix"
            )
        if output_tag:
            output_tag = self.verify_output_source_type(
                identifier, output_tag, "output_tag"
            )

        otp = OutputNode(
            self,
            identifier=identifier,
            datatype=get_instantiated_type(datatype),
            source=(node, tag),
            output_tag=output_tag,
            output_prefix=output_prefix,
        )
        self.nodes[identifier] = otp
        self.output_nodes[identifier] = otp
        return otp

    def verify_output_source_type(
        self,
        identifier,
        out: Union[
            str,
            InputSelector,
            ConnectionSource,
            List[Union[str, InputSelector, ConnectionSource]],
        ],
        outtype: str,
    ):
        if isinstance(out, list):
            return [self.verify_output_source_type(o) for o in out]

        if isinstance(out, str):
            return out

        if isinstance(out, tuple):
            # ConnectionSource tuple
            out = out[0]

        if isinstance(out, Node):
            if not isinstance(out, InputNode):
                raise Exception(
                    f"The source for the {outtype} '{identifier}' was a {out.__class__.__name__} and must be an Input"
                )

            return InputSelector(out.identifier)

        if isinstance(out, InputSelector):
            keys = set(self.input_nodes.keys())
            if out.input_to_select not in keys:
                raise Exception(
                    f"Couldn't find the input {out.input_to_select} in the workflow, expected one of: "
                    + ", ".join(keys)
                )
            return out

        raise Exception(f"Invalidate type for {outtype}: {out.__class__.__name__}")

    def step(
        self,
        identifier: str,
        tool: Tool,
        scatter: Union[str, ScatterDescription] = None,
        ignore_missing=False,
    ):
        """
        Method to do
        :param identifier:
        :param tool:
        :param scatter: Indicate whether a scatter should occur, on what, and how
        :type scatter: Union[str, ScatterDescription]
        :param ignore_missing: Don't throw an error if required params are missing from this function
        :return:
        """

        self.verify_identifier(identifier, tool.id())

        if scatter is not None and isinstance(scatter, str):

            scatter = ScatterDescription(
                [scatter] if isinstance(scatter, str) else scatter
            )

        # verify scatter
        if scatter:
            ins = set(tool.inputs_map().keys())
            fields = set(scatter.fields)
            if any(f not in ins for f in fields):
                # if there is a field not in the input map, we have a problem
                extra_keys = ", ".join(f"'{f}'" for f in (fields - ins))
                raise Exception(
                    f"Couldn't scatter the field(s) for step '{identifier}' "
                    f"(tool: '{tool}') {extra_keys} as they were not found in the input map"
                )

        tool.workflow = self
        inputs = tool.inputs_map()

        connections = tool.connections

        provided_keys = set(connections.keys())
        all_keys = set(inputs.keys())
        required_keys = set(
            i for i, v in inputs.items() if not v.input_type.optional and not v.default
        )

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

        stp = StepNode(self, identifier=identifier, tool=tool, scatter=scatter)

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

            si = e.finish.sources[e.ftag] if e.ftag else first_value(e.finish.sources)
            self.has_multiple_inputs = self.has_multiple_inputs or si.multiple_inputs

        self.has_scatter = self.has_scatter or scatter is not None
        self.has_subworkflow = self.has_subworkflow or isinstance(tool, Workflow)
        self.nodes[identifier] = stp
        self.step_nodes[identifier] = stp

        return stp

    def __getattr__(self, item):
        if item in self.__dict__ or item == "nodes":
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
        tool_to_console=False,
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
            tool_to_console=tool_to_console,
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

    def generate_inputs_override(
        self, additional_inputs=None, with_resource_overrides=False, hints=None
    ):
        """
        Generate the overrides to be used with Janis. Although it may work with
        other
        :return:
        """
        ad = additional_inputs or {}

        d = {
            i.id(): ad.get(i.id(), i.value or i.default)
            for i in self.input_nodes.values()
            if i.id() in ad or i.default or i.value or not i.datatype.optional
        }

        if with_resource_overrides:
            d.update(translations.CwlTranslator().build_resources_input(self, hints))

        return d

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


class WorkflowBuilder(Workflow):
    def __init__(self, identifier: str = None, friendly_name: str = None):
        self._identifier = identifier
        self._name = friendly_name

        super().__init__()

    def id(self):
        return self._identifier

    def friendly_name(self):
        return self._name

    def __call__(self, **connections):
        self.connections = connections
        return self

    def constructor(self):
        """
        Empty placeholder as users will construct their workflow manually, and not as part of this class
        :return:
        """
        return self


if __name__ == "__main__":
    # from janis_unix import Echo

    wf = WorkflowBuilder("simple")

    from janis_bioinformatics.tools.cutadapt import CutAdapt_1_18 as Cutadapt
    from janis_bioinformatics.tools.common import BwaMem_SamToolsView
    from janis_bioinformatics.tools.gatk4 import Gatk4SortSam_4_0
    from janis_bioinformatics.data_types import FastaWithDict, Fastq

    # Inputs
    wf.input("sampleName", str)
    wf.input("reference", FastaWithDict)
    wf.input("fastq", Fastq)

    # Steps
    wf.step("cutadapt", Cutadapt(fastq=wf.fastq))
    wf.step(
        "bwamem",
        BwaMem_SamToolsView(
            reads=wf.cutadapt.out, sampleName=wf.sampleName, reference=wf.reference
        ),
    )
    wf.step(
        "sortsam",
        Gatk4SortSam_4_0(
            bam=wf.bwamem.out,
            sortOrder="coordinate",
            createIndex=True,
            validationStringency="SILENT",
            maxRecordsInRam=5000000,
            tmpDir=".",
        ),
    )

    # outputs
    wf.output("out", source=wf.sortsam)

    wf.translate("wdl")
