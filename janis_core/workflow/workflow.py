import copy
import os
from abc import abstractmethod, ABC
from inspect import isclass
from typing import List, Union, Optional, Dict, Tuple, Any, Set, Iterable, Type

from janis_core.graph.node import Node, NodeType
from janis_core.graph.steptaginput import StepTagInput
from janis_core.operators import (
    InputSelector,
    Operator,
    StringFormatter,
    StepOutputSelector,
    InputNodeSelector,
    Selector,
    AliasSelector,
)
from janis_core.operators.logical import AndOperator, NotOperator, or_prev_conds
from janis_core.operators.standard import FirstOperator
from janis_core.tool.commandtool import CommandTool
from janis_core.tool.documentation import (
    InputDocumentation,
    OutputDocumentation,
    InputQualityType,
    DocumentationMeta,
)
from janis_core.tool.tool import Tool, ToolType, TInput, TOutput
from janis_core.translationdeps.exportpath import ExportPathKeywords
from janis_core.translationdeps.supportedtranslations import SupportedTranslation
from janis_core.types import (
    DataType,
    ParseableType,
    get_instantiated_type,
    Array,
    Filename,
)
from janis_core.types.data_types import is_python_primitive
from janis_core.utils import first_value, fully_qualify_filename
from janis_core.utils.logger import Logger
from janis_core.utils.metadata import WorkflowMetadata
from janis_core.utils.scatter import ScatterDescription, ScatterMethod
from janis_core.utils.validators import Validators

ConnectionSource = Union[Node, StepOutputSelector, Tuple[Node, str]]


def verify_or_try_get_source(
    source: Union[ConnectionSource, List[ConnectionSource]]
) -> Union[StepOutputSelector, InputNodeSelector, List[StepOutputSelector], Operator]:

    if isinstance(source, StepOutputSelector):
        return source
    elif isinstance(source, InputNodeSelector):
        return source
    elif isinstance(source, AliasSelector):
        return source
    elif isinstance(source, list):
        return [verify_or_try_get_source(s) for s in source]

    if isinstance(source, Operator):
        # return [verify_or_try_get_source(r) for r in rt] if len(rt) > 1 else rt
        return source
    node, tag = None, None
    if isinstance(source, tuple):
        node, tag = source
    elif isinstance(source, Node):
        node = source
    else:
        raise Exception(f"Unrecognised source type: {source} ({type(source).__name__})")

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

    return StepOutputSelector(node, tag)


class InputNode(Node):
    def __init__(
        self,
        wf,
        identifier: str,
        datatype: DataType,
        default: any,
        value: any,
        doc: InputDocumentation = None,
    ):
        super().__init__(wf, NodeType.INPUT, identifier)
        self.datatype = datatype
        self.default = default
        self.doc: Optional[InputDocumentation] = doc
        self.value = value

    def as_operator(self):
        return InputNodeSelector(self)

    def outputs(self) -> Dict[Optional[str], TOutput]:
        # Program will just grab first value anyway
        return {None: TOutput(self.identifier, self.datatype)}

    def inputs(self):
        return None


class StepNode(Node):
    def __init__(
        self,
        wf,
        identifier,
        tool: Tool,
        doc: DocumentationMeta = None,
        scatter: ScatterDescription = None,
        when: Operator = None,
    ):
        super().__init__(wf, NodeType.STEP, identifier)
        self.tool = tool
        self.doc = doc
        self.scatter = scatter
        self.when = when

        self.parent_has_conditionals = False
        self.has_conditionals = when is not None

    def inputs(self) -> Dict[str, TInput]:
        ins = self.tool.inputs_map()

        # if self.parent_has_conditionals:
        #     q = {}
        #     for iv in ins.values():
        #         intype = copy.copy(iv.intype)
        #         intype.optional = True
        #         q[iv.id()] = TInput(
        #             iv.id(), intype=intype, doc=iv.doc, default=iv.default
        #         )
        #
        #     ins = q

        return ins

    def outputs(self) -> Dict[str, TOutput]:
        outs = self.tool.outputs_map()

        if self.has_conditionals:
            q = {}
            for ov in outs.values():
                outtype = copy.copy(ov.outtype)
                outtype.optional = True
                q[ov.id()] = TOutput(ov.id(), outtype=outtype, doc=ov.doc)

            outs = q

        return outs

    def _add_edge(self, tag: str, source: ConnectionSource):
        stepoperator = verify_or_try_get_source(source)
        # node, outtag = stepoperator.node, stepoperator.tag

        # Todo: readd this when dealing with conditionals, essentially,
        #  need to determine whether an upstream is conditional (and hence this will have optional out)
        # if isinstance(node, StepNode):
        #     if node.has_conditionals or node.parent_has_conditionals:
        #         self.parent_has_conditionals = True

        if tag not in self.sources:
            self.sources[tag] = StepTagInput(self, tag)

        # If tag is in scatter.fields, then we can
        scatter = self.scatter and tag in self.scatter.fields

        return self.sources[tag].add_source(source, should_scatter=scatter)

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]

        return self.get_item(item)

    def __getitem__(self, item) -> StepOutputSelector:
        return self.get_item(item)

    def get_item(self, item) -> StepOutputSelector:
        ins = self.inputs()
        if item in ins:
            return StepOutputSelector(self, item)
        outs = self.outputs()
        if item in outs:
            return StepOutputSelector(self, item)

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
        source: Union[List[ConnectionSource], ConnectionSource],
        doc: OutputDocumentation = None,
        output_folder: Union[
            str, InputSelector, List[Union[str, InputSelector]]
        ] = None,
        output_name: Union[str, InputSelector] = None,
        extension: Optional[str] = None,
        skip_typecheck=False,
    ):
        super().__init__(wf, NodeType.OUTPUT, identifier)
        self.datatype = datatype

        sources = source if isinstance(source, list) else [source]
        single_source = sources[0]

        if isinstance(single_source, StepOutputSelector):
            snode = single_source.node
            stype = snode.outputs()[single_source.tag].outtype
            if snode.scatter:
                stype = Array(stype)
        elif isinstance(single_source, Selector):
            stype = single_source.returntype()
        else:
            raise Exception("Unsupported output source type " + str(single_source))

        if not skip_typecheck and not datatype.can_receive_from(stype):
            Logger.critical(
                f"Mismatch of types when joining to output node '{source.node.id()}.{source.tag}' to '{identifier}' "
                f"({stype.id()} -/→ {datatype.id()})"
            )

        self.source = verify_or_try_get_source(source)
        self.doc = (
            doc
            if isinstance(doc, OutputDocumentation)
            else OutputDocumentation(doc=doc)
        )
        self.output_folder = output_folder
        self.output_name = output_name
        self.extension = extension

    def inputs(self) -> Dict[str, TInput]:
        # Program will just grab first value anyway
        return {None: TInput(self.identifier, self.datatype)}

    def outputs(self):
        return None


class WorkflowBase(Tool, ABC):
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

    @abstractmethod
    def friendly_name(self):
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
        doc: Union[str, InputDocumentation, Dict[str, any]] = None,
    ):
        """
        Create an input node on a workflow
        :return:
        """

        self.verify_identifier(identifier, repr(datatype))

        datatype = get_instantiated_type(datatype)
        if default is not None:
            datatype.optional = True

        inp = InputNode(
            self,
            identifier=identifier,
            datatype=datatype,
            default=default,
            doc=InputDocumentation.try_parse_from(doc),
            value=value,
        )
        self.nodes[identifier] = inp
        self.input_nodes[identifier] = inp
        return InputNodeSelector(inp)

    def output(
        self,
        identifier: str,
        datatype: Optional[ParseableType] = None,
        source: Union[
            List[Union[Selector, ConnectionSource]], Union[Selector, ConnectionSource]
        ] = None,
        output_folder: Union[str, Selector, List[Union[str, Selector]]] = None,
        output_name: Union[bool, str, Selector, ConnectionSource] = True,
        extension: Optional[str] = None,
        doc: Union[str, OutputDocumentation] = None,
    ):
        """
        Create an output on a workflow

        :param identifier: The identifier for the output
        :param datatype: Optional data type of the output to check. This will be automatically inferred if not provided.
        :param source: The source of the output, must be an output to a step node
        :param output_folder: Decides the output folder(s) where the output will reside. If a list is passed, it
            represents a structure of nested directories, the first element being the root directory.
                - None (default): the assistant will copy to the root of the output directory
                - Type[Selector]: will be resolved before the workflow is run, this means it may only depend on the inputs
            NB: If the output_source is an array, a "shard_n" will be appended to the output_name UNLESS the output_source
            also resolves to an array, which the assistant can unwrap multiple dimensions of arrays ONLY if the number
            of elements in the output_scattered source and the number of resolved elements is equal.

        :param output_name: Decides the name of the output (without extension) that an output will have:
                - True (default): the assistant will choose an output name based on output identifier (tag),
                - None / False: the assistant will use the original filename (this might cause filename conflicts)
                - Type[Selector]: will be resolved before the workflow is run, this means it may only depend on the inputs
            NB: If the output_source is an array, a "shard_n" will be appended to the output_name UNLESS the output_source
                also resolves to an array, which the assistant can unwrap multiple dimensions of arrays.
        :param extension: The extension to use if janis renames the output. By default, it will pull the extension
            from the inherited data type (eg: CSV -> ".csv"), or it will attempt to pull the extension from the file.
        :return: janis.WorkflowOutputNode
        """
        self.verify_identifier(identifier, repr(datatype))

        if source is None:
            raise Exception("Output source must not be 'None'")

        sourceoperator = verify_or_try_get_source(source)
        skip_typecheck = False
        # if isinstance(sourceoperator, list):
        #     sourceoperator = sourceoperator[0]
        # else:
        #     sourceoperator = sourceoperator

        # if isinstance(sourceoperator, StepOutputSelector):
        #     node, tag = sourceoperator.node, sourceoperator.tag
        # elif isinstance(sourceoperator, InputNodeSelector):
        #     node = sourceoperator.input_node
        #     tag = None
        # else:
        #     raise Exception("Unsupported output source type: " + str(sourceoperator))

        if not datatype:
            while isinstance(sourceoperator, list):
                sourceoperator: Selector = sourceoperator[0]

            datatype: DataType = copy.copy(
                get_instantiated_type(sourceoperator.returntype()).received_type()
            )
            if (
                isinstance(sourceoperator, InputNodeSelector)
                and sourceoperator.input_node.default is not None
            ):
                datatype.optional = False

            elif isinstance(sourceoperator, StepNode) and sourceoperator.scatter:
                datatype = Array(datatype)

            skip_typecheck = True

        if output_name is not None:
            if isinstance(output_name, list):
                raise Exception("An output_name cannot be of type 'list'")
            output_name = self.verify_output_source_type(
                identifier, output_name, "output_name"
            )
        if output_folder is not None:
            ot = output_folder if isinstance(output_folder, list) else [output_folder]
            output_folder = self.verify_output_source_type(
                identifier, ot, "output_folder"
            )
        doc = (
            doc
            if isinstance(doc, OutputDocumentation)
            else OutputDocumentation(doc=doc)
        )

        otp = OutputNode(
            self,
            identifier=identifier,
            datatype=get_instantiated_type(datatype),
            source=sourceoperator,
            output_folder=output_folder,
            output_name=output_name,
            extension=extension,
            doc=doc,
            skip_typecheck=skip_typecheck,
        )
        self.nodes[identifier] = otp
        self.output_nodes[identifier] = otp
        return otp

    def forward_inputs_from_tool(
        self,
        tool: Union[Tool, Type[Tool]],
        inputs_to_forward: Optional[Iterable[str]] = None,
        inputs_to_ignore: Iterable[str] = None,
        input_prefix: str = "",
    ) -> Dict[str, InputNodeSelector]:
        """
        Usage:
            yourstp_inputs = self.forward_inputs_from_tool(YourTool, ["inp1", "inp2"])
                OR
            yourstp_inputs = self.forward_inputs_from_tool(YourTool, inputs_to_ignore=["inp2"])

            self.step("yourstp", YourTool(**yourstp_inputs))

        :param tool: The tool for which to forward the inputs for
        :param inputs_to_forward: List of inputs to forward. You MUST specify ALL the inputs you want to forward.
        :param inputs_to_ignore: You can choose _ALL_ inputs, except those specified here. This option is IGNORED if inputs_to_forward is defined
        :param input_prefix: Add a prefix when forwarding it to the parent workflow (self)
        :return:
        """

        itool = tool() if isclass(tool) else tool
        qualified_inputs_to_forward = inputs_to_forward
        tinps: Dict[str, TInput] = {t.id(): t for t in itool.tool_inputs()}

        if inputs_to_forward is None and inputs_to_ignore is None:
            raise Exception(
                f"You must specify ONE of inputs_to_forward OR inputs_to_ignore when "
                f"calling 'forward_inputs_from_tool' with tool {tool.id()})"
            )
        elif not inputs_to_forward:
            qualified_inputs_to_forward = [
                inpid for inpid in tinps.keys() if inpid not in inputs_to_ignore
            ]

        d = {}
        for inp in qualified_inputs_to_forward:
            if inp not in tinps:
                raise Exception(
                    f"Couldn't find the input {inp} in the tool {itool.id()}"
                )
            tinp = tinps[inp]
            d[inp] = self.input(input_prefix + inp, tinp.intype, doc=tinp.doc)

        return d

    def capture_outputs_from_step(
        self,
        step: StepNode,
        output_prefix=None,
        default_output_name: Union[bool, str, Selector, ConnectionSource] = True,
    ):
        op = output_prefix or ""

        tool = step.tool
        input_ids = set(self.input_nodes.keys())

        # Selector with 3 possibilities
        #   1. Valid - no requirements -> return true
        #   2. Almost valid - requires a transformation -> {}
        #   3. Invalid - uses an invalid type -> return false

        def check_selector_and_get_transformation(
            selector, output_id: str
        ) -> Union[Selector, bool]:
            """
            :return:
                - True if no transformation is required
                - None if the selector is invalid
                - Selector if the transformation is required
            """

            if is_python_primitive(selector):
                return True
            elif isinstance(selector, InputSelector):
                if selector.input_to_select not in input_ids:
                    Logger.warn(
                        f"Couldn't port the through selector for for {step.id()}.{output_id} as the input "
                        f"'{selector.input_to_select}' was not found in the current inputs"
                    )
                    return False
                return True

            elif isinstance(selector, InputNodeSelector):
                if selector.id() not in input_ids:
                    Logger.warn(
                        f"Couldn't port the through selector for for {step.id()}.{output_id} as an input node with ID"
                        f"'{selector.id()}' was not found in the current inputs"
                    )
                    return False
                return InputSelector(selector.id())

            Logger.warn(
                f"The selector {selector} ({type(selector)}) could not be transformed and cannot be used as an output/name folder for {self.id()}.{output_id}"
            )
            return False

        def get_transformed_selector(selector, output_id: str):

            if isinstance(selector, list):
                return [get_transformed_selector(sel, output_id) for sel in selector]

            if isinstance(selector, Operator):
                transformation = {}
                for sel in selector.get_leaves():
                    tr = check_selector_and_get_transformation(sel, output_id)
                    if tr is True:
                        continue
                    if tr is False:
                        # something was invalid and couldn't be transformed
                        # it's been logged though
                        return None
                    transformation[sel] = tr

                if len(transformation) > 0:
                    return selector.rewrite_operator(transformation)
                return selector

            tr = check_selector_and_get_transformation(selector, output_id)
            if tr is False:
                return None
            elif tr is True:
                return selector
            else:
                return tr

        for out in tool.tool_outputs():
            output_folders = None
            output_name = default_output_name
            ext = None
            if isinstance(tool, Workflow):
                outnode = tool.output_nodes[out.id()]
                ext = outnode.extension

                if outnode.output_folder is not None:
                    output_folders = get_transformed_selector(
                        outnode.output_folder, out.id()
                    )

                if outnode.output_name is not None:
                    output_name = get_transformed_selector(
                        outnode.output_name, out.id()
                    )

            self.output(
                op + out.id(),
                source=step[out.id()],
                output_name=output_name,
                output_folder=output_folders or [],
                extension=ext,
            )

    def all_input_keys(self):
        from janis_core.translations.translationbase import TranslatorBase

        return super().all_input_keys() + list(
            TranslatorBase.build_resources_input(tool=self, hints={}).keys()
        )

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
            return [self.verify_output_source_type(identifier, o, outtype) for o in out]
        if isinstance(out, bool):
            return out
        if isinstance(out, (str, float, int)):
            return str(out)

        if isinstance(out, tuple):
            # ConnectionSource tuple
            out = out[0]

        if isinstance(out, Node):
            if not isinstance(out, InputNode):
                raise Exception(
                    f"The source for the {outtype} '{identifier}' was a {out.__class__.__name__} and must be an Input"
                )

            return InputSelector(out.identifier)

        if isinstance(out, InputNodeSelector):
            return InputSelector(out.input_node.id())

        if isinstance(out, InputSelector):
            keys = set(self.input_nodes.keys())
            if out.input_to_select not in keys:
                raise Exception(
                    f"Couldn't find the input {out.input_to_select} in the workflow, expected one of: "
                    + ", ".join(keys)
                )
            return out

        if isinstance(out, Operator):
            fixed_types = [
                self.verify_output_source_type(identifier, o, outtype)
                for o in out.get_leaves()
            ]
            return out

        raise Exception(f"Invalid type for {outtype}: {out.__class__.__name__}")

    def step(
        self,
        identifier: str,
        tool: Tool,
        scatter: Union[str, List[str], ScatterDescription] = None,
        when: Optional[Operator] = None,
        ignore_missing=False,
        doc: str = None,
    ):
        """
        Construct a step on this workflow.

        :param identifier: The identifier of the step, unique within the workflow.
        :param tool: The tool that should run for this step.
        :param scatter: Indicate whether a scatter should occur, on what, and how.
        :type scatter: Union[str, ScatterDescription]
        :param when: An operator / condition that determines whether the step should run
        :type when: Optional[Operator]
        :param ignore_missing: Don't throw an error if required params are missing from this function
        :return:
        """

        self.verify_identifier(identifier, tool.id())

        if scatter is not None and not isinstance(scatter, ScatterDescription):

            fields = None
            if isinstance(scatter, str):
                fields = [scatter]
            elif isinstance(scatter, list):
                fields = scatter
            else:
                raise Exception(
                    f"Couldn't scatter with field '{scatter}' ({type(scatter)}"
                )

            scatter = ScatterDescription(fields, method=ScatterMethod.dot)

        # verify scatter
        if scatter:
            ins = set(tool.inputs_map().keys())
            fields = set(scatter.fields)
            if any(f not in ins for f in fields):
                # if there is a field not in the input map, we have a problem
                extra_keys = ", ".join(f"'{f}'" for f in (fields - ins))
                raise Exception(
                    f"Couldn't scatter the field(s) {extra_keys} for step '{identifier}' "
                    f"as they are not inputs to the tool '{tool.id()}'"
                )

        tool.workflow = self
        inputs = tool.inputs_map()

        connections = tool.connections

        provided_keys = set(connections.keys())
        all_keys = set(inputs.keys())
        required_keys = set(
            # The input is optional if it's optional or has default)
            i
            for i, v in inputs.items()
            if not (v.intype.optional or v.default is not None)
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

        d = doc if isinstance(doc, DocumentationMeta) else DocumentationMeta(doc=doc)
        stp = StepNode(
            self, identifier=identifier, tool=tool, scatter=scatter, when=when, doc=d
        )

        added_edges = []
        for (k, v) in connections.items():

            isfilename = isinstance(v, Filename)
            if is_python_primitive(v) or isfilename:
                inp_identifier = f"{identifier}_{k}"
                referencedtype = copy.copy(inputs[k].intype) if not isfilename else v
                parsed_type = get_instantiated_type(v)

                if parsed_type and not referencedtype.can_receive_from(parsed_type):
                    raise TypeError(
                        f"The type {parsed_type.id()} inferred from the value '{v}' is not "
                        f"compatible with the '{identifier}.{k}' type: {referencedtype.id()}"
                    )

                referencedtype.optional = True

                indoc = inputs[k].doc
                indoc.quality = InputQualityType.configuration

                v = self.input(
                    inp_identifier,
                    referencedtype,
                    default=v.generated_filename() if isfilename else v,
                    doc=indoc,
                )
            if v is None:
                inp_identifier = f"{identifier}_{k}"
                doc = copy.copy(InputDocumentation.try_parse_from(inputs[k].doc))
                doc.quality = InputQualityType.configuration
                v = self.input(inp_identifier, inputs[k].intype, default=v, doc=doc)

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
        self.has_subworkflow = self.has_subworkflow or isinstance(tool, WorkflowBase)
        self.nodes[identifier] = stp
        self.step_nodes[identifier] = stp

        return stp

    def conditional(
        self, stepid: str, conditions: List[Union[Tuple[Operator, Tool], Tool]]
    ):
        if len(conditions) <= 1:
            raise Exception("A switch statement must include at least 2 conditions")
        # validate tools
        tools = []
        for i in range(len(conditions)):
            con = conditions[i]
            if isinstance(con, tuple):
                tools.append(con[1])
            else:
                if i < (len(conditions) - 1):
                    raise Exception(
                        "A default statement (no condition) must be the last element in the switch"
                    )
                tools.append(con)

        # check output schema
        compare_schema = tools[0].outputs_map()
        in_keys = set(compare_schema.keys())
        non_matching_tools = []
        for i in range(1, len(tools)):
            tool = tools[i]
            outs = tool.outputs_map()
            extra_params = set(outs.keys()) - in_keys
            non_matching_els = list(
                k
                for k, v in compare_schema.items()
                if k not in outs or not v.outtype.can_receive_from(outs[k].outtype)
            )
            if len(non_matching_els) > 0:
                non_matching_tools.append(
                    tool.id() + ": " + ", ".join(non_matching_els + list(extra_params))
                )

        if non_matching_tools:
            raise Exception(
                "Not all tools in the switch statement shared the output schema: "
                + " || ".join(non_matching_tools)
            )

        # should be good to build the workflow:
        w = wrap_steps_in_workflow(stepid, conditions)
        self.step(stepid, w)

    def __getattr__(self, item):
        if item in self.__dict__ or item == "nodes":
            return self.__dict__.get(item)

        try:
            return self.__getitem__(item)
        except KeyError:
            pass

        raise AttributeError(
            f"AttributeError: '{type(self).__name__}' object has no attribute '{item}'"
        )

    def __getitem__(self, item):

        if self.nodes and item in self.nodes:
            node = self.nodes[item]
            if isinstance(node, InputNode):
                return node.as_operator()
            return node

        raise KeyError(f"KeyError: '{type(self).__name__}' object has no node '{item}'")

    @classmethod
    def type(cls) -> ToolType:
        return ToolType.Workflow

    def tool_inputs(self) -> List[TInput]:
        """
        List of ToolInputs of the workflow, we can toss out most of the metadata
        about positioning, prefixes, etc that the ToolInput class uses
        """
        return [
            TInput(i.id(), i.datatype, default=i.default, doc=i.doc)
            for i in self.input_nodes.values()
        ]

    def tool_outputs(self) -> List[TOutput]:
        """
        Similar to inputs, return a list of ToolOutputs of the workflow
        """
        return [
            TOutput(o.id(), o.datatype, doc=o.doc) for o in self.output_nodes.values()
        ]

    def translate(
        self,
        translation: Union[str, SupportedTranslation],
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
        max_duration=None,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        from janis_core.translations import translate_workflow

        return translate_workflow(
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
            max_duration=max_duration,
            allow_empty_container=allow_empty_container,
            container_override=container_override,
        )

    def generate_inputs_override(
        self,
        additional_inputs=None,
        with_resource_overrides=False,
        hints=None,
        include_defaults=True,
        values_to_ignore: Set[str] = None,
        quality_type: List[InputQualityType] = None,
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
            if (
                i.id() in ad
                or i.value
                or not i.datatype.optional
                or (i.default and include_defaults)
            )
            and not (values_to_ignore and i.id() in values_to_ignore)
            and (not (i.doc and quality_type) or i.doc.quality in quality_type)
        }

        if with_resource_overrides:
            from janis_core.translations import CwlTranslator

            d.update(CwlTranslator().build_resources_input(self, hints))

        return d

    def generate_resources_file(
        self,
        translation: Union[str, SupportedTranslation],
        hints: Dict[str, Any] = None,
        to_console=True,
    ):
        from janis_core.translations import build_resources_input

        tr = build_resources_input(self, translation, hints)
        if to_console:
            print(tr)
        return tr

    def get_tools(self) -> Dict[str, CommandTool]:
        tools: Dict[str, CommandTool] = {}
        for t in self.step_nodes.values():
            tl = t.tool
            if isinstance(tl, WorkflowBase):
                tools.update(tl.get_tools())
            elif t.id() not in tools:
                tools[tl.id()] = tl
        return tools

    def get_subworkflows(self):
        tools: Dict[str, WorkflowBase] = {}
        for t in self.step_nodes.values():
            tl = t.tool
            if not isinstance(tl, WorkflowBase):
                continue
            tools[self.versioned_id()] = self
            tools.update(tl.get_subworkflows())

        return tools

    def containers(self) -> Dict[str, str]:
        tools: Dict[str, str] = {}
        for t in self.step_nodes.values():
            tools.update(t.tool.containers())
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

    @staticmethod
    def get_step_ids_from_selector(selector: Selector) -> Set[str]:
        if isinstance(selector, StepOutputSelector):
            return {selector.node.id()}
        elif isinstance(selector, Operator):
            from itertools import chain

            return set(
                chain.from_iterable(
                    Workflow.get_step_ids_from_selector(s)
                    for s in selector.get_leaves()
                )
            )
        return set()

    @staticmethod
    def get_dot_plot_internal(
        tool,
        graph: Optional = None,
        default_base_connection=None,
        prefix="",
        expand_subworkflows=True,
        depth=0,
    ):

        if graph is None:
            from graphviz import Digraph

            graph = Digraph(
                name=tool.id(),
                comment=tool.friendly_name(),
                node_attr={"shape": "record"},
            )

        add_later: Dict[str, Set[str]] = {}

        pref = f"{prefix}_" if prefix else ""

        for stp in tool.step_nodes.values():
            tool = stp.tool

            fn = stp.id()
            if tool.friendly_name():
                fn += f" ({tool.friendly_name()})"
            elif stp.doc and stp.doc.doc:
                fn += f" ({stp.doc.doc})"
            is_subworkflow = isinstance(tool, WorkflowBase)
            if expand_subworkflows and is_subworkflow:
                subid = pref + stp.id()
                bgcolor = f"grey{(9 - depth) * 10}"
                with graph.subgraph(name="cluster_" + subid, comment=stp.doc.doc) as g:
                    g.attr(color=bgcolor, label=fn, style="filled")
                    # if prefix:
                    g.node(subid, shape="Msquare")
                    WorkflowBase.get_dot_plot_internal(
                        tool=tool,
                        graph=g,
                        default_base_connection=subid,
                        prefix=subid,
                        depth=depth + 1,
                    )

            else:
                bgcolor = "grey80" if is_subworkflow else None
                graph.node(
                    pref + stp.id(),
                    fn,
                    style="filled" if is_subworkflow else None,
                    color=bgcolor,
                )

            if stp.sources:
                to_add = set()
                for srcId, steptaginput in stp.sources.items():
                    sti: StepTagInput = steptaginput
                    src = sti.source()
                    if src is None:
                        continue
                    if isinstance(src, list):
                        for s in src:
                            to_add.update(Workflow.get_step_ids_from_selector(s.source))
                    else:

                        to_add.update(Workflow.get_step_ids_from_selector(src.source))
                # if len(to_add) == 0 and default_base_connection is not None:
                #     to_add.add(default_base_connection)
                if to_add:
                    if stp.id() in add_later:
                        add_later[stp.id()].update(to_add)
                    else:
                        add_later[stp.id()] = to_add

        for (src, finals) in add_later.items():
            for f in finals:
                graph.edge(pref + f, pref + src)

        return graph

    def get_dot_plot(
        self,
        show=False,
        log_to_stdout=True,
        expand_subworkflows=False,
        persist_subworkflows=False,
        output_directory: Optional[str] = None
        # these options are primarily for the recu
    ):

        tools = [self]
        if persist_subworkflows:
            tools = [self] + [t for t in self.get_subworkflows().values()]

        Logger.info(f"Generating graphs for {len(tools)} workflows")
        if output_directory:
            output_directory = fully_qualify_filename(output_directory)
            Logger.info(f"Persisting to '{output_directory}'")

        graphs = {}
        for tool in tools:
            graph = self.get_dot_plot_internal(
                tool, expand_subworkflows=expand_subworkflows
            )
            graphs[tool.versioned_id()] = graph

            if output_directory:
                pb = os.path.join(output_directory, tool.versioned_id()) + ".dot"
                Logger.debug(f"Outputting workflow to '{pb}'")
                graph.render(filename=pb, format="png", view=False)

        primary_graph = graphs[self.versioned_id()]
        if log_to_stdout:
            print(primary_graph.source)
        if show:
            primary_graph.render(view=True)

        return graphs

    def version(self):
        meta: WorkflowMetadata = self.bind_metadata() or self.metadata
        if meta and meta.version:
            return meta.version

    def apply_input_documentation(
        self,
        inputs: Dict[str, Union[InputDocumentation, str, Dict[str, any]]],
        should_override=False,
        strict=False,
    ):
        """
        Apply a dictionary of input documentation to a number of input nodes

        :param inputs: Dict[InputNode_ID, Union[InputDocumentation, str, Dict]]
        :param should_override: Should override the doc on an input node.
        :param strict: Ensure every key in the inputs dictionary is in the workflow, otherwise throw an error.
        :return: None
        """
        missing, skipped = set(), set()
        innodes = self.input_nodes
        for inpid, doc in inputs.items():
            if inpid not in innodes:
                if strict:
                    missing.add(inpid)
                continue
            node = innodes[inpid]
            existing_doc = node.doc and node.doc.doc
            if existing_doc is None or should_override:
                node.doc = InputDocumentation.try_parse_from(doc)
            else:
                skipped.add(inpid)

        if missing:
            raise Exception(
                "Couldn't find the following inputs to update: " + ", ".join(missing)
            )

        if skipped:
            Logger.log(
                "Skipped updating fields as they already had documentation: "
                + ", ".join(skipped)
            )


class Workflow(WorkflowBase):
    def __init__(self, **connections):
        super().__init__(**connections)
        # Now that we've initialised everything, we can "construct" the workflows for that subclass this class
        # else, for the WorkflowBuilder it will do nothing and they'll add workflows later
        self.constructor()

    @abstractmethod
    def constructor(self):
        """
        A place to construct your workflows. This is called directly after initialisation.
        :return:
        """
        pass


class DynamicWorkflow(WorkflowBase):
    def __init__(self, **connections):
        super().__init__(**connections)
        self.has_constructed = False

    @abstractmethod
    def constructor(self, inputs: Dict[str, any], hints: Dict[str, str]):
        """
        A place to construct your workflows. This is called after inputs are initially processed

        :param inputs: Dictionary of input values
        :param hints: Dictionary of hints that should be applied
        """
        self.has_constructed = True

    def modify_inputs(self, inputs, hints: Dict[str, str]):
        """
        :param inputs: Dictionary of input values
        :param hints: Dictionary of hints that should be applied
        :return: The modified inputs dictionary
        """
        return inputs


class WorkflowBuilder(Workflow):
    def __init__(
        self,
        identifier: str = None,
        friendly_name: str = None,
        version: str = None,
        metadata: WorkflowMetadata = None,
        tool_provider: str = None,
        tool_module: str = None,
    ):
        self._identifier = identifier
        self._name = friendly_name
        self._version = version
        self._metadata = metadata
        self._tool_provider = tool_provider
        self._tool_module = tool_module

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

    def version(self):
        return self._version

    def tool_provider(self):
        return self._tool_provider

    def tool_module(self):
        return self._tool_module

    def bind_metadata(self):
        return self._metadata

    def __str__(self):
        return f'WorkflowBuilder("{self._identifier}")'


def wrap_steps_in_workflow(
    stepId: str, steps: List[Union[Tuple[Operator, Tool], Tool]]
):
    # skip validation of steps

    # need to resolve connections and try to map them to the same

    w = WorkflowBuilder(stepId)
    workflow_connection_map = {}

    def rebuild_condition(operator: Operator):
        if operator is None:
            return None

        if not isinstance(operator, Selector):
            return operator

        # traverse through the operator, rebuilding and swapping out InputOperator / StepOperator
        # ensure traversal through StringFormatters to looking for the same thing

        if isinstance(operator, StepOutputSelector):
            # pipe in the input of the step_operator
            identifier = f"cond_{operator.node.id()}_{operator.tag}"
            if identifier in w.step_nodes:
                return w[identifier]
            else:
                workflow_connection_map[identifier] = operator
                return w.input(
                    identifier, operator.node.outputs()[operator.tag].outtype
                )
        if isinstance(operator, InputNodeSelector):
            identifier = f"cond_{operator.input_node.id()}"
            if identifier in w.input_nodes:
                return w[identifier]
            else:
                innode: InputNode = operator.input_node
                workflow_connection_map[identifier] = operator
                return w.input(identifier, first_value(innode.outputs()).outtype)

        if isinstance(operator, StringFormatter):
            return StringFormatter(
                operator._format,
                **{k: rebuild_condition(v) for k, v in operator.kwargs.items()},
            )

        if isinstance(operator, Operator):
            return operator.__class__(*[rebuild_condition(t) for t in operator.args])

        # if isinstance(operator, SingleValueOperator):
        #     return operator.__class__(rebuild_condition(operator.internal))
        #
        # if isinstance(operator, TwoValueOperator):
        #     return operator.__class__(
        #         rebuild_condition(operator.lhs), rebuild_condition(operator.rhs)
        #     )

        return operator

    prevconds: List[Operator] = []

    for i in range(len(steps)):
        step = steps[i]
        stepid = "switch_case_" + str(i + 1)

        if isinstance(step, tuple):
            newcond = rebuild_condition(step[0])
            tool = step[1]
            prevcond = or_prev_conds(prevconds)
            cond = (
                AndOperator(newcond, NotOperator(prevcond))
                if prevcond is not None
                else newcond
            )
            prevconds.append(newcond)
        else:
            tool = step
            cond = NotOperator(or_prev_conds(prevconds))

        toolinputs = tool.inputs_map()
        connection_map = {}
        for (k, v) in tool.connections.items():

            # unique
            identifier = stepid + "_" + k

            if is_python_primitive(v):
                # the add step will create the literal for us
                connection_map[k] = v
            else:
                # generate an input for this step, of whatever type our tool wants
                toolin = toolinputs[k]
                workflow_connection_map[identifier] = v
                connection_map[k] = w.input(identifier, toolin.intype, doc=toolin.doc)

        w.step(stepid, tool(**connection_map), when=cond)

    steps = list(w.step_nodes.values())
    # They should all have exact same contract
    outputs = steps[0].tool.outputs_map().values()
    for o in outputs:
        out_source = [stp[o.id()] for stp in steps]
        w.output(o.id(), source=FirstOperator(out_source))

    return w(**workflow_connection_map)
