from typing import List, Type, Union, Optional, Dict, Tuple

from janis_core.graph.node import Node, NodeTypes
from janis_core.graph.stepinput import StepInput
from janis_core.translations import ExportPathKeywords
from janis_core.types import DataType, ParseableType, get_instantiated_type
from janis_core.tool.tool import Tool, ToolType, ToolTypes, ToolInput, ToolOutput
import janis_core.translations as translations

ConnectionSource = Union[Node, Tuple[Node, str]]


def verify_or_try_get_source(source: ConnectionSource):
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
    def __init__(self, wf, identifier: str, datatype: DataType, default: any):
        super().__init__(wf, NodeTypes.INPUT, identifier)
        self.datatype = datatype
        self.default = default

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {"": ToolOutput(self.identifier, self.datatype)}

    def inputs(self):
        return None


class StepNode(Node):
    def __init__(self, wf, identifier, tool: Tool):
        super().__init__(wf, NodeTypes.STEP, identifier)
        self.tool = tool
        self.sources: Dict[str, StepInput] = {}

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
        print("Searching for " + item)
        if item in self.__dict__:
            return self.__dict__[item]
        print("Not found in __dict__: " + item)

        if "tool" in self.__dict__:
            ins = self.inputs()
            if item in ins:
                return self, item
            outs = self.outputs()
            if item in outs:
                return self, item

            tags = ", ".join(
                [f"in.{i}" for i in ins.keys()] + [f"out.{o}" for o in outs.keys()]
            )

            raise AttributeError(
                f"Step '{self.id()}' with tool '{self.tool.id()}' has no identifier '{item}' ({tags})"
            )
        raise AttributeError(
            f"Step '{self.id()}' with tool '{self.tool.id()}' has no identifier '{item}'"
        )

    # "always set attributes keys"
    always_set = {"tool", "sources", "node_type", "identifier"}

    def __setitem__(self, key, value):
        print("Attempting to set: " + key)
        if key not in StepNode.always_set and key in self.inputs():
            return self._add_edge(key, value)

        self.__dict__[key] = value


class OutputNode(Node):
    def __init__(
        self, wf, identifier: str, datatype: DataType, source: ConnectionSource
    ):
        super().__init__(wf, NodeTypes.INPUT, identifier)
        self.datatype = datatype
        self.source = verify_or_try_get_source(source)

    def inputs(self) -> Dict[str, ToolInput]:
        # Program will just grab first value anyway
        return {"": ToolInput(self.identifier, self.datatype)}

    def outputs(self):
        return None


class Workflow2(Tool):
    def __init__(self, identifier: str):
        super().__init__()
        self.identifier = identifier

        # { nodeId: Node }
        self.nodes = {}
        # { stepId: StepInput }
        self.connections = {}

    def input(self, identifier: str, datatype: ParseableType, default: any = None):
        """
        Create an input node on a workflow
        :return:
        """

        if identifier in self.nodes:
            existing = self.nodes[identifier]
            raise Exception(
                f"There already exists a node (and component) with id '{identifier}'. The added "
                f"component ('{repr(datatype)}') clashes with '{repr(existing)}')."
            )

        inp = InputNode(
            self,
            identifier=identifier,
            datatype=get_instantiated_type(datatype),
            default=default,
        )
        self.nodes[identifier] = inp
        return inp

    def output(
        self,
        identifier: str,
        datatype: Optional[ParseableType] = None,
        source: Union[StepNode, ConnectionSource] = None,
    ):
        if identifier in self.nodes:
            existing = self.nodes[identifier]
            raise Exception(
                f"There already exists a node (and component) with id '{identifier}'. The added "
                f"component ('{repr(datatype)}') clashes with '{repr(existing)}')."
            )

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
        return otp

    def step(self, identifier: str, tool: Union[Tool, Type[Tool]], **connections):
        if identifier in self.nodes:
            existing = self.nodes[identifier]
            raise Exception(
                f"There already exists a node (and component) with id '{identifier}'. The added "
                f"component ('{repr(tool)}') clashes with '{repr(existing)}')."
            )

        if issubclass(tool, Tool):
            tool = tool()

        tool.workflow = self
        inputs = tool.inputs_map()

        provided_keys = set(connections.keys())
        required_keys = set(i for i, v in inputs.items() if not v.input_type.optional)

        if not required_keys.issubset(provided_keys):
            missing = ", ".join(required_keys - provided_keys)
            raise Exception(
                f"Missing the parameters {missing} when creating '{identifier}' ({tool.id()})"
            )

        stp = StepNode(self, identifier=identifier, tool=tool)

        for (k, v) in connections.items():
            stp._add_edge(k, verify_or_try_get_source(v))

        self.nodes[identifier] = stp
        return stp

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__.get(item)

        if item in self.nodes:
            return self.nodes[item]

        raise AttributeError(
            f"AttributeError: '{type(self).__name__}' object has no attribute '{item}'"
        )

    @classmethod
    def type(cls) -> ToolType:
        return ToolTypes.Workflow

    def id(self) -> str:
        return self.identifier

    def inputs(self) -> List[ToolInput]:
        pass

    def outputs(self) -> List[ToolOutput]:
        pass

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
        should_validate=False,
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
            should_validate=should_validate,
            merge_resources=merge_resources,
            hints=hints,
            allow_null_if_not_optional=allow_null_if_not_optional,
            additional_inputs=additional_inputs,
            max_cores=max_cores,
            max_mem=max_mem,
        )


if __name__ == "__main__":
    from janis_unix import Echo

    w = Workflow2("simple")
    w.input("inputId", str)
    w.step("echo", Echo, inp=w.inputId)
    w.output("out", source=w.echo)
    w.translate("wdl")
