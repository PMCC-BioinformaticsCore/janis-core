"""
This set of code is used for building ONE-WAY transformations between types.
We can use this to build up a set of operations to
"""
from typing import Optional, List, Dict

from janis_core.tool.tool import Tool
from janis_core.types import get_instantiated_type, ParseableType, DataType, File
from janis_core.utils.logger import Logger
from janis_core.workflow.workflow import Workflow, WorkflowBuilder


class JanisTransformation:
    def id(self):
        return f"{self.type1}>{self.type2}"

    def __init__(
        self,
        start_type: ParseableType,
        finish_type: ParseableType,
        tool: Tool,
        relevant_tool_input: Optional[str] = None,
        relevant_tool_output: Optional[str] = None,
    ):
        self.type1 = get_instantiated_type(start_type)
        self.type2 = get_instantiated_type(finish_type)
        self.tool = tool

        connection_type = f"`{self.type1} -> {self.type2}`"

        Logger.log(
            f"Building transformation for {connection_type} using tool '{tool.id()}"
        )

        self.relevant_tool_input = self.evaluate_tool_input(relevant_tool_input)
        self.relevant_tool_output = self.evaluate_tool_output(relevant_tool_output)

    def evaluate_tool_input(self, relevant_tool_input: Optional[str]):

        connection_type = f"`{self.type1} -> {self.type2}`"

        ins = self.tool.tool_inputs()
        ins_keys = set(i.id() for i in ins)

        if not relevant_tool_input:
            relevant_types = [t for t in ins if isinstance(t.intype, type(self.type1))]
            if len(relevant_types) == 0:
                raise Exception(
                    f"Couldn't detect an input of '{self.tool.id()}' of "
                    f"type '{self.type1}' for {connection_type} JanisTransformation"
                )
            elif len(relevant_types) > 1:
                raise Exception(
                    f"There were too many relevant inputs of '{self.tool.id()}' of type '{self.type1}' for "
                    f"JanisTransformation to automatically select the relevant input, you should "
                    f"provide a 'relevant_tool_input' in the {connection_type} JanisTransformation."
                )
            relevant_tool_input = relevant_types[0].id()

        if relevant_tool_input not in ins_keys:
            raise Exception(
                f"Couldn't find input '{relevant_tool_input}' in '{self.tool.id()}' for "
                f"{connection_type} JanisTransformation"
            )

        return relevant_tool_input

    def evaluate_tool_output(self, relevant_tool_output: Optional[str]):

        connection_type = f"`{self.type1} -> {self.type2}`"

        outs = self.tool.tool_outputs()
        outs_keys = set(o.id() for o in outs)

        if not relevant_tool_output:
            relevant_types = [
                t for t in outs if isinstance(t.outtype, type(self.type2))
            ]
            if len(relevant_types) == 0:
                raise Exception(
                    f"Couldn't detect an output of '{repr(self.tool)}' of "
                    f"type '{self.type2}' for {connection_type} JanisTransformation"
                )
            elif len(relevant_types) > 1:
                raise Exception(
                    f"There were too many outputs of '{repr(self.tool)}' of type '{self.type2}' for "
                    f"JanisTransformation to automatically select the relevant input, you should "
                    f"provide a 'relevant_tool_output' in the {connection_type} JanisTransformation."
                )
            relevant_tool_output = relevant_types[0].id()

        if relevant_tool_output not in outs_keys:
            raise Exception(
                f"Couldn't find output '{relevant_tool_output}' in '{self.tool.id()}' for "
                f"{connection_type} JanisTransformation"
            )

        return relevant_tool_output

    @staticmethod
    def convert_transformations_to_workflow(transformations) -> Workflow:

        transformations: List[JanisTransformation] = transformations

        initial_tr: JanisTransformation = transformations[0]
        final_tr: JanisTransformation = transformations[-1]

        w = WorkflowBuilder(
            f"convert_{initial_tr.type1.name().lower()}_to_{final_tr.type2.name().lower()}"
        )

        prev_input = w.input(f"inp_{initial_tr.type1.name().lower()}", initial_tr.type1)

        for transform in transformations:
            stpid = f"transform_{transform.type1.name().lower()}_to_{transform.type2.name().lower()}"
            stp_inputs = {
                **(transform.tool.connections or {}),
                transform.relevant_tool_input: prev_input,
            }

            prev_input = w.step(stpid, transform.tool(**stp_inputs))[
                transform.relevant_tool_output
            ]

        w.output("out", final_tr.type2, source=prev_input, output_name=False)

        return w


class JanisTransformationGraph:
    def __init__(self):

        self._edges: Dict[str, List[JanisTransformation]] = {}

    def build_workflow_to_translate(
        self, source_dt: ParseableType, desired_dt: ParseableType
    ) -> Optional[Workflow]:
        transformations = self.find_connection(source_dt, desired_dt)

        if len(transformations) == 0:
            return None

        return JanisTransformation.convert_transformations_to_workflow(transformations)

    def add_edges(self, edges: List[JanisTransformation]):
        for edge in edges:
            dt_id = edge.type1.name()
            if dt_id in self._edges:
                self._edges[dt_id].append(edge)
            else:
                self._edges[dt_id] = [edge]

    def find_connection(
        self, source_dt: ParseableType, desired_dt: ParseableType
    ) -> List[JanisTransformation]:

        from inspect import getmro

        source = get_instantiated_type(source_dt)
        desired = get_instantiated_type(desired_dt)

        if desired.can_receive_from(source):
            return []

        types = getmro(type(source))

        for T in types:
            if not issubclass(T, DataType) or T == DataType:
                continue

            transformation = self.find_connection_inner(T, desired)
            if transformation is not None:
                return transformation

        raise Exception(
            f"There's no transformation that can satisfy {source.name()} -> {desired.name()}"
        )

    def find_connection_inner(
        self, source_dt: ParseableType, desired_dt: ParseableType
    ) -> Optional[List[JanisTransformation]]:

        source = get_instantiated_type(source_dt)
        desired = get_instantiated_type(desired_dt)

        if desired.can_receive_from(source):
            return []

        queue: List[JanisTransformation] = []
        parent_mapping: Dict[str, JanisTransformation] = {}

        desired_dt_name = desired.name()

        queue.extend(self._edges.get(source.name(), []))

        parent_mapping[source.name()] = None

        while len(queue) > 0:
            edge = queue.pop(0)

            end_name = edge.type2.name()

            if end_name in parent_mapping:
                # we've already got a parent (hence we've already seen it), so let's skip it
                continue

            parent_mapping[end_name] = edge

            if end_name == desired_dt_name:
                return JanisTransformationGraph.trace(parent_mapping, end_name)

            queue.extend(self._edges.get(end_name, []))

        return None

    @staticmethod
    def trace(
        parent_mapping: Dict[str, JanisTransformation], start: str
    ) -> List[JanisTransformation]:
        pathway = []
        key = start
        while parent_mapping[key] is not None:
            edge = parent_mapping[key]
            key = edge.type1.name()
            pathway.insert(0, edge)

        return pathway
