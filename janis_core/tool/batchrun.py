from inspect import isclass
from typing import Union, List, Type

from janis_core import WorkflowBase
from janis_core.utils.scatter import ScatterDescription, ScatterMethod

from janis_core.types import Array

from janis_core.workflow.workflow import Tool, Workflow


class BatchRunTool(Workflow):
    def __init__(
        self,
        innertool: Union[Tool, Type[Tool]],
        fields: List[str],
        group_by: str = None,
    ):
        self.inner: Tool = innertool() if isclass(innertool) else innertool
        self.fields = set(fields)
        self.group_by = group_by

        super().__init__()

    def friendly_name(self):
        return "Multiplied " + self.inner.friendly_name()

    def constructor(self):
        ins = self.inner.tool_inputs()

        inkeys = set(i.id() for i in ins)
        invalid_keys = self.fields - inkeys
        if len(invalid_keys) > 0:
            raise Exception(
                f"Couldn't create BatchRunTool from fields {', '.join(invalid_keys)} "
                f"as they do not exist on '{self.inner.id()}'"
            )

        innode_map = {}

        for i in ins:
            intype = i.intype
            default = i.default
            if i.id() in self.fields:
                intype = Array(intype)
                default = None
            innode_map[i.id()] = self.input(i.id(), intype, default=default, doc=i.doc)

        self.step(
            self.inner.id(),
            self.inner(**innode_map),
            scatter=ScatterDescription(list(self.fields), ScatterMethod.dot),
        )

        if isinstance(self.inner, WorkflowBase):
            # We can do special output_folder grouping
            for oid, o in self.output_nodes.items():
                folders = o.output_folder
                if self.group_by:
                    if not folders:
                        folders = []
                    folders.append(self[self.group_by])

                self.output(
                    oid, o.datatype, output_folder=folders, output_name=o.output_name
                )
        else:
            for o in self.inner.tool_outputs():
                f = self[self.group_by] if o else None
                self.output(o.id(), o.outtype, output_folder=f)

    def id(self) -> str:
        return "multiple_" + self.inner.id()
