from textwrap import indent
from typing import Optional, Union, List

from janis_core.translations.nfgen.common import NFBase, filter_null
from . import settings

class WorkflowInput(NFBase):
    def __init__(self, name: str, as_param: Optional[str] = None):
        self.name = name
        self.as_param = as_param

    def get_string(self) -> str:
        return self.name


class WorkflowOutput(NFBase):
    def __init__(self, name, expression: Optional[str] = None):
        self.name = name
        self.expression = expression

    def get_string(self) -> str:
        if self.expression is not None:
            return f"{self.name} = {self.expression}"

        return self.name


class WorkflowPublish(NFBase):
    def __init__(self, name: str, to: str):
        self.name = name
        self.to = to

    def get_string(self) -> str:
        return f"{self.name} to: {self.to}"


class Workflow(NFBase):
    def __init__(
        self,
        name: Optional[str],
        main: Optional[Union[str, List[str]]],
        take: Optional[List[WorkflowInput]] = None,
        emit: Optional[List[WorkflowOutput]] = None,
        publish: Optional[List[WorkflowPublish]] = None,
    ):
        self.name = name
        self.main = main
        self.take = take or []
        self.emit = emit or []
        self.publish = publish or []

    @property
    def inputs(self):
        return self.take

    def prepare_main(self):
        main = "\n".join(self.main) if isinstance(self.main, list) else self.main

        if self.take or self.emit or self.publish:
            main = "main:\n" + indent(main, settings.NEXTFLOW_INDENT)

        return indent(main, settings.NEXTFLOW_INDENT)

    def prepare_take(self):
        if not self.take:
            return None
        return indent(
            "take:\n" + "\n".join(prefix + i.get_string() for i in self.take), 
            settings.NEXTFLOW_INDENT
        )

    def prepare_emit(self):
        if not self.emit:
            return None
        return indent(
            "emit:\n" + "\n".join(prefix + i.get_string() for i in self.emit), 
            settings.NEXTFLOW_INDENT
        )

    def prepare_publish(self):
        if not self.publish:
            return None
        return indent(
            "publish:\n" + "\n".join(prefix + i.get_string() for i in self.publish),
            settings.NEXTFLOW_INDENT
        )

    def get_string(self) -> str:
        components = filter_null(
            [
                self.prepare_take(),
                self.prepare_main(),
                self.prepare_emit(),
                self.prepare_publish(),
            ]
        )
        name = self.name or ""
        components_str = '\n'.join(components)

        return f"""\
workflow {name} {{

{components_str}

}}
"""
