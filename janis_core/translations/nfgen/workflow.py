from textwrap import indent
from typing import Optional, Union, List

from janis_core.translations.nfgen.common import NFBase, filter_null


class WorkflowInput(NFBase):
    def __init__(self, name: str, as_param: Optional[str] = None):
        self.name = name
        self.as_param = as_param

    def get_string(self):
        return self.name


class WorkflowOutput(NFBase):
    def __init__(self, name, expression: Optional[str] = None):
        self.name = name
        self.expression = expression

    def get_string(self):
        if self.expression is not None:
            return f"{self.name} = {self.expression}"

        return self.name


class WorkflowPublish(NFBase):
    def __init__(self, name: str, to: str):
        self.name = name
        self.to = to

    def get_string(self):
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

    def prepare_main(self, prefix="  "):
        main = "\n".join(self.main) if isinstance(self.main, list) else self.main

        if self.take or self.emit or self.publish:
            main = "main:\n" + indent(main, 2 * " ")

        return indent(main, prefix)

    def prepare_take(self, prefix="  "):
        if not self.take:
            return None
        return indent(
            "take:\n" + "\n".join(prefix + i.get_string() for i in self.take), "  "
        )

    def prepare_emit(self, prefix="  "):
        if not self.emit:
            return None
        return indent(
            "emit:\n" + "\n".join(prefix + i.get_string() for i in self.emit), "  "
        )

    def prepare_publish(self, prefix="  "):
        if not self.publish:
            return None
        return indent(
            "publish:\n" + "\n".join(prefix + i.get_string() for i in self.publish),
            "  ",
        )

    def get_string(self):
        nl = "\n"

        components = filter_null(
            [
                self.prepare_take(),
                self.prepare_main(),
                self.prepare_emit(),
                self.prepare_publish(),
            ]
        )
        name = self.name or ""
        components_str = (2 * nl).join(components)

        return f"""\
workflow {name} 
{{

{components_str}

}}
"""
