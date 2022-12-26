from textwrap import indent
from typing import Optional

from . import settings
from . import naming


def filter_null(iterable):
    if iterable is None:
        return None
    if not isinstance(iterable, list):
        raise Exception(
            f"filter_null: Expected iterable ('{str(iterable)}') to be of type list, received '{type(iterable)}'"
        )

    return [el for el in iterable if el is not None]


class WorkflowTake:
    """
    A nextflow workflow can accept channels as inputs.
    These channels are assigned new names in the workflow `take:` section.
    The only thing binding these channels is their argument position. 
    eg 
        // main wf
        workflow {
            my_pipeline(Channel.of(params.greeting))
        }
        // subwf
        workflow my_pipeline {
            take:
            ch_greeting
        }
    """
    def __init__(self, name: str, as_param: Optional[str] = None):
        self.name = name
        self.as_param = as_param

    def get_string(self) -> str:
        return self.name


class WorkflowEmit:
    """
    A nextflow workflow can broadcast channels as outputs. 
    For translation, we always use named outputs.
    eg
        // main wf
        workflow {
            my_pipeline(Channel.of(params.greeting))
            my_pipeline.out.mydata.view()
        }
        // sub wf
        workflow my_pipeline {
            take:
            ch_greeting

            main:
            CONVERTTOUPPER(ch_greeting)

            emit:
            mydata = CONVERTTOUPPER.out.upper
        }
    """
    def __init__(self, name: str, expression: Optional[str] = None):
        self.name = name
        self.expression = expression

    def get_string(self) -> str:
        if self.expression is not None:
            return f"{self.name} = {self.expression}"

        return self.name


# class WorkflowPublish:
#     def __init__(self, name: str, to: str):
#         self.name = name
#         self.to = to

#     def get_string(self) -> str:
#         return f"{self.name} to: {self.to}"


class Workflow:
    def __init__(
        self,
        name: str,
        main: list[str],
        take: Optional[list[WorkflowTake]]=None,
        emit: Optional[list[WorkflowEmit]]=None,
        is_subworkflow: bool=False
    ):
        self.name = naming.gen_varname_workflow(name)
        self.main = main
        self.take = take or []
        self.emit = emit or []
        self.is_subworkflow = is_subworkflow

    @property
    def main_block(self) -> Optional[str]:
        main = "\n".join(self.main)
        if self.is_subworkflow:
            main = "main:\n" + main
        return indent(main, settings.NF_INDENT)

    @property
    def take_block(self) -> Optional[str]:
        if not self.take:
            return None
        return indent(
            "take:\n" + "\n".join(i.get_string() for i in self.take) + '\n', 
            settings.NF_INDENT
        )

    @property
    def emit_block(self) -> Optional[str]:
        if not self.emit:
            return None
        return indent(
            "emit:\n" + "\n".join(i.get_string() for i in self.emit),
            settings.NF_INDENT
        )

    def get_string(self) -> str:
        components = filter_null(
            [
                self.take_block,
                self.main_block,
                self.emit_block,
            ]
        )
        components_str = '\n'.join(components)

        return f"""\
workflow {self.name} {{

{components_str}

}}
"""
