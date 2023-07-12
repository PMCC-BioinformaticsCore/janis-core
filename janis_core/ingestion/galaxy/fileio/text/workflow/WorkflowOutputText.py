


from typing import Tuple

from janis_core.ingestion.galaxy.internal_model.workflow.step.outputs import StepOutput

from ..TextRender import TextRender
from .. import ordering
from .. import formatting

from janis_core.ingestion.galaxy import tags  # remove if possible


class WorkflowOutputText(TextRender):
    def __init__(self, entity: StepOutput):
        super().__init__()
        self.entity = entity

    @property
    def imports(self) -> list[Tuple[str, str]]:
        jtype = self.entity.tool_output.datatype
        imports: list[Tuple[str, str]] = []
        imports.append((jtype.import_path, jtype.classname))

        # TODO opportunity for decorator
        if self.entity.tool_output.array:
            imports.append(('janis_core', 'Array'))

        # TODO opportunity for decorator
        imports = list(set(imports))
        return ordering.order_imports(imports)

    def render(self) -> str:
        step_tag = tags.get(self.entity.step_uuid)
        out_tag = self.entity.tool_output.tag
        datatype_str = formatting.format_typestr(self.entity.tool_output)
        out_str: str = ''
        out_str += 'w.output(\n'
        out_str += f'\t"{step_tag}_{out_tag}",\n'
        out_str += f'\t{datatype_str},\n'
        out_str += f'\tsource=(w.{step_tag}, "{out_tag}")'
        # out_str += f',\n\toutput_folder="{output_folder}"' if output_folder else ''
        # out_str += f',\n\toutput_name="{output_name}"' if output_name else ''
        # out_str += f',\n\textension="{extension}"' if extension else ''
        # out_str += f',\n\tdoc="{doc}"' if doc else ''
        out_str += '\n)\n'
        return out_str