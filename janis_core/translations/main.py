

from typing import Optional, Any
from inspect import isclass

from janis_core import settings
from janis_core import CodeTool, CommandToolBuilder, WorkflowBase, WorkflowBuilder
from janis_core import Tool
from janis_core.utils import lowercase_dictkeys
from janis_core.translation_deps.supportedtranslations import SupportedTranslation
from janis_core.translations.common import to_builders
from janis_core.translations.common import prune_workflow
from janis_core.translations.common import wrap_tool_in_workflow
from .translationbase import TranslatorBase


def translate(
    entity: Tool,
    dest_fmt: str,
    mode: Optional[str] = None,

    # file io
    to_disk: Optional[bool] = None,
    export_path: Optional[str] = None,
    as_workflow: Optional[bool] = None,
    should_zip: Optional[bool] = None,   
    to_console: Optional[bool] = None,
    tool_to_console: Optional[bool] = None,
    write_inputs_file: Optional[bool] = None,
    source_files: Optional[list[str]] = None,
    
    # inputs
    additional_inputs: Optional[dict[str, str]] = None,
    hints: Optional[dict[str, str]] = None,
    
    # containers
    with_container: Optional[bool] = None,
    allow_empty_container: Optional[bool] = None,
    container_override: Optional[str | dict[str, Any]] = None,
    
    # resouces
    with_resource_overrides: Optional[bool] = None,
    merge_resources: Optional[bool] = None,
    max_cores: Optional[int] = None,
    max_mem: Optional[int] = None,
    max_duration: Optional[int] = None,
    
    # misc
    render_comments: Optional[bool] = None,
    should_validate: Optional[bool] = None,

) -> Any:  
    
    # settings 
    settings.translate.DEST = dest_fmt             # set translate dest
    settings.validation.STRICT_IDENTIFIERS = False
    settings.validation.VALIDATE_STRINGFORMATTERS = False
    
    if mode is not None:
        settings.translate.MODE = mode
    if export_path:
        settings.translate.EXPORT_PATH = export_path
        settings.translate.TO_DISK = True
    if to_disk:
        settings.translate.TO_DISK = True
    if allow_empty_container is not None:
        settings.translate.ALLOW_EMPTY_CONTAINER = allow_empty_container
    if merge_resources is not None:
        settings.translate.MERGE_RESOURCES = merge_resources
    if render_comments is not None:
        settings.translate.RENDER_COMMENTS = render_comments
    if should_validate is not None:
        settings.translate.SHOULD_VALIDATE = should_validate
    if should_zip is not None:
        settings.translate.SHOULD_ZIP = should_zip
    if to_console is not None:
        settings.translate.TO_CONSOLE = to_console
    if tool_to_console is not None:
        settings.translate.TOOL_TO_CONSOLE = tool_to_console
    if with_container is not None:
        settings.translate.WITH_CONTAINER = with_container
    if with_resource_overrides is not None:
        settings.translate.WITH_RESOURCE_OVERRIDES = with_resource_overrides
    if write_inputs_file is not None:
        settings.translate.WRITE_INPUTS_FILE = write_inputs_file
    if source_files is not None:
        settings.translate.SOURCE_FILES = source_files
    if additional_inputs is not None:
        settings.translate.ADDITIONAL_INPUTS = additional_inputs
    if container_override is not None:
        if isinstance(container_override, str):
            container_override = {entity.id().lower(): container_override}
        settings.translate.CONTAINER_OVERRIDE = lowercase_dictkeys(container_override)
    if hints is not None:
        settings.translate.HINTS = hints
    if max_cores is not None:
        settings.translate.MAX_CORES = max_cores
    if max_duration is not None:
        settings.translate.MAX_DURATION = max_duration
    if max_mem is not None:
        settings.translate.MAX_MEM = max_mem

    # preprocessing
    entity = to_builders(entity)
    if settings.translate.MODE in ['skeleton', 'regular'] and isinstance(entity, WorkflowBuilder):
        assert(isinstance(entity, WorkflowBuilder))
        prune_workflow(entity)

    if as_workflow and isinstance(entity, CommandToolBuilder):
        entity = wrap_tool_in_workflow(entity)

    # select the translation unit 
    translator = get_translator(dest_fmt)

    # do translation 
    if isinstance(entity, WorkflowBuilder):
        return translator.translate_workflow(entity)
    elif isinstance(entity, CommandToolBuilder):
        return translator.translate_tool(entity)
    elif isinstance(entity, CodeTool):
        return translator.translate_code_tool(entity)
    else:
        name = entity.__name__ if isclass(entity) else entity.__class__.__name__
        raise Exception("Unsupported tool type: " + name)

def get_translator(translation: str | SupportedTranslation) -> TranslatorBase:
    if not isinstance(translation, SupportedTranslation):
        translation = SupportedTranslation(translation)
    return translation.get_translator()


# standalone funcs for just building inputs, rather than doing full translation. 
# doesn't really work or make sense. 
def build_resources_input(
    workflow: WorkflowBase,
    translation: str | SupportedTranslation,
    hints: Optional[dict[str, str]] = None,
    max_cores: Optional[int] = None,
    max_mem: Optional[int] = None,
    max_duration: Optional[int] = None,
) -> dict[str, Any]:
    if hints is not None:
        settings.translate.HINTS = hints
    if max_cores is not None:
        settings.translate.MAX_CORES = max_cores
    if max_duration is not None:
        settings.translate.MAX_DURATION = max_duration
    if max_mem is not None:
        settings.translate.MAX_MEM = max_mem

    translator = get_translator(translation)
    inputs_dict = translator._build_resources_dict(workflow)
    return inputs_dict

def build_resources_file(
    workflow: WorkflowBase,
    translation: str | SupportedTranslation,
    hints: Optional[dict[str, str]] = None,
    max_cores: Optional[int] = None,
    max_mem: Optional[int] = None,
    max_duration: Optional[int] = None,
) -> str:
    translator = get_translator(translation)
    inputs_dict = build_resources_input(
        workflow,
        translation,
        hints=hints,
        max_cores=max_cores,
        max_mem=max_mem,
        max_duration=max_duration
    )
    inputs_file = translator.stringify_translated_inputs(inputs_dict)
    return inputs_file