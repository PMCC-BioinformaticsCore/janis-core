

from typing import Any
from abc import ABC, abstractmethod
from . import rules as rules


class FormattingStrategy(ABC):

    @abstractmethod
    def format(self, starting_text: str, entity: Any) -> str:
        ...

# default
class GenericFormattingStrategy(FormattingStrategy):
    def format(self, starting_text: str, entity: Any) -> str:
        tag = starting_text
        tag = rules.numeric(tag, entity)
        tag = rules.numeric_start(tag, entity)
        tag = rules.non_alphanumeric(tag, entity)
        tag = rules.short_tag(tag, entity)
        # tag = rules.camelify(tag)
        # tag = rules.encode(tag)
        tag = rules.replace_keywords(tag, entity)
        return tag

class WorkflowInputFormattingStrategy(FormattingStrategy):
    def format(self, starting_text: str, entity: Any) -> str:
        tag = starting_text
        if not entity.is_runtime:
            tag = rules.numeric(tag, entity)
            tag = rules.numeric_start(tag, entity)
            tag = rules.non_alphanumeric(tag, entity)
            tag = rules.short_tag(tag, entity)
            # tag = rules.camelify(tag)
            # tag = rules.encode(tag)
            tag = rules.replace_keywords(tag, entity)
        return tag

# capitalisation is allowed
class ToolNameStrategy(FormattingStrategy):
    def format(self, starting_text: str, entity: Any) -> str:
        tag = starting_text
        tag = rules.non_alphanumeric(tag, entity)
        tag = rules.numeric(tag, entity)
        tag = rules.numeric_start(tag, entity)
        tag = rules.short_tag(tag, entity)
        # tag = rules.encode(tag)
        # tag = rules.camelify(tag)
        tag = rules.replace_keywords(tag, entity)
        return tag


STRATEGIES = {
    'Workflow': GenericFormattingStrategy(),
    'WorkflowInput': WorkflowInputFormattingStrategy(),
    'WorkflowStep': GenericFormattingStrategy(),
    'ITool': ToolNameStrategy(),
    'Positional': GenericFormattingStrategy(),
    'Flag': GenericFormattingStrategy(),
    'Option': GenericFormattingStrategy(),
    'RedirectOutput': GenericFormattingStrategy(),
    'WildcardOutput': GenericFormattingStrategy(),
    'InputOutput': GenericFormattingStrategy(),
} 


def format_tag(starting_text: str, entity: Any) -> str:
    strategy = STRATEGIES[entity.__class__.__name__]
    return strategy.format(starting_text, entity)



