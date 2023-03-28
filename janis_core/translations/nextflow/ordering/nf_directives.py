


from typing import Type
from abc import ABC, abstractmethod

from ..model.process.directives import (
    NFProcessDirective,
    NFContainerDirective,
    NFDebugDirective,
    NFCpusDirective,
    NFDiskDirective,
    NFMemoryDirective,
    NFPublishDirDirective,
    NFTimeDirective
)

directive_priorities: dict[Type[NFProcessDirective], int] = {
    NFDebugDirective: 0,
    NFContainerDirective: 1,
    NFPublishDirDirective: 2,
    NFCpusDirective: 9,
    NFDiskDirective: 9,
    NFMemoryDirective: 9,
    NFTimeDirective: 9,
}

class DirectiveOrderer(ABC):
    @abstractmethod
    def order(self, directives: list[NFProcessDirective]) -> list[NFProcessDirective]:
        ...

class PriorityDirectiveOrderer(DirectiveOrderer):
    def order(self, directives: list[NFProcessDirective]) -> list[NFProcessDirective]:
        return sorted(directives, key=lambda x: directive_priorities[type(x)])

class AlphabeticalDirectiveOrderer(DirectiveOrderer):
    def order(self, directives: list[NFProcessDirective]) -> list[NFProcessDirective]:
        return sorted(directives, key=lambda x: type(x).__name__)


def order_nf_directives(directives: list[NFProcessDirective]) -> list[NFProcessDirective]:
    for orderer in [AlphabeticalDirectiveOrderer, PriorityDirectiveOrderer]:
        directives = orderer().order(directives)
    return directives