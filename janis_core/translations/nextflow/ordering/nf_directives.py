


from typing import Type
from abc import ABC, abstractmethod

from ..process.directives import (
    ProcessDirective,
    ContainerDirective,
    DebugDirective,
    CpusDirective,
    DiskDirective,
    MemoryDirective,
    PublishDirDirective,
    TimeDirective
)

directive_priorities: dict[Type[ProcessDirective], int] = {
    DebugDirective: 0,
    ContainerDirective: 1,
    PublishDirDirective: 2,
    CpusDirective: 9,
    DiskDirective: 9,
    MemoryDirective: 9,
    TimeDirective: 9,
}

class DirectiveOrderer(ABC):
    @abstractmethod
    def order(self, directives: list[ProcessDirective]) -> list[ProcessDirective]:
        ...

class PriorityDirectiveOrderer(DirectiveOrderer):
    def order(self, directives: list[ProcessDirective]) -> list[ProcessDirective]:
        return sorted(directives, key=lambda x: directive_priorities[type(x)])

class AlphabeticalDirectiveOrderer(DirectiveOrderer):
    def order(self, directives: list[ProcessDirective]) -> list[ProcessDirective]:
        return sorted(directives, key=lambda x: type(x).__name__)


def order_nf_directives(directives: list[ProcessDirective]) -> list[ProcessDirective]:
    for orderer in [AlphabeticalDirectiveOrderer, PriorityDirectiveOrderer]:
        directives = orderer().order(directives)
    return directives