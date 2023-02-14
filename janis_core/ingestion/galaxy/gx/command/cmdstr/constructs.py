

from janis_core.ingestion.galaxy.logs import logging
from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional, Tuple
from uuid import uuid4


CH_OPEN_CONDITIONAL = set(['#if ', '#unless '])
CH_WITHIN_CONDITIONAL = set(['#else if ', '#elif ', '#else'])
CH_CLOSE_CONDITIONAL = set(['#end if', '#end unless'])
CH_ALL_CONDITIONAL = CH_OPEN_CONDITIONAL | CH_WITHIN_CONDITIONAL | CH_CLOSE_CONDITIONAL

CH_OPEN_LOOP = set(['#for ', '#while '])
CH_CLOSE_LOOP = set(['#end for', '#end while'])
CH_ALL_LOOP = CH_OPEN_LOOP | CH_CLOSE_LOOP

CH_TRY_EXCEPT = set(['#try ', '#except', '#end try'])

CH_OPEN_FUNC = set(['#def '])
CH_CLOSE_FUNC = set(['#end def'])
CH_ALL_FUNC = CH_OPEN_FUNC | CH_CLOSE_FUNC

CH_ENV = set(['#set ', '#import ', '#from ', '#silent ', '#echo '])

LINUX_ALIAS = set(['set ', 'ln ', 'cp ', 'mv ', 'export '])
LINUX_CMD = set(['mkdir ', 'tar ', 'ls ', 'head ', 'wget ', 'grep ', 'awk ', 'cut ', 'sed ', 'gzip ', 'gunzip ', 'cd ', 'echo ', 'trap ', 'touch '])

CH_CONSTRUCTS = CH_ALL_CONDITIONAL | CH_ALL_LOOP | CH_ALL_FUNC | CH_TRY_EXCEPT | CH_ENV


"""small module which holds <command> section constructs of importance"""


class ConstructType(Enum):
    CONDITIONAL = auto()
    LOOP        = auto()
    FUNCTION    = auto()


@dataclass
class Construct:
    subtype: ConstructType
    
    def __post_init__(self):
        self.lines: list[str] = []
        self.uuid: str = str(uuid4())

# construct subtypes and their opening / closing cheetah tags (for identification in text)
constructs_opens_closes: list[Tuple[ConstructType, set[str], set[str]]] = [
    (ConstructType.CONDITIONAL, CH_OPEN_CONDITIONAL, CH_CLOSE_CONDITIONAL),
    (ConstructType.LOOP, CH_OPEN_LOOP, CH_CLOSE_LOOP),
    (ConstructType.FUNCTION, CH_OPEN_FUNC, CH_CLOSE_FUNC),
]


class ConstructStack:
    def __init__(self):
        self.stack: list[Construct] = []
    
    def add(self, construct: Construct) -> None:
        self.stack.append(construct)
        if construct.subtype == ConstructType.LOOP:
            logging.has_cheetah_loop()
        elif construct.subtype == ConstructType.FUNCTION:
            logging.has_cheetah_function()
    
    def pop(self, subtype: ConstructType) -> None:
        if self.depth > 0:
            assert(self.current_construct)
            assert(self.current_construct.subtype == subtype)
            self.stack.pop()

    @property
    def depth(self) -> int:
        return len(self.stack)
    
    @property
    def current_construct(self) -> Optional[Construct]:
        if len(self.stack) > 0:
            return self.stack[-1]
        return None

    def within(self, ctype: ConstructType) -> bool:
        for construct in self.stack:
            if construct.subtype == ctype:
                return True
        return False


class ConstructTracker:

    def __init__(self):
        self.stack: ConstructStack = ConstructStack()
    
    # methods to manage construct stack
    def update(self, line: str):
        self.handle_entering_constructs(line)
        self.handle_line(line)
        self.handle_leaving_constructs(line)

    def handle_entering_constructs(self, line: str) -> None:
        for subtype, openings, _ in constructs_opens_closes:
            if any ([line.startswith(openstr) for openstr in openings]):
                construct = Construct(subtype)
                self.stack.add(construct)

    def handle_line(self, line: str) -> None:
        if self.stack.current_construct: 
            self.stack.current_construct.lines.append(line)

    def handle_leaving_constructs(self, line: str) -> None:
        for subtype, _, closings in constructs_opens_closes:
            if any ([line.startswith(closestr) for closestr in closings]):
                self.stack.pop(subtype)

    @property
    def within_conditional(self) -> bool:
        if self.stack.within(ConstructType.CONDITIONAL):
            return True
        return False
    
    @property
    def within_loop(self) -> bool:
        if self.stack.within(ConstructType.LOOP):
            return True
        return False
    
    @property
    def within_banned_segment(self) -> bool:
        banned_constructs = [ConstructType.LOOP, ConstructType.FUNCTION]
        for construct in banned_constructs: 
            if self.stack.within(construct):
                return True
        return False

    def active_is_boundary(self, line: str) -> bool:
        if any ([line.startswith(kw) for kw in CH_CONSTRUCTS]):
            return True
        return False

