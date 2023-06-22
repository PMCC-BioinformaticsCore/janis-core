

from typing import Optional, Tuple

import re

from .. import utils

from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.expressions.patterns import (
    CHEETAH_SET,
    LINUX_LN,
    LINUX_MV,
    LINUX_CP
)


def resolve_aliases(cmdstr: str) -> str:
    resolver = AliasResolver()
    return resolver.resolve(cmdstr)

def extract_alias(line: str) -> Tuple[Optional[str], Optional[str]]:
    source, dest = None, None
    for func in [
        get_set,
        get_symlink, 
        #get_move, 
        #get_copy, 
    ]:
        source, dest = func(line)
        if source and dest:
            break
    return source, dest

def is_alias_line(line: str) -> bool:
    source, dest = extract_alias(line)
    if source and dest:
        return True
    return False

def get_set(line: str) -> Tuple[Optional[str], Optional[str]]:
    matches = expressions.get_matches(line, CHEETAH_SET)
    if matches:
        m = matches[0]
        return m.group(1), m.group(2)
    return None, None

def get_symlink(line: str) -> Tuple[Optional[str], Optional[str]]:
    matches = expressions.get_matches(line, LINUX_LN)
    if matches:
        m = matches[0]
        return m.group(2), m.group(1)
    return None, None

def get_copy(line: str) -> Tuple[Optional[str], Optional[str]]:
    matches = expressions.get_matches(line, LINUX_CP)
    if matches:
        m = matches[0]
        return m.group(2), m.group(1)
    return None, None

def get_move(line: str) -> Tuple[Optional[str], Optional[str]]:
    matches = expressions.get_matches(line, LINUX_MV)
    if matches:
        m = matches[0]
        return m.group(1), m.group(2)
    return None, None




# class because it may change in future. 
# want to keep the .get() and .add() methods stable
class AliasRegister:
    def __init__(self):
        self.aliases: dict[str, list[str]] = {}
    
    def get(self, source: str) -> Optional[str]:
        if source in self.aliases:
            destinations = self.aliases[source]
            if len(destinations) == 1:  # won't return if ambiguity regarding dest
                return destinations[0]
        return None
    
    def add(self, source: str, dest: str) -> None:
        if source not in self.aliases:
            self.aliases[source] = []    
        self.aliases[source].append(dest)


class AliasResolver:
    def __init__(self):
        self.register: AliasRegister = AliasRegister()

    def resolve(self, cmdstr: str) -> str:
        self.register_aliases(cmdstr)
        return self.resolve_aliases(cmdstr)

    def register_aliases(self, cmdstr: str) -> None:
        lines = utils.split_lines(cmdstr)
        for line in lines:
            self.register_line(line)

    def register_line(self, line: str) -> None:
        source, dest = extract_alias(line)
        if source and dest and source != dest:
            self.register.add(source, dest)

    def resolve_aliases(self, cmdstr: str) -> str:
        resolved_lines: list[str] = []
        for line in utils.split_lines_whitespace(cmdstr):
            line = self.resolve_line(line)
            resolved_lines.append(line)
        cmdstr = utils.join_lines(resolved_lines)
        return cmdstr

    def resolve_line(self, line: str) -> str:
        if not is_alias_line(line):
            for source in self.register.aliases:
                dest = self.register.get(source)
                if dest:
                    # ensure ends/starts with whitespace or quotes or forwardslash or is curly brackets
                    matches = self.get_line_matches(source, line)
                    for match in matches:
                        line = self.perform_replacement(match, dest, line)
        return line

    def get_line_matches(self, source: str, line: str) -> list[re.Match[str]]:
        pattern = re.escape(source)
        raw_matches = re.finditer(pattern, line)
        delims = ['"', "'", ' ', '}']
        matches: list[re.Match[str]] = []

        for m in raw_matches:
            if m.end() == len(line):
                matches.append(m)
            elif line[m.end()] in delims:
                matches.append(m)
            elif m.end()+1 < len(line) and line[m.end()+1] in delims:
                matches.append(m)
        
        return matches

    def perform_replacement(self, match: re.Match[str], dest: str, line: str) -> str:
        line = line[:match.start()] + dest + line[match.end():]
        return line
        
        