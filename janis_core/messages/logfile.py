
from typing import Optional
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class LogLine:
    level: str
    category: Optional[str]
    uuid: Optional[str]
    message: str


class LogFile:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.lines: list[LogLine] = []
        self.load()

    def add(self, level: str, category: Optional[str], uuid: Optional[str], message: str) -> None:
        # add to file
        message = message.replace('\n', '')
        message = fr'{message}'
        with open(self.filepath, 'a') as fp:
            fp.write(f'{level}\t{category}\t{uuid}\t{message}\n')
        # add to in-memory
        logline = LogLine(level, category, uuid, message)
        self.lines.append(logline)

    def load(self) -> None:
        # check file exists
        path = Path(self.filepath)
        if not path.exists():
            path.touch()

        # load messages in file
        self.lines = []
        with open(self.filepath, 'r') as fp:
            for line in fp.readlines():
                level, category, uuid, message = line.strip('\n').split('\t')
                logline = LogLine(level, category, uuid, message)
                self.lines.append(logline)

    