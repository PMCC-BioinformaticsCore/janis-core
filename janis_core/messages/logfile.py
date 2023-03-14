

from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class LogLine:
    level: str
    uuid: str
    message: str


class LogFile:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.messages: dict[str, list[LogLine]] = defaultdict(list)
        self.load()

    def add(self, level: str, uuid: str, message: str) -> None:
        # add to file
        message = fr'{message}'
        with open(self.filepath, 'a') as fp:
            fp.write(f'{level}\t{uuid}\t{message}\n')
        # add to in-memory
        logline = LogLine(level, uuid, message)
        self.messages[uuid].append(logline)

    def load(self) -> None:
        # check file exists
        path = Path(self.filepath)
        if not path.exists():
            path.touch()

        # load messages in file
        with open(self.filepath, 'r') as fp:
            for line in fp.readlines():
                level, uuid, message = line.strip('\n').split('\t')
                logline = LogLine(level, uuid, message)
                self.messages[uuid].append(logline)

    