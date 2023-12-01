
from typing import Optional
from pathlib import Path
from dataclasses import dataclass
from .enums import ErrorCategory


@dataclass
class LogLine:
    message: str
    category: Optional[ErrorCategory]
    entity_uuid: Optional[str]

    def __str__(self) -> str:
        if self.category is not None:
            level, cat = self.category.value
        else:
            level, cat = None, None
        return f'{level}\t{cat}\t{self.entity_uuid}\t{self.message}'


class LogFile:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.lines: list[LogLine] = []
        self.load()

    def add(
        self, 
        category: Optional[ErrorCategory], 
        entity_uuid: Optional[str], 
        msg: str, 
        ) -> None:

        # format message 
        message = msg.replace('\n', '')
        message = fr'{message}'

        # create logline
        logline = LogLine(
            message=message, 
            category=category, 
            entity_uuid=entity_uuid, 
        )
        
        # write to file
        with open(self.filepath, 'a') as fp:
            fp.write(f'{str(logline)}\n')

        # add to in-memory?
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
                logline = self.string_to_logline(line)
                self.lines.append(logline)

    def string_to_logline(self, line: str) -> LogLine:
        str_level, str_cat, str_uuid, str_message = line.strip('\n').split('\t')
        if str_cat == 'None':
            category = None
        else:
            category = ErrorCategory.from_str(str_cat)
        entity_uuid = None if str_uuid == 'None' else str_uuid
        return LogLine(str_message, category, entity_uuid)

    