
from typing import Optional
from pathlib import Path
from dataclasses import dataclass
from .enums import ErrorLevel
from .enums import ErrorCategory


@dataclass
class LogLine:
    level: ErrorLevel
    message: str
    category: Optional[ErrorCategory]
    tool_uuid: Optional[str]
    subsection: Optional[str]

    def __str__(self) -> str:
        level = self.level.value
        category = self.category.value if self.category is not None else None
        return f'{level}\t{category}\t{self.tool_uuid}\t{self.subsection}\t{self.message}'


class LogFile:
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.lines: list[LogLine] = []
        self.load()

    def add(
        self, 
        level: ErrorLevel, 
        category: Optional[ErrorCategory], 
        tool_uuid: Optional[str], 
        msg: str, 
        subsection: Optional[str]
        ) -> None:

        # format message 
        message = msg.replace('\n', '')
        message = fr'{message}'

        # create logline
        logline = LogLine(
            level=level, 
            message=message, 
            category=category, 
            tool_uuid=tool_uuid, 
            subsection=subsection
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
        str_level, str_cat, str_tool_uuid, str_subsection, str_message = line.strip('\n').split('\t')
        level = ErrorLevel(str_level)
        category = None if str_cat == 'None' else ErrorCategory(str_cat)
        tool_uuid = None if str_tool_uuid == 'None' else str_tool_uuid
        subsection = None if str_subsection == 'None' else str_subsection
        return LogLine(level, str_message, category, tool_uuid, subsection)

    