


from dataclasses import dataclass


@dataclass
class Logline:
    """
    represents a single log line within a .log file
    structure = '[%(levelname)s] [%(name)s] %(message)s'
    """
    level: str
    logger: str
    message: str
    contents: str


@dataclass
class Logfile:
    """represents a single .log file"""
    subtype: str    
    path: str
    lines: list[Logline]

    @property
    def messages(self) -> list[str]:
        messages = [line.message for line in self.lines]
        return list(set(messages))

