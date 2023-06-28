

from enum import Enum, auto

class Stream(Enum):
    STDIN = auto()
    STDOUT = auto()
    STDERR = auto()
    BOTH = auto()


class StreamMerge:
    def __init__(self, text: str):
        self.source: Stream = self.extract_source(text)
        self.destination: Stream = self.extract_dest(text)

    def extract_source(self, text: str) -> Stream:
        if text[0] == '2':
            return Stream.STDERR
        return Stream.STDOUT

    def extract_dest(self, text: str) -> Stream:
        if text[-1] == '2':
            return Stream.STDERR
        return Stream.STDOUT
