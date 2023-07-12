
from typing import Any
from dataclasses import dataclass

@dataclass 
class XMLTest:
    name: str
    inputs: dict[str, Any]

class XMLTestRegister:
    def __init__(self, gxtests: list[XMLTest]):
        self.tests: dict[str, XMLTest] = {t.name: t for t in gxtests}

    def list(self) -> list[XMLTest]:
        return list(self.tests.values())