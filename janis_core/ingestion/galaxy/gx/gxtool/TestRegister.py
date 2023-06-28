




from janis_core.tool.test_classes import TTestCase


class TestRegister:
    def __init__(self, gxtests: list[TTestCase]):
        self.tests: dict[str, TTestCase] = {t.name: t for t in gxtests}

    def list(self) -> list[TTestCase]:
        return list(self.tests.values())



