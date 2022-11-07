


from janis_core.types import File


class SecondaryTestType(File):
    def __init__(self, optional=False):
        super().__init__(optional, extension=".txt")

    @staticmethod
    def secondary_files():
        return ["^.txt"]


class AppendedSecondaryTestType(File):
    def __init__(self, optional=False):
        super().__init__(optional, extension=".bam")

    @staticmethod
    def secondary_files():
        return [".bai"]


class ReplacedSecondaryTestType(File):
    def __init__(self, optional=False):
        super().__init__(optional, extension=".bam")

    @staticmethod
    def secondary_files():
        return ["^.bai"]


class NonEscapedSecondaryTestType(File):
    @staticmethod
    def secondary_files():
        return [".txt"]