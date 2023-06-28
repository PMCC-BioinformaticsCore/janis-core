


from ..gxtool.tool import XMLToolDefinition
from .factory import CommandFactory
from .Command import Command


def gen_command(xmltool: XMLToolDefinition) -> Command:
    factory = CommandFactory(xmltool)
    return factory.create()



