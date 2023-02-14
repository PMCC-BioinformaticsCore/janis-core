


from ..gxtool.tool import XMLToolDefinition
from ..gxtool.text import load_command

from .cmdstr.CommandString import CommandString
from .cmdstr.generate import gen_command_string
from .annotation import ArgumentCommandAnnotator
from .annotation import CmdstrCommandAnnotator
from .Command import Command


class CommandFactory:
    def __init__(self, xmltool: XMLToolDefinition):
        self.xmltool = xmltool
        self.xmlcmdstr = self.gen_cmdstr_from_xml()
        self.command = Command(self.xmlcmdstr)

    def create(self) -> Command:
        self.update_command_via_arguments()
        self.update_command_via_cmdstrs()
        return self.command

    def update_command_via_arguments(self) -> None:
        """uses galaxy params with an 'argument' attribute to update command"""
        annotator = ArgumentCommandAnnotator(self.command, self.xmltool)
        annotator.annotate()
    
    def update_command_via_cmdstrs(self) -> None:
        """
        uses valid command line strings from tests, and the tool XML <command> section
        to further identify the structure and options of the underling software tool
        """
        # create command strings (from evaluated tests, simplified xml <command>)
        cmdstrs = self.gen_cmdstrs()
        annotator = CmdstrCommandAnnotator(self.command, self.xmltool, cmdstrs)
        annotator.annotate()

    def gen_cmdstrs(self) -> list[CommandString]:
        # note ordering: xml then test
        return [self.xmlcmdstr]

    def gen_cmdstr_from_xml(self) -> CommandString:
        text = load_command()
        return gen_command_string(source='xml', the_string=text, xmltool=self.xmltool)


