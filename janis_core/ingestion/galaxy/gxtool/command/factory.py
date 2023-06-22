

from typing import Optional, Any

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxworkflow.tool_state.load import load_tool_state

from .load import load_vanilla_command_str
from .load import load_templated_command_str
from .cmdstr.CommandString import CommandString
from .cmdstr.CommandString import CommandStringSource
from .cmdstr.generate import gen_command_string
from .annotation import ArgumentCommandAnnotator
from .annotation import CmdstrCommandAnnotator
from .Command import Command


"""
Generates a Command object from a Galaxy XML tool definition.
The Command object stores the positional, optional, and flag arguments of the software tool. 
"""

class CommandFactory:
    def __init__(self, xmltool: XMLTool, gxstep: Optional[dict[str, Any]]=None):
        self.xmltool = xmltool
        self.gxstep = gxstep
        self.command = Command()

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
        uses valid command line strings from tests, tool state, and the tool XML <command> section
        to further identify the structure and options of the underling software tool
        """
        # create command strings (from evaluated tests, tool state, simplified xml <command>)
        cmdstrs = self.gen_cmdstrs()
        # update self.command with info from cmdstrs
        annotator = CmdstrCommandAnnotator(self.command, self.xmltool, cmdstrs)
        annotator.annotate()
    
    def gen_cmdstrs(self) -> list[CommandString]:
        cmdstrs: list[CommandString] = []
        
        # vanilla xml
        text = load_vanilla_command_str()
        print(text)
        cmdstr = gen_command_string(source=CommandStringSource.XML, text=text, xmltool=self.xmltool)
        cmdstrs.append(cmdstr)
        
        # templated tests
        for test in self.xmltool.tests.list():
            text = load_templated_command_str(test.inputs)
            print(text)
            cmdstr = gen_command_string(source=CommandStringSource.TEST, text=text, xmltool=self.xmltool)
            cmdstrs.append(cmdstr)
        
        # templated tool state  NOTE unsure on ordering - tests first, or tool state first?
        if self.gxstep:
            inputs_dict = load_tool_state(
                self.gxstep, 
                additional_filters=[
                    'ReplaceNullWithVarname'
                    'ReplaceConnectedWithVarname',
                    'ReplaceRuntimeWithVarname',
                ]
            )
            text = load_templated_command_str(inputs_dict)
            print(text)
            cmdstr = gen_command_string(source=CommandStringSource.TOOL_STATE, text=text, xmltool=self.xmltool)
            cmdstrs.append(cmdstr)

        return cmdstrs






