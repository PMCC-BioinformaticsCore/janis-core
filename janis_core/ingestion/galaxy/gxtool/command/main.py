

from typing import Optional, Any

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxworkflow.tool_state.load import load_tool_state

from .loading import load_vanilla_command_str
from .loading import load_templated_command_str
from .cmdstr.CommandString import CommandString
from .cmdstr.CommandString import CommandStringSource
from .cmdstr.generate import gen_command_string
from .annotation import SimpleInlineBoolAnnotator
from .annotation import SimpleMultilineBoolAnnotator
from .annotation import SimpleSelectAnnotator
from .annotation import OptionParamAnnotator
from .annotation import LocalCmdstrAnnotator
from .annotation import GlobalCmdstrAnnotator
from .Command import Command

"""
Generates a Command object from a Galaxy XML tool definition.
The Command object stores the positional, optional, and flag arguments of the software tool. 
"""

def gen_command(xmltool: XMLTool, gxstep: Optional[dict[str, Any]]=None) -> Command:
    factory = CommandFactory(xmltool, gxstep)
    return factory.create()


class CommandFactory:
    def __init__(self, xmltool: XMLTool, gxstep: Optional[dict[str, Any]]=None):
        self.xmltool = xmltool
        self.gxstep = gxstep
        self.command = Command()

    def create(self) -> Command:
        SimpleInlineBoolAnnotator(self.command, self.xmltool).annotate()
        SimpleMultilineBoolAnnotator(self.command, self.xmltool).annotate()
        SimpleSelectAnnotator(self.command, self.xmltool).annotate()
        OptionParamAnnotator(self.command, self.xmltool).annotate()
        LocalCmdstrAnnotator(self.command, self.xmltool).annotate()
        GlobalCmdstrAnnotator(self.command, self.xmltool, self.gen_cmdstrs()).annotate()
        return self.command
    
    def gen_cmdstrs(self) -> list[CommandString]:
        cmdstrs: list[CommandString] = []
        
        # vanilla xml
        text = load_vanilla_command_str()
        cmdstr = gen_command_string(source=CommandStringSource.XML, text=text, xmltool=self.xmltool)
        cmdstrs.append(cmdstr)
        
        # templated tests
        for test in self.xmltool.tests.list():
            text = load_templated_command_str(test.inputs)
            cmdstr = gen_command_string(source=CommandStringSource.TEST, text=text, xmltool=self.xmltool)
            cmdstrs.append(cmdstr)
        
        # templated tool state  NOTE unsure on ordering - tests first, or tool state first?
        if self.gxstep:
            # TODO HERE
            inputs_dict = load_tool_state(
                self.gxstep, 
                additional_filters=[
                    'ReplaceNullWithVarname'
                    'ReplaceConnectedWithVarname',
                    'ReplaceRuntimeWithVarname',
                ]
            )
            text = load_templated_command_str(inputs_dict)
            cmdstr = gen_command_string(source=CommandStringSource.TOOL_STATE, text=text, xmltool=self.xmltool)
            cmdstrs.append(cmdstr)

        return cmdstrs






