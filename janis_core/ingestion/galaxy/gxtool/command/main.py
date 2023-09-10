

from typing import Optional, Any

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxworkflow.tool_state.load import load_tool_state

from .loading import load_vanilla_command_str
from .loading import load_templated_command_str
from .cmdstr.DynamicCommandStatement import DynamicCommandStatement
# from .cmdstr.CommandString import CommandString
# from .cmdstr.CommandString import CommandStringSource
# from .cmdstr.generate import gen_command_string
from .annotation import SimpleInlineBoolAnnotator
from .annotation import SimpleMultilineBoolAnnotator
from .annotation import SimpleSelectAnnotator
from .annotation import OptionParamAnnotator
from .annotation import LocalCmdstrAnnotator
from .annotation import GlobalCmdstrAnnotator
from .tokenise import tokenise_text
from .Command import Command

"""
Generates a Command object from a Galaxy XML tool definition.
The Command object stores the positional, optional, and flag arguments of the software tool. 
"""

def gen_command(
    xmltool: XMLTool, 
    gxstep: Optional[dict[str, Any]]=None,
    annotators: Optional[list[str]]=None
    ) -> Command:
    factory = CommandFactory(xmltool, gxstep, annotators)
    return factory.create()

class CommandFactory:
    DEFAULT_ANNOTATORS = [
        'SimpleInlineBoolAnnotator',
        'SimpleMultilineBoolAnnotator',
        'SimpleSelectAnnotator',
        'OptionParamAnnotator',
        'LocalCmdstrAnnotator',
        'GlobalCmdstrAnnotator',
    ]

    def __init__(
        self, 
        xmltool: XMLTool, 
        gxstep: Optional[dict[str, Any]]=None,
        annotators: Optional[list[str]]=None
        ) -> None:
        self.xmltool = xmltool
        self.gxstep = gxstep
        self.annotators = annotators if annotators else self.DEFAULT_ANNOTATORS
        self.command = Command()

    def create(self) -> Command:
        # split into pre, main, post
        # remove post 
        # all annotators except Local/GlobalCmdstrAnnotator - supply main as text
        # Local/GlobalCmdstrAnnotator - supply supply pre & main, only start greedy search from main
        cmdtext = load_vanilla_command_str(self.xmltool)
        mainstmt_text = cmdtext.split('__JANIS_MAIN__')[1]

        if 'SimpleInlineBoolAnnotator' in self.annotators:
            SimpleInlineBoolAnnotator(self.command, mainstmt_text, self.xmltool).annotate()
        if 'SimpleMultilineBoolAnnotator' in self.annotators:
            SimpleMultilineBoolAnnotator(self.command, mainstmt_text, self.xmltool).annotate()
        if 'SimpleSelectAnnotator' in self.annotators:
            SimpleSelectAnnotator(self.command, mainstmt_text, self.xmltool).annotate()
        if 'OptionParamAnnotator' in self.annotators:
            OptionParamAnnotator(self.command, mainstmt_text, self.xmltool).annotate()
        if 'LocalCmdstrAnnotator' in self.annotators:
            LocalCmdstrAnnotator(self.command, mainstmt_text, self.xmltool).annotate()
        if 'GlobalCmdstrAnnotator' in self.annotators:
            stmts_dynamic = self.get_dynamic_statements()
            GlobalCmdstrAnnotator(self.command, self.xmltool, stmts_dynamic).annotate()
        return self.command
    
    def get_dynamic_statements(self) -> list[DynamicCommandStatement]:
        stmts_dynamic = []

        # cheetah templating if galaxy step tool state present
        if self.gxstep:
            inputs_dict = load_tool_state(
                self.xmltool,
                self.gxstep, 
                additional_filters=[
                    'ReplaceNullWithVarname',
                    'ReplaceBoolWithValue',
                    'ReplaceConnectedWithVarname',
                    'ReplaceRuntimeWithVarname',
                ]
            )
            cmdtext = load_templated_command_str(self.xmltool, inputs_dict)
            mainstmt_text = cmdtext.split('__JANIS_MAIN__')[1]    
            mainstmt_tokens = tokenise_text(mainstmt_text, self.xmltool)
            mainstmt_dynamic = DynamicCommandStatement(mainstmt_text, mainstmt_tokens)
            stmts_dynamic.append(mainstmt_dynamic)

        # vanilla xml
        cmdtext = load_vanilla_command_str(self.xmltool)
        mainstmt_text = cmdtext.split('__JANIS_MAIN__')[1]    
        mainstmt_tokens = tokenise_text(mainstmt_text, self.xmltool)
        mainstmt_dynamic = DynamicCommandStatement(mainstmt_text, mainstmt_tokens)
        stmts_dynamic.append(mainstmt_dynamic)
        
        return stmts_dynamic







