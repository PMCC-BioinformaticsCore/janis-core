
from copy import deepcopy

from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.expressions.patterns import LINUX_STATEMENT_DELIMS

from .CommandString import CommandString
from .CommandString import CommandStringSource
from .DynamicCommandStatement import DynamicCommandStatement
from .MainStatementInferrer import MainStatementInferrer
from .RealisedTokenValues import RealisedTokenFactory


def gen_command_string(source: CommandStringSource, text: str, xmltool: XMLTool) -> CommandString:
    return CommandStringGenerator(source, text, xmltool).generate()


class CommandStringGenerator:

    def __init__(self, source: CommandStringSource, text: str, xmltool: XMLTool) -> None:
        self.source = source
        self.text = text
        self.xmltool = xmltool

    def generate(self) -> CommandString:
        # get all statements in <command>
        statements = self.gen_command_statements()

        # split statements into pre-main-post
        if self.xmltool.metadata.main_requirement is not None:
            requirement = self.xmltool.metadata.main_requirement.name
        else:
            requirement = self.xmltool.metadata.id
        statement_dict = self.split_pre_main_post_statements(statements, requirement)

        # create command string & return
        return self.init_command_string(statement_dict)

    def gen_command_statements(self) -> list[DynamicCommandStatement]:
        statements: list[DynamicCommandStatement] = []
        for cmdstmt in self.split_text_statements():
            factory = RealisedTokenFactory(self.xmltool)
            realised_tokens = factory.try_tokenify(cmdstmt)
            dynamicstmt = DynamicCommandStatement(cmdstmt, realised_tokens)
            statements.append(dynamicstmt)
        return statements

    def split_text_statements(self) -> list[str]:
        text = deepcopy(self.text)
        statements: list[str] = []
        
        delim_matches = expressions.get_matches(text, LINUX_STATEMENT_DELIMS)
        quoted_sections = expressions.get_quoted_sections(text)

        # has to be reverse order otherwise m.start() and m.end() are out of place
        for m in sorted(delim_matches, key=lambda x: x.start(), reverse=True): 
            if quoted_sections[m.start()] == False and quoted_sections[m.end()] == False:
                left_split = text[:m.start()]
                right_split = text[m.end():]
                statements = [right_split] + statements # prepend
                text = left_split

        statements = [text] + statements # prepend the final split (left this time)
        return statements

    def split_pre_main_post_statements(self, statements: list[DynamicCommandStatement], requirement: str) -> dict[str, list[DynamicCommandStatement]]:
        out: dict[str, list[DynamicCommandStatement]] = {'pre': [], 'main': [], 'post': []}
        inferrer = MainStatementInferrer(statements, self.source, requirement)
        main_index = inferrer.infer()
        
        for i, statement in enumerate(statements):
            if i < main_index:
                out['pre'].append(statement)
            elif i == main_index:
                out['main'].append(statement)
            else:
                out['post'].append(statement)

        return out

    def init_command_string(self, statement_dict: dict[str, list[DynamicCommandStatement]]) -> CommandString:
        if statement_dict['pre']:
            pass
            # logging.has_preprocessing()
        if statement_dict['post']:
            pass
            # logging.has_postprocessing()
        return CommandString(
            self.source,
            main=statement_dict['main'][0], 
            preprocessing=statement_dict['pre'], 
            postprocessing=statement_dict['post'],
        )

