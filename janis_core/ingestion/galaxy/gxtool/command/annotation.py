

from janis_core.ingestion.galaxy import expressions

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxtool.model import XMLParam
from janis_core.ingestion.galaxy.gxtool.model import XMLBoolParam
from janis_core.ingestion.galaxy.gxtool.model import XMLSelectParam

from .epath.ExecutionPath import ExecutionPath
from .epath.ExecutionPathAnnotator import GreedyExecutionPathAnnotator

from .Command import Command 
from .cmdstr.CommandString import CommandString
from .components.inputs.Flag import Flag
from .components.inputs.Option import Option
from .components import factory
from .update import update_command


# Utility functions to check aspects of a gxparam argument attribute

def argument_exists(gxparam: XMLParam, xmltool: XMLTool, command: Command) -> bool:
    if hasattr(gxparam, 'argument') and gxparam.argument is not None:  #type: ignore
        return True
    return False

def argument_format(gxparam: XMLParam, xmltool: XMLTool, command: Command) -> bool:
    banned_argument_chars = [' ', '/', '\\']
    if any([char in gxparam.argument for char in banned_argument_chars]):  #type: ignore
        return False
    return True

def argument_component_not_exists(gxparam: XMLParam, xmltool: XMLTool, command: Command) -> bool:
    # check whether the argument is already known as a prefix 
    prefix_components: list[Flag | Option] = []
    prefix_components += command.get_flags()
    prefix_components += command.get_options()
    for component in prefix_components:
        if gxparam.argument == component.prefix: # type: ignore
            return False
    return True

def argument_has_command_presence(gxparam: XMLParam, xmltool: XMLTool, command: Command) -> bool:
    # prefix appears in the cmd text
    if gxparam.argument in xmltool.raw_command: # type: ignore
        return True
    return False

def select_is_flag(gxparam: XMLSelectParam) -> bool:
    """checks if a galaxy XMLSelectParam is being used as a Flag() tool input"""
    values = gxparam.get_all_values(nonempty=True)
    if len(values) == 1:
        if values[0].startswith('-'):
            return True
    return False  

# def bool_is_option(gxparam: XMLBoolParam, xmltool: XMLTool) -> bool:
#     """
#     checks if a galaxy XMLBoolParam is being used as a Option() tool input.
#     eg cutadapt: $info_file is a XMLBoolParam, but in the command we see 
#     '--info-file=$info_file'. '--info-file $info_file' would also count as an Option()
#     """
#     if isinstance(gxparam, XMLBoolParam):
#         variable_fmt1: str = f'${gxparam.name}'
#         variable_fmt2: str = f'${{{gxparam.name}}}'
#         if variable_fmt1 in text or variable_fmt2 in text:
#             pass
#         #return True
#     text: str = xmltool.raw_command
#     prefix: str = gxparam.argument

#     pattern = r'(?<=\s|^)' + prefix + r'[:=][\'"]?\${?([\w_.])*?' + gxparam.name + r'}?[\'"]?(?=\s|$)'
#     if expressions.get_matches(text, pattern):
#         return True
#     return False


# Annotator classes

class ArgumentCommandAnnotator:
    def __init__(self, command: Command, xmltool: XMLTool):
        self.command = command
        self.xmltool = xmltool

    def annotate(self) -> None:
        # add any gxparams which hint they are components (or update the component)
        for gxparam in self.xmltool.inputs.list():
            if self.should_update_command_components(gxparam):
                gxparam = self.refine_argument(gxparam)
                self.update_command_components(gxparam)
    
    def should_update_command_components(self, gxparam: XMLParam) -> bool:
        checks = [
            argument_exists,
            argument_format,
            argument_component_not_exists,
            argument_has_command_presence,
        ]
        for check in checks:
            if not check(gxparam, self.xmltool, self.command):
                return False
        return True

    def refine_argument(self, gxparam: XMLParam) -> XMLParam:
        """
        'argument' attribute of gxparam not always written with 
        correct number of preceeding dashes. this aims to discover
        the correct amount by looking in the <command> section for the
        argument
        """
        old_argument: str = gxparam.argument # type: ignore
        matches = expressions.get_preceeding_dashes(
            search_term=old_argument,
            text=self.xmltool.raw_command
        )
        if matches:
            num_dashes = max(len(dashes) for dashes in matches) 
            gxparam.argument = '-' * num_dashes + old_argument # type: ignore
        return gxparam

    def update_command_components(self, gxparam: XMLParam) -> None:
        """
        gxparam definitely has 'argument' field
        the below are assumptions - galaxy XML can be written any way you want. 
        will need finer tuning / more disgression for use cases later?
        """
        assert(gxparam.argument) # type: ignore
        match gxparam:
            case XMLBoolParam():
                pass
                #self.handle_bool_param(gxparam)
            case XMLSelectParam():
                self.handle_select_param(gxparam)
            case _:
                self.handle_generic_param(gxparam)

    # def handle_bool_param(self, gxparam: XMLBoolParam) -> None:
    #     assert(gxparam.argument)
    #     if bool_is_option(gxparam, self.xmltool):
    #         option = factory.option(prefix=gxparam.argument, gxparam=gxparam)
    #         self.command.update(option) 
    #     else:
    #         flag = factory.flag(prefix=gxparam.argument, gxparam=gxparam)
    #         self.command.update(flag)

    def handle_select_param(self, gxparam: XMLSelectParam) -> None:
        assert(gxparam.argument)
        if select_is_flag(gxparam):
            flag = factory.flag(prefix=gxparam.argument, gxparam=gxparam)
            update_command(self.command, flag)
        else:
            values = [opt.value for opt in gxparam.options]
            option = factory.option(prefix=gxparam.argument, gxparam=gxparam, values=values)
            update_command(self.command, option)
        
    def handle_generic_param(self, gxparam: XMLParam) -> None:
        option = factory.option(prefix=gxparam.argument, gxparam=gxparam)
        update_command(self.command, option)



class CmdstrCommandAnnotator:

    def __init__(self, command: Command, xmltool: XMLTool, cmdstrs: list[CommandString]):
        self.command = command
        self.xmltool = xmltool
        self.cmdstrs = cmdstrs
        self.epath_count: int = 0

    def annotate(self) -> None:
        for cmdstr in self.cmdstrs:
            for epath in cmdstr.main.get_execution_paths():
                # logging.runtime_data(str(epath))
                self.extract_components(epath)

    def extract_components(self, epath: ExecutionPath) -> None:
        epath.id = self.epath_count
        epath = self.assign_epath_components(epath)
        for component in epath.get_components():
            update_command(self.command, component)
            self.epath_count += 1
    
    def assign_epath_components(self, epath: ExecutionPath) -> ExecutionPath:
        annotator = GreedyExecutionPathAnnotator(epath, self.xmltool, self.command)
        return annotator.annotate()
