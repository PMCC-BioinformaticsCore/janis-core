
import regex as re
from typing import Optional, Tuple
from janis_core.ingestion.galaxy import expressions

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxtool.model import XMLParam
from janis_core.ingestion.galaxy.gxtool.model import XMLBoolParam
from janis_core.ingestion.galaxy.gxtool.model import XMLSelectParam

from .epath.ExecutionPath import ExecutionPath
from .epath.ExecutionPathAnnotator import GreedyExecutionPathAnnotator

from .Command import Command 
from .cmdstr.CommandString import CommandString
from .cmdstr.DynamicCommandStatement import DynamicCommandStatement
from .cmdstr.RealisedTokenValues import RealisedTokens
from .components import CommandComponent
from .components import ComponentConfidence
from .components.inputs.InputComponent import InputComponent
from .components.inputs.Flag import Flag
from .components.inputs.Option import Option
from .components import factory

from .cmdstr import analysis
from .cmdstr.analysis import CmdstrReference
from .cmdstr.analysis import CmdstrReferenceType

from .update import update_command
from .tokenise import tokenise_line


"""
This module exists to analyse the command string of a galaxy tool, and identify the components.
This is the main orchestrator of turning a galaxy tool XML into an internal software 'Command' object. 

InlineBoolParamAnnotator:
    - identifies components from boolean params
    - only annotates if the boolean param has single appearence in XML <command> 
    - only looks at the local context of where the param is used in tool XML <command>

MultilineBoolParamAnnotator:
    - identifies components from boolean params
    - only annotates if the boolean param has single appearence in XML <command> matching pattern
    - only annotates if the use pattern seen in the XML <command> is a multiline bool pattern
    - only looks at the local context of where the param is used in tool XML <command>

SelectParamAnnotator:
    - identifies components from select params
    - only annotates if the select param has single appearence in XML <command>
    - only looks at the local context of where the param is used in tool XML <command>

ArgumentParamAnnotator:
    - identifies components from non-boolean, non-select params with an 'argument' attribute
    - only annotates if the param has single appearence in XML <command>
    - only looks at the local context of where the param is used in tool XML <command>

LocalCmdstrAnnotator:
    - identifies components locally by looking at where a param is used in the XML <command> string

GlobalCmdstrAnnotator:
    - identifies components directly from the XML <command> string
    - looks at the entire XML <command> string, and identifies components based on the order of arguments
"""

### HELPER FUNCS & PATTERNS - TODO move to .analysis.py ###




### IDENTIFYING COMPONENTS VIA BOOLEAN PARAMS ###

class SimpleInlineBoolAnnotator:
    def __init__(self, command: Command, xmltool: XMLTool):
        self.command = command
        self.xmltool = xmltool

    def annotate(self) -> None:
        components: list[InputComponent] = []
        available_params = [p for p in self.xmltool.inputs.list() if not self.command.gxparam_is_attached(p)]
        available_params = [p for p in available_params if isinstance(p, XMLBoolParam)]
        
        for param in available_params:
            # ignore params with more than 1 appearance in <command>
            if not analysis.single_inline_plaintext_appearence(self.xmltool, param):
                continue

            # set confidence based on whether param has argument and it looks right
            # meh this isn't particularly good
            confidence = 'low'
            if hasattr(param, 'argument') and param.argument is not None:
                if analysis.argument_resembles_prefix(param.argument, param.truevalue):
                    confidence = 'high'
                elif analysis.argument_resembles_prefix(param.argument, param.falsevalue):
                    confidence = 'high'

            # check if boolean param is flag (normal) 
            if self.looks_like_simple_flags(param):
                comps = self.handle_as_simple_flags(param)
                for comp in comps:
                    comp.set_confidence(confidence)
                components += comps
            
            # or option (weird)
            elif self.looks_like_simple_option(param):
                comp = self.handle_as_simple_option(param)
                comp.set_confidence(confidence)
                components.append(comp)

        for comp in components:
            update_command(self.command, comp)
    
    def looks_like_simple_flags(self, param: XMLBoolParam) -> bool:
        """
        <param argument="--no-indels" type="boolean" value="False" truevalue="--no-indels" falsevalue=""...
        <param name="force_se" argument="-S" type="boolean" truevalue="-S" falsevalue="" checked="False" label="Treat as single-end"/>
        <param name="check_names" argument="--check-names" type="boolean" checked="False" truevalue="--paired --check-names" falsevalue="--paired" label="Verify read names match"/>
        <param name="ignore_overlaps" argument="-x/--ignore-overlaps" type="boolean" truevalue="-x" falsevalue="" checked="False" label="Disable read-pair overlap detection" />
        <param name="skip_anomalous_read_pairs" argument="-A/--count-orphans" type="boolean" truevalue="-A" falsevalue="" checked="False" label="Do not discard anomalous read pairs" />
        """
        appearence = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)[0]

        # ensure that cmdstr appearence is a single word, and that word is a variable
        if not analysis.is_simple_variable(appearence.text):
            return False
        
        # sorting into blank and non-blank items
        items = [param.truevalue, param.falsevalue]
        blank_phrases = [i for i in items if analysis.is_blank(i)]
        value_phrases = [i for i in items if not analysis.is_blank(i)]

        # ensure that at least one item has a value
        if blank_phrases == 2:
            return False
        
        # ensure that all non-blank items are simple flags
        for phrase in value_phrases:
            if not analysis.is_simple_flags(phrase):
                return False
        
        return True

    def handle_as_simple_flags(self, param: XMLBoolParam) -> list[Flag]:
        flags: list[Flag] = []
        items = [param.truevalue, param.falsevalue]
        value_phrases = [i for i in items if not analysis.is_blank(i)]
        
        for value in value_phrases:
            prefixes = value.strip().split()
            for prefix in prefixes:
                flag = factory.flag(prefix=prefix, gxparam=param)
                flags.append(flag)
        
        return flags

    def looks_like_simple_option(self, param: XMLBoolParam) -> bool:
        """
        this is made up, but something like this
        <param argument="--reference" type="boolean" value="False" truevalue="mm10" falsevalue="hg38"...
        """
        appearence = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)[0]

        # param has argument & both options have values
        if not analysis.is_blank(param.truevalue) and not analysis.is_blank(param.falsevalue):
            raise NotImplementedError
            
        # param has argument & truevalue has value
        elif not analysis.is_blank(param.truevalue) and analysis.is_simple_phrases(param.truevalue):
            raise NotImplementedError

        # param has argument & falsevalue has value
        elif not analysis.is_blank(param.falsevalue) and analysis.is_simple_phrases(param.falsevalue):
            raise NotImplementedError

        # if argument, also check argument in cmdstr
        return False
    
    def handle_as_simple_option(self, param: XMLBoolParam) -> InputComponent:
        raise NotImplementedError



### IDENTIFYING COMPONENTS VIA PARAMS WHICH APPEAR IN <COMMAND> AS PREDEFINED PATTERNS ###

class SimpleMultilineBoolAnnotator:
    def __init__(self, command: Command, xmltool: XMLTool):
        self.command = command
        self.xmltool = xmltool

    def annotate(self) -> None:
        # components we will extract
        components: list[InputComponent] = []
        available_params = [p for p in self.xmltool.inputs.list() if not self.command.gxparam_is_attached(p)]
        available_params = [p for p in available_params if isinstance(p, XMLBoolParam)]
        
        # gather components for params
        for param in available_params:
            if self.is_simple_multiline_bool(param):
                component = self.handle_simple_multiline_bool(param)
                components.append(component)

        for component in components:
            component.set_confidence('high')
            update_command(self.command, component)
    
    def is_simple_multiline_bool(self, param: XMLParam) -> bool:
        appearences = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.MULTILINE_BOOL)
        for appearence in appearences:
            prefix = self.get_simple_prefix(appearence.text)
            if prefix:
                return True
        return False

    def handle_simple_multiline_bool(self, param: XMLParam) -> InputComponent:
        appearences = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.MULTILINE_BOOL)
        for appearence in appearences:
            prefix = self.get_simple_prefix(appearence.text)
            if prefix:
                component = factory.flag(prefix=prefix, gxparam=param)
                return component
        raise RuntimeError
         
    def get_simple_prefix(self, text: str) -> Optional[str]:
        textlines = text.split('\n')
        statement = textlines[1].strip()
        # ensure flag format
        if statement.startswith('-'):
            # ensure single word and correct pattern
            if re.match(r'^[A-Za-z0-9_-]+$', statement):
                return statement
        return None


### IDENTIFYING COMPONENTS VIA SELECT PARAMS ###

class SimpleSelectAnnotator:
    """
    for any param with which is the 'select' type, attempt to update components. 
    only looking for select params where each item is blank or a flag. 
    TODO add all other things Select param can be - see cutadapt.xml $output_selector
    """
    def __init__(self, command: Command, xmltool: XMLTool):
        self.command = command
        self.xmltool = xmltool

    def annotate(self) -> None:
        components: list[InputComponent] = []
        available_params = [p for p in self.xmltool.inputs.list() if not self.command.gxparam_is_attached(p)]
        available_params = [p for p in available_params if isinstance(p, XMLSelectParam)]

        for param in available_params:
            # ignore params with more than 1 appearance in <command>
            if not analysis.single_inline_plaintext_appearence(self.xmltool, param):
                continue
            
            # all select options are blank or flags
            if self.looks_like_simple_flag_selector(param):
                flags = self.handle_as_simple_flag_selector(param)
                for flag in flags:
                    flag.set_confidence('high')
                components += flags
                
            # select options are values, <command> reference preceeded by prefix
            elif self.looks_like_simple_option_selector(param):
                option = self.handle_as_simple_option_selector(param)
                option.set_confidence('high')
                components.append(option)
            
            # select options are prefixes and values
            elif self.looks_like_complex_option_selector(param):
                options = self.handle_as_complex_option_selector(param)
                for opt in options:
                    opt.set_confidence('low')
                components += options
        
        for comp in components:
            update_command(self.command, comp)

    def looks_like_simple_flag_selector(self, param: XMLSelectParam) -> bool:
        num_blank_options = 0
        num_flag_options = 0
        
        for option in param.options:
            if analysis.is_blank(option.value):
                num_blank_options += 1
            elif analysis.is_simple_flags(option.value):
                num_flag_options += 1
        
        # ensure all options are blank or flags, and max 1 blank option
        if num_blank_options <= 1 and num_blank_options + num_flag_options == len(param.options):
            return True
        return False
    
    def handle_as_simple_flag_selector(self, param: XMLSelectParam) -> list[Flag]:
        flags: list[Flag] = []
        
        phrases = [opt.value for opt in param.options if not analysis.is_blank(opt.value)]
        for phrase in phrases:
            prefixes = phrase.strip().split()
            for prefix in prefixes:
                flag = factory.flag(prefix=prefix, gxparam=param)
                flags.append(flag)

        return flags
    
    def looks_like_simple_option_selector(self, param: XMLSelectParam) -> bool:
        appearence = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)[0]
        if analysis.is_simple_option(appearence.text):
            return True
        return False
    
    def handle_as_simple_option_selector(self, param: XMLSelectParam) -> Option:
        appearence = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)[0]
        prefix, separator, value = analysis.extract_simple_option(appearence.text)
        option = factory.option(prefix=prefix, separator=separator, gxparam=param)
        option.values.add(value)
        return option
    
    def looks_like_complex_option_selector(self, param: XMLSelectParam) -> bool:
        # appearence = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)[0]
        raise NotImplementedError
    
    def handle_as_complex_option_selector(self, param: XMLSelectParam) -> list[Option]:
        raise NotImplementedError
        


### IDENTIFYING COMPONENTS VIA PARAMS WITH ARGUMENT ATTRIBUTE ###

class OptionParamAnnotator:
    """
    for any param which looks like option, attempt to update components. 
    must identify its use in the command, & check whether it looks like an simple option.
    """
    def __init__(self, command: Command, xmltool: XMLTool):
        self.command = command
        self.xmltool = xmltool

    def annotate(self) -> None:
        components: list[InputComponent] = []
        available_params = [p for p in self.xmltool.inputs.list() if not self.command.gxparam_is_attached(p)]
        available_params = [p for p in available_params if not isinstance(p, XMLBoolParam)]
        available_params = [p for p in available_params if not isinstance(p, XMLSelectParam)]
        
        for param in available_params:
            if self.looks_like_simple_option(param):
                component, confidence = self.handle_as_simple_option(param)
                component.set_confidence(confidence)
                components.append(component)
        
        for comp in components:
            update_command(self.command, comp)

    def looks_like_simple_option(self, param: XMLParam) -> bool:
        appearences = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)
        
        # ensure single inline plain text appearence (all appearences are the same)
        # need to do it this way due to macros or reuse of the same param in different logic blocks
        lines_set = set([x.text.strip() for x in appearences])
        if not len(lines_set) == 1:
            return False
        
        # ensure its "--prefix $param.name" or "--prefix=$param.name" etc
        appearence = appearences[0]
        if not analysis.is_simple_option(appearence.text):
            return False

        return True

    def handle_as_simple_option(self, param: XMLParam) -> Tuple[InputComponent, str]:
        # create option
        appearence = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)[0]
        prefix, separator, value = analysis.extract_simple_option(appearence.text)
        option = factory.option(prefix=prefix, separator=separator, gxparam=param)
        option.values.add(value)
        
        # confidence based on argument resembling prefix
        confidence = 'low'
        if hasattr(param, 'argument') and param.argument is not None:
            if analysis.argument_resembles_prefix(param.argument, prefix):  # type: ignore
                confidence = 'high'

        return option, confidence

    

### IDENTIFYING COMPONENTS VIA LOCAL PARAM CONTEXT IN COMMAND ###

class LocalCmdstrAnnotator:

    def __init__(self, command: Command, xmltool: XMLTool):
        self.command = command
        self.xmltool = xmltool

    def annotate(self) -> None:
        identified_components: list[CommandComponent] = []
        available_params = [p for p in self.xmltool.inputs.list() if not self.command.gxparam_is_attached(p)]

        # get INLINE_PLAIN_TEXT uses of params
        # for each, prepend the previous word (not including control structure lines)(token)
        # this forms the local text to parse
        # turn local text into list of epaths
        # for each epath, extract components & pool those annotated with the param
        # for each component in the pool, if they all agree (are the same component type,
        # have the same attributes eg prefix & type etc), select first & add to command

        for param in available_params:
            pooled_components: list[CommandComponent] = []
            appearences = analysis.get_cmdstr_appearences(self.xmltool.raw_command, param, filter_to=CmdstrReferenceType.INLINE_PLAIN_TEXT)
            
            for appearence in appearences:
                epath_count = 0
                epaths = self.get_epaths(appearence)
                for epath in epaths:
                    components = self.extract_components(epath, epath_count)
                    pooled_components += self.get_param_linked_components(components, param)
                    epath_count += 1

            if self.components_agree(pooled_components):
                pooled_components[0].set_confidence('high')
                identified_components.append(pooled_components[0])

        for comp in identified_components:
            update_command(self.command, comp)

    def get_epaths(self, appearence: CmdstrReference) -> list[ExecutionPath]:
        # turn local text into annotated epaths
        realised_tokens = tokenise_line(appearence.text, self.xmltool)
        dynamicstmt = DynamicCommandStatement(appearence.text, realised_tokens)
        epaths = list(dynamicstmt.get_execution_paths())
        return epaths
    
    def extract_components(self, epath: ExecutionPath, epath_count: int) -> list[CommandComponent]:
        epath.id = epath_count
        annotator = GreedyExecutionPathAnnotator(epath, self.xmltool, self.command)
        annotator.annotate()
        return epath.get_components()

    def get_param_linked_components(self, components: list[CommandComponent], param: XMLParam) -> list[CommandComponent]:
        # for each epath, extract components annotated with the provided param
        out: list[CommandComponent] = []
        for comp in components:
            if comp.gxparam and comp.gxparam.name == param.name:
                out.append(comp)
        return out
    
    def components_agree(self, components: list[CommandComponent]) -> bool:
        # for each epath, extract components annotated with the provided param
        if not components:
            return False
        
        ctypes_set = set([c.__class__.__name__ for c in components])
        dtypes_set = set([c.datatype.classname for c in components])

        # component types must all be the same
        if not len(ctypes_set) == 1:
            return False
        # data types must all be the same
        if not len(dtypes_set) == 1:
            return False
        # if option, all prefixes must be the same
        if all([isinstance(c, Option) for c in components]):
            prefixes_set = set([c.prefix for c in components])
            if not len(prefixes_set) == 1:
                return False
        
        return True
            



### IDENTIFYING COMPONENTS VIA GLOBAL COMMAND LINE ARGUMENT ORDER ###

class GlobalCmdstrAnnotator:

    def __init__(self, command: Command, xmltool: XMLTool, cmdstrs: list[CommandString]):
        self.command = command
        self.xmltool = xmltool
        self.cmdstrs = cmdstrs
        self.epath_count: int = 0

    def annotate(self) -> None:
        for cmdstr in self.cmdstrs:
            for epath in cmdstr.main.get_execution_paths():
                self.extract_components(epath)

    def extract_components(self, epath: ExecutionPath) -> None:
        final_pass = True if self.epath_count == len(self.cmdstrs) - 1 else False
        epath.id = self.epath_count
        epath = self.assign_epath_components(epath)
        for component in epath.get_components():
            update_command(self.command, component, final_pass)
        self.epath_count += 1
    
    def assign_epath_components(self, epath: ExecutionPath) -> ExecutionPath:
        annotator = GreedyExecutionPathAnnotator(epath, self.xmltool, self.command)
        return annotator.annotate()
