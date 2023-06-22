

from typing import Optional, Any, Tuple
from collections import defaultdict

from janis_core import settings

from janis_core.ingestion.galaxy.gxtool.command import Command
from janis_core.ingestion.galaxy.gxtool.command.components import InputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import OutputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import RedirectOutput
from janis_core.ingestion.galaxy.gxtool.command.components import Flag
from janis_core.ingestion.galaxy.gxtool.command.components import factory as component_factory
from janis_core.ingestion.galaxy.gxworkflow import load_tool_state
from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxtool.model import (
    XMLParam, 
    XMLInputParam,
    XMLBoolParam,
    XMLOutputParam
)

# this module imports
from .tool import ITool


def gen_tool(
    xmltool: XMLTool, 
    command: Command, 
    container: Optional[str], 
    gxstep: Optional[dict[str, Any]]=None
    ) -> ITool:
    tfactory = ToolFactory(xmltool, command, container, gxstep)
    tool = tfactory.create()
    return tool


class ToolFactory:
    def __init__(
        self, 
        xmltool: XMLTool, 
        command: Command, 
        container: Optional[str],
        gxstep: Optional[dict[str, Any]]=None
        ) -> None:
    
        self.xmltool = xmltool
        self.command = command
        self.container = container
        self.gxstep = gxstep

    def create(self) -> ITool:
        tool = ITool(
            metadata=self.xmltool.metadata,
            container=self.container,
            base_command=self.get_base_command(),
            # gxparam_register=self.xmltool.inputs,
            configfiles=self.xmltool.configfiles,
            scripts=self.xmltool.scripts,
        )
        self.supply_inputs(tool)
        self.supply_outputs(tool)
        return tool

    def supply_inputs(self, tool: ITool) -> None:
        extractor = InputExtractor(self.xmltool, self.command, self.gxstep)
        inputs = extractor.extract()
        for inp in inputs:
            tool.add_input(inp)

    def supply_outputs(self, tool: ITool) -> None:
        extractor = OutputExtractor(self.xmltool, self.command, self.gxstep)
        outputs = extractor.extract()
        for out in outputs:
            tool.add_output(out)

    def get_base_command(self) -> list[str]:
        positionals = self.command.get_base_positionals()
        if not positionals:
            pass
            # logging.no_base_cmd()
        return [p.default_value for p in positionals]
    

class InputExtractor:
    def __init__(
        self, 
        xmltool: XMLTool, 
        command: Command, 
        gxstep: Optional[dict[str, Any]]=None
        ) -> None:
        self.xmltool = xmltool
        self.command = command
        self.gxstep = gxstep
        self.inputs: list[InputComponent] = []

    def extract(self) -> list[InputComponent]:
        # getting the inputs is a bit of a mess, so we'll do it in stages
        # inputs we identified in the tool command
        self.add_captured_inputs()
        # inputs we didn't identify in the tool command
        self.add_uncaptured_inputs()
        return self.inputs

    def add_captured_inputs(self) -> None:
        # get the inputs we identified in the tool command
        self.command.set_cmd_positions()
        self.inputs += self.command.list_inputs(include_base_cmd=False)
    
    def add_uncaptured_inputs(self) -> None:
        # get the tool xml inputs we didn't capture by analyzing the tool command
        # ie we know they're tool inputs, but don't know how they wire to args in the tool command
        if not self.gxstep:
            return None
        
        tool_state = load_tool_state(self.gxstep, additional_filters=['Flatten', 'DeNestClass'])
        for param_name, param_value in tool_state.items():
            if self.should_create_uncaptured_input(param_name, param_value): # type: ignore
                param = self.xmltool.inputs.get(param_name, strategy="lca")
                inp = self.create_uncaptured_input(param) # type: ignore
                self.inputs.append(inp)

    def should_create_uncaptured_input(self, param_name: str, param_value: str) -> bool:
        param = self.xmltool.inputs.get(param_name, strategy="lca")
        if not param:
            return False
        
        for inp in self.inputs:
            if inp.gxparam and inp.gxparam.name == param.name:
                return False
            
        return True
        
        # if param_value in ['ConnectedValue', 'RuntimeValue']:
        #     return True
        
        # return False
    
    def create_uncaptured_input(self, gxparam: XMLInputParam) -> InputComponent:
        if isinstance(gxparam, XMLBoolParam):
            return self.create_component_for_boolparam(gxparam)
        elif gxparam.argument:
            return component_factory.option(prefix=gxparam.argument, gxparam=gxparam)
        else:
            return component_factory.positional(gxparam=gxparam)
        
    def create_component_for_boolparam(self, gxparam: XMLBoolParam) -> InputComponent:
        truevalue_prefix = gxparam.truevalue.strip().split(' ')[0]
        falsevalue_prefix = gxparam.falsevalue.strip().split(' ')[0]

        # prefix held in truevalue
        if truevalue_prefix.startswith('-') and falsevalue_prefix == '':
            return component_factory.flag(prefix=truevalue_prefix, gxparam=gxparam)
        
        # prefix held in falsevalue
        elif falsevalue_prefix.startswith('-') and truevalue_prefix == '':
            return component_factory.flag(prefix=falsevalue_prefix, gxparam=gxparam)

        # same prefix appears in both truevalue / falsevalue
        elif (
            truevalue_prefix.startswith('-') 
            and falsevalue_prefix.startswith('-')
            and truevalue_prefix == falsevalue_prefix
        ):
            return component_factory.flag(prefix=truevalue_prefix, gxparam=gxparam)
        
        # different prefix appears in both truevalue / falsevalue
        # treat as text instead of bool
        elif (
            truevalue_prefix.startswith('-') 
            and falsevalue_prefix.startswith('-')
            and truevalue_prefix != falsevalue_prefix
        ):
            return component_factory.positional(gxparam=gxparam)
        
        # doesn't seem to be a prefix in either.
        # treat as text instead of bool
        else:
            return component_factory.positional(gxparam=gxparam)



class OutputExtractor:
    def __init__(
        self, 
        xmltool: XMLTool, 
        command: Command, 
        gxstep: Optional[dict[str, Any]]=None
        ) -> None:
        self.xmltool = xmltool
        self.command = command
        self.gxstep = gxstep
        self.whitelisted_outputs: list[XMLParam] = []
        self.tool_outputs: dict[str, list[OutputComponent]] = defaultdict(list)

    def extract(self) -> list[OutputComponent]:
        outputs: list[OutputComponent] = []
        self.define_whitelisted_outputs()
        self.gather_input_outputs()
        self.gather_wildcard_outputs()
        self.gather_redirect_outputs()
        self.gather_uncertain_outputs()
        outputs = self.prioritise_outputs()
        return outputs

    def prioritise_outputs(self) -> list[OutputComponent]:
        # sorry this func is last minute and horrendous
        prioritised: list[OutputComponent] = []

        # set up data_structure so we can look up the possible outputs by name
        data_structure: dict[str, list[Tuple[str, OutputComponent]]] = defaultdict(list)
        for otype, outputs in self.tool_outputs.items():
            for out in outputs:
                assert(out.gxparam)
                data_structure[out.gxparam.name].append((otype, out))
        
        priorities = {
            'redirect': 0,
            'input': 1,
            'wildcard': 2,
            'uncertain': 3,
        }

        for out in self.whitelisted_outputs:
            if out.name not in data_structure:
                print()
            possible = data_structure[out.name]
            possible_sorted = sorted(possible, key=lambda x: priorities[x[0]])
            prioritised.append(possible_sorted[0][1])
            
        return prioritised
        
    def define_whitelisted_outputs(self) -> None:
        if settings.translate.MODE in ['skeleton', 'regular'] and self.gxstep:
            for out in self.gxstep['outputs']:
                param = self.xmltool.outputs.get(out['name']) 
                if param:
                    self.whitelisted_outputs.append(param)
        else:
            self.whitelisted_outputs = self.xmltool.outputs.list()
        
    def gather_wildcard_outputs(self) -> None:
        # verified vs unverified:
        # only if no post-processing! otherwise, wildcard outputs may 
        # have come from post-processing.
        
        # outputs which were not identified in the command
        # usually just because they have a file collection strategy
        # like from_work_dir or a <discover_datatsets> tag as a child
        for gxparam in self.whitelisted_outputs:
            if self.should_create_wildcard_output(gxparam):
                output = component_factory.wildcard_output(gxparam)
                self.tool_outputs['wildcard'].append(output)
    
    def gather_redirect_outputs(self) -> None:
        # redirect outputs (stdout) already identified when creating Command()
        # need to ensure they're linked to a gxparam
        # if not, try to link to dataset collector, else just ignore as the 
        # redirect seems to be dropped.
        redirects: list[RedirectOutput] = self.command.list_outputs()
        for r in redirects:
            self.attempt_redirect_gxparam_link(r)
            if r.gxparam is not None:
                self.tool_outputs['redirect'].append(r)

    def attempt_redirect_gxparam_link(self, r: RedirectOutput) -> None:
        if not r.gxparam: 
            for query_param in self.whitelisted_outputs:
                if query_param.discover_pattern is not None:
                    if query_param.discover_pattern == r.values.most_common_value:
                        r.gxparam = query_param
    
    def gather_input_outputs(self) -> None:
        # can be identified by looking at the input components which
        # have attached gxparams which are outputs. 
        # example case: toolname input.fastq -o $outfile
        # the -o option would be picked up as a command component, and the
        # gxparam referred to by $outfile is stored on the command component
        for component in self.command.list_inputs(include_base_cmd=False):
            if self.should_create_input_output(component):
                output = component_factory.input_output(component)
                self.tool_outputs['input'].append(output)
    
    def gather_uncertain_outputs(self) -> None:
        for gxparam in self.whitelisted_outputs:
            output = component_factory.uncertain_output(gxparam)
            self.tool_outputs['uncertain'].append(output)

    # CHECKS
    def should_create_input_output(self, component: OutputComponent) -> bool:
        if not isinstance(component, Flag):
            if isinstance(component.gxparam, XMLInputParam):
                return True
            elif isinstance(component.gxparam, XMLOutputParam):
                if component.gxparam.name in [x.name for x in self.whitelisted_outputs]:
                    return True
        return False

    def should_create_wildcard_output(self, gxparam: XMLParam) -> bool:
        """test to see if this *galaxy output param* should spawn WildcardOutput"""
        if hasattr(gxparam, 'from_work_dir') and gxparam.from_work_dir is not None: # type: ignore
            return True
        elif hasattr(gxparam, 'discover_pattern') and gxparam.discover_pattern is not None: # type: ignore
            return True
        return False

    

