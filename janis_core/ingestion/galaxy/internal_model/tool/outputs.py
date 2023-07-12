

from typing import Optional, Any, Tuple
from collections import defaultdict
from janis_core.ingestion.galaxy.gxtool.command import Command
from janis_core.ingestion.galaxy.gxtool.command.components import CommandComponent
from janis_core.ingestion.galaxy.gxtool.command.components import RedirectOutput
from janis_core.ingestion.galaxy.gxtool.command.components import Flag
from janis_core.ingestion.galaxy.gxtool.command.components import factory
from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxtool.model import XMLParam
from janis_core.ingestion.galaxy.gxtool.model import XMLInputParam
from janis_core.ingestion.galaxy.gxtool.model import XMLOutputParam
from janis_core import settings


def extract_outputs(
    xmltool: XMLTool, 
    command: Command, 
    gxstep: Optional[dict[str, Any]]=None
    ) -> list[CommandComponent]:
    extractor = OutputExtractor(xmltool, command, gxstep)
    return extractor.extract()


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
        self.tool_outputs: dict[str, list[CommandComponent]] = defaultdict(list)

    def extract(self) -> list[CommandComponent]:
        outputs: list[CommandComponent] = []
        self.define_whitelisted_outputs()
        self.gather_input_outputs()
        self.gather_wildcard_outputs()
        self.gather_redirect_outputs()
        self.gather_uncertain_outputs()
        outputs = self.prioritise_outputs()
        return outputs

    def prioritise_outputs(self) -> list[CommandComponent]:
        # sorry this func is last minute and horrendous
        prioritised: list[CommandComponent] = []

        # set up data_structure so we can look up the possible outputs by name
        data_structure: dict[str, list[Tuple[str, CommandComponent]]] = defaultdict(list)
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
                output = factory.wildcard_output(gxparam)
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
                output = factory.input_output(component)
                self.tool_outputs['input'].append(output)
    
    def gather_uncertain_outputs(self) -> None:
        for gxparam in self.whitelisted_outputs:
            output = factory.uncertain_output(gxparam)
            self.tool_outputs['uncertain'].append(output)

    # CHECKS
    def should_create_input_output(self, component: CommandComponent) -> bool:
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
