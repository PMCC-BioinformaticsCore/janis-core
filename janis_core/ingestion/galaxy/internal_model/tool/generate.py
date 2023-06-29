

from typing import Optional, Any, Tuple
from collections import defaultdict
from dataclasses import dataclass

from janis_core import settings
from janis_core.ingestion.galaxy.gxtool.command.cmdstr.analysis import CmdstrReferenceType
from janis_core.ingestion.galaxy.gxtool.command.cmdstr.analysis import get_cmdstr_appearences
from janis_core.ingestion.galaxy.gxtool.command import Command
from janis_core.ingestion.galaxy.gxtool.command.components import InputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import OutputComponent
from janis_core.ingestion.galaxy.gxtool.command.components import RedirectOutput
from janis_core.ingestion.galaxy.gxtool.command.components import Flag
from janis_core.ingestion.galaxy.gxtool.command.components import Option
from janis_core.ingestion.galaxy.gxtool.command.components import Positional
from janis_core.ingestion.galaxy.gxtool.command.components import factory as component_factory
from janis_core.ingestion.galaxy.gxworkflow import load_tool_state
from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.gxtool.model import (
    XMLParam, 
    XMLInputParam,
    XMLDataParam,
    XMLDataCollectionParam,
    XMLBoolParam,
    XMLOutputParam,
    XMLTextParam
)

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
    

@dataclass
class ComponentMetrics:
    component: InputComponent

    @property
    def score(self) -> int:
        score = 0
        score += 1 if self.component.confidence.value == 3 else 0
        score += 1 if self.component.gxparam is not None else 0
        return score


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
        self.add_captured_inputs()
        self.add_uncaptured_inputs()
        # self.resolve_duplicates()
        return self.inputs

    def add_captured_inputs(self) -> None:
        # inputs we identified in the tool command
        self.command.set_cmd_positions()
        self.inputs += self.command.list_inputs(include_base_cmd=False)
    
    def add_uncaptured_inputs(self) -> None:
        # inputs we identified in the tool command
        if not self.gxstep:
            return None
        
        unlinked_params: list[XMLParam] = []
        uncaptured_inputs: list[InputComponent] = []

        tool_state = load_tool_state(self.gxstep, additional_filters=['Flatten', 'DeNestClass'])
        for pname, pvalue in tool_state.items():
            param = self.xmltool.inputs.get(pname)
            if param:
                if not self.param_unlinked(param):
                    continue
                if not self.should_create(param, pvalue):
                    continue

                unlinked_params.append(param)
                component = self.create_uncaptured_input(param, pvalue)
                uncaptured_inputs.append(component)

        if unlinked_params:
            print('\n--- UNLINKED PARAMS ---')
            for param in unlinked_params:
                print(f'{param.name}: {param.__class__.__name__}')

        if uncaptured_inputs:
            print('\n--- UNCAPTURED INPUTS ---')
            for component in uncaptured_inputs:
                assert(component.gxparam)
                print(f'{component.__class__.__name__}: {component.gxparam.name}')

        if uncaptured_inputs:
            if self.can_link_uncaptured_positionals(uncaptured_inputs):
                uncaptured_inputs = self.link_uncaptured_positionals(uncaptured_inputs)
        
        for component in uncaptured_inputs:
            self.inputs.append(component)

    def should_create(self, param: XMLParam, pvalue: Any) -> bool:
        # data
        if isinstance(param, XMLDataParam) or isinstance(param, XMLDataCollectionParam):
            return True
        
        # text
        elif isinstance(param, XMLTextParam):
            BLANK_TEXTPARAM_VALUES = set(['', 'null'])
            if pvalue not in BLANK_TEXTPARAM_VALUES:
                return True
        
        # all other types
        elif self.param_is_used(param):
            return True

        return False

    def param_unlinked(self, param: XMLParam) -> bool:
        for inp in self.inputs:
            if inp.gxparam and inp.gxparam.name == param.name:
                return False
        return True

    def param_is_used(self, param: XMLParam) -> bool:
        # ignore if only appears in <command> section control structures, ignore
        appearences = get_cmdstr_appearences(
            self.xmltool.raw_command, 
            param, 
            filter_to=[
                CmdstrReferenceType.INLINE_CHEETAH_LOOP,
                CmdstrReferenceType.INLINE_PLAIN_TEXT
            ])
        if not appearences:
            return False
        return True
    
    def create_uncaptured_input(self, gxparam: XMLParam, pvalue: Any) -> InputComponent:
        if isinstance(gxparam, XMLBoolParam):
            component = self.create_component_for_boolparam(gxparam)
        elif gxparam.argument:  # this should never happen as ArgumentAnnotator should have identified it
            component = component_factory.option(prefix=gxparam.argument, gxparam=gxparam)
            if pvalue is not None:
                component.values.add(pvalue)
        else:
            component = component_factory.positional(gxparam=gxparam)
            if pvalue is not None:
                component.values.add(pvalue)
        return component
        
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
        
    def can_link_uncaptured_positionals(self, uncaptured_inputs: list[InputComponent]) -> bool:
        unlinked_positionals = [x for x in self.inputs if isinstance(x, Positional) and not x.gxparam]
        uncaptured_positionals = [x for x in uncaptured_inputs if isinstance(x, Positional)]
        if len(uncaptured_positionals) == len(unlinked_positionals):
            print('\nLINKING UNCAPTURED POSITIONALS POSSIBLE!')
            return True
        return False

    def link_uncaptured_positionals(self, uncaptured_inputs: list[InputComponent]) -> list[InputComponent]:
        unlinked_positionals = [x for x in self.inputs if isinstance(x, Positional) and not x.gxparam]
        uncaptured_positionals = [x for x in uncaptured_inputs if isinstance(x, Positional)]
        uncaptured_nonpositionals = [x for x in uncaptured_inputs if not isinstance(x, Positional)]
        
        if self.can_link_via_numbering(uncaptured_positionals, unlinked_positionals):
            self.do_positional_link_via_numbering(uncaptured_positionals, unlinked_positionals)
        else:
            self.do_positional_link_arbitrary(uncaptured_positionals, unlinked_positionals)
        
        return uncaptured_nonpositionals

    def can_link_via_numbering(self, uncaptured_positionals: list[Positional], unlinked_positionals: list[Positional]) -> bool:
        numbers = list(range(1, len(uncaptured_positionals) + 1))
        for num in numbers:
            if sum([x.name.endswith(str(num)) for x in uncaptured_positionals]) != 1:
                if sum([x.name.endswith(str(num)) for x in unlinked_positionals]) != 1:
                    return False
        return True

    def do_positional_link_via_numbering(self, uncaptured_positionals: list[Positional], unlinked_positionals: list[Positional]) -> None:
        numbers = list(range(1, len(uncaptured_positionals) + 1))
        for num in numbers:
            uncaptured = [x for x in uncaptured_positionals if x.name.endswith(str(num))][0]
            unlinked = [x for x in unlinked_positionals if x.name.endswith(str(num))][0]
            unlinked.gxparam = uncaptured.gxparam

    def do_positional_link_arbitrary(self, uncaptured_positionals: list[Positional], unlinked_positionals: list[Positional]) -> None:
        for i in range(len(uncaptured_positionals)):
            unlinked_positionals[i].gxparam = uncaptured_positionals[i].gxparam
        


    ### DEPRECATED

    def resolve_duplicates(self) -> None:
        """
        resolves situations where 2 components were extracted for seemingly the same CLI argument.
        an example is where a flag was identified with prefix = '-F' an option with the same prefix was also identified.
        TODO this should probably happen automatically in the CmdstrCommandAnnotator class.
        benefit is that we could be more aware of the structure of the command and avoid this situation.
        """
        self.resolve_duplicate_gxparam_components()
        self.resolve_duplicate_prefix_components()

    def resolve_duplicate_gxparam_components(self) -> None:
        # TODO i can't remember if this can happen, but if it can, we should resolve it.
        gxparam_dict = defaultdict(list) 
        for comp in self.inputs:
            if comp.gxparam:
                gxparam_dict[comp.gxparam.name].append(comp)
        
        duplicates = True if any([len(components) > 1 for components in gxparam_dict.values()]) else False
        if duplicates:
            print('\n---DUPLICATE GXPARAM ---')
            for gxparam, components in gxparam_dict.items():
                if len(components) > 1:
                    for component in components:
                        if isinstance(component, Flag):
                            details = f'[flag] {component.prefix}'
                        elif isinstance(component, Option):
                            details = f'[option] {component.prefix}'
                        elif isinstance(component, Positional):
                            details = f'[positional] {component.cmd_pos}'
                        elif isinstance(component, OutputComponent):
                            details = f'[output] {component.name}'
                        print(f'{gxparam}: {details}')

    def resolve_duplicate_prefix_components(self) -> None:
        flags = set([x.prefix for x in self.inputs if isinstance(x, Flag)])
        options = set([x.prefix for x in self.inputs if isinstance(x, Option)])
        intersection = flags & options
        if intersection:
            print('\n--- DUPLICATE FLAG / OPTION PREFIXES ---')
            for prefix in intersection:
                print(prefix) 

        # whitelisted_uuids = [x.uuid for x in self.inputs]
        
        # # generate data structure
        # prefix_dict = defaultdict(list) 
        # for comp in self.inputs:
        #     if isinstance(comp, Flag | Option): 
        #         prefix_dict[comp.prefix].append(comp)
        #     else:
        #         # ignore positionals - automatically whitelisted
        #         whitelisted_uuids.append(comp.uuid)
        
        # # resolve
        # for prefix, components in prefix_dict.items():
        #     if len(components) == 0:
        #         raise RuntimeError('should not happen - debugging')
        #     elif len(components) == 1:
        #         component = components[0]
        #         whitelisted_uuids.append(component.uuid)
        #     elif len(components) >= 2:
        #         component = self.select_best_component(components)
        #         whitelisted_uuids.append(component.uuid)

        # self.inputs = [x for x in self.inputs if x.uuid in whitelisted_uuids]

    # def select_best_component(self, components: list[InputComponent]) -> InputComponent:
    #     components_metrics = [ComponentMetrics(x) for x in components]
    #     components_metrics.sort(key=lambda x: x.score, reverse=True)

    #     # select clear best 
    #     if components_metrics[0].score > components_metrics[1].score:
    #         return components_metrics[0].component
        
    #     # remove those which are not equal best
    #     components = [x.component for x in components_metrics if x.score == components_metrics[0].score]
        
    #     # if all have gxparam or all have high confidence, select equal best via type
    #     if all([x.confidence.value == 3 for x in components]) or all([x.gxparam is not None for x in components]):
    #         if any([isinstance(x, Option) for x in components]):
    #             return [x for x in components if isinstance(x, Option)][0]
    #         else:
    #             return components[0]
        
    #     # if some have gxparam and some have high confidence, select best via confidence
    #     elif any([x.confidence.value == 3 for x in components]) and any([x.gxparam is not None for x in components]):
    #         return [x for x in components if x.confidence.value == 3][0]
        
    #     # if none have high confidence or gxparam, select best via type
    #     elif any([isinstance(x, Option) for x in components]):
    #         return [x for x in components if isinstance(x, Option)][0]
        
    #     # all options with low confidence & no gxparam
    #     else:
    #         return components[0]
        


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

    

