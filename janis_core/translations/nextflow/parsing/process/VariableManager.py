

from typing import Optional, Any
from enum import Enum, auto
from dataclasses import dataclass, field

from janis_core import CommandTool, PythonTool
from janis_core.workflow.workflow import InputNode
from janis_core import translation_utils as utils

from ... import data_sources
from ...scope import Scope



class VariableType(Enum): 
    ProcessInput    = auto()
    ParamInput      = auto()
    Static          = auto()
    Ignored         = auto()
    Local           = auto()

@dataclass
class Variable:
    vtype: VariableType
    value: Optional[str | list[str]]

@dataclass
class VariableHistory:
    tinput_id: str
    items: list[Variable] = field(default_factory=list)

    @property
    def original(self) -> Variable:
        return self.items[0]
    
    @property
    def current(self) -> Variable:
        return self.items[-1]
    
    @property
    def all(self) -> list[Variable]:
        return self.items


class VariableManager:
    """
    within a process, for a given TInput, get the variable which currently represents that TInput.
    eg 

    input:
    1 path indexed_bam_flat
    2 path cath_crossmapped, stageAs: 'cath_crossmapped'

    script:
    3 def in_bams = get_primary_files(indexed_bam_flat, 2)
    4 def in_bams_joined = in_bams.join(' ')
    5 def cath_crossmapped = cath_crossmapped ? "-cx ${cath_crossmapped}" : ""
    '''
    6 echo \
    7 ${in_bams_joined} \
    8 ${cath_crossmapped} \
    9 --promoter=${params.promoter_bp} \
    '''

    TInput('inBams')
    line 1 -> indexed_bam_flat
    line 3 -> in_bams
    line 4 -> in_bams_joined
    line 6 -> in_bams_joined

    TInput('cathCrossmapped')
    line 2 -> cath_crossmapped
    line 8 -> cath_crossmapped

    TInput('promoterBp')
    line 9 -> params.promoter_bp

    data_structure description:
    {
        tinput_id: [
            first,  (original)
            second, 
            third,  (current)
            ..
        ]
    }

    """
    def __init__(self, scope: Scope) -> None:
        self.scope = scope
        self.data_structure: dict[str, VariableHistory] = {}
        self.process_inputs = data_sources.process_inputs(self.scope)
        self.param_inputs = data_sources.param_inputs(self.scope)
        self.internal_inputs = data_sources.internal_inputs(self.scope)

    def get(self, tinput_id: str) -> VariableHistory:
        return self.data_structure[tinput_id]
        
    def update(self, tinput_id: str, category: str, value: Optional[str | list[str]]) -> None:
        # create new variable history to track tinput_id if not exists
        if tinput_id not in self.data_structure:
            new_history = VariableHistory(tinput_id)
            self.data_structure[tinput_id] = new_history
        
        # gen new variable & update history
        history = self.data_structure[tinput_id]
        new_variable = self.gen_variable(category, value)
        history.items.append(new_variable)

    def update_for_tool(self, tool: CommandTool | PythonTool, sources: Optional[dict[str, Any]]=None) -> None:
        for tinput in tool.inputs():
            value = data_sources.get_variable(self.scope, tinput)
            # value is how the variable will be initially referred to inside a process
            # will be a process input, or param input.

            # if not fed via process input or param input, is a static input, or ignored input. 
            # check sources to see if there is an edge to an input node. if so, use its default. 
            if value is None:
                
                # workflow translations will use step.sources
                if sources:
                    value = self.get_value_from_input_node_default(tinput.id(), sources)
                    print()
                
                # tool translations will use tool.connections
                else:
                    value = self.get_value_from_input_node_default(tinput.id(), tool.connections)
                    print()
            
            category = self.get_category(tinput.id(), value)
            self.update(tinput.id(), category, value)

    def get_value_from_input_node_default(self, tinput_id: str, sources: dict[str, Any]) -> Any:
        if tinput_id in sources:
            src = sources[tinput_id]
            node = utils.resolve_node(src)
            if isinstance(node, InputNode):
                return node.default
        return None
    
    def get_value_from_connection(self, tinput_id: str, connections: dict[str, Any]) -> Any:
        if tinput_id in connections:
            return connections[tinput_id]
        return None

    def gen_variable(self, category: str, value: Optional[str | list[str]]) -> Variable:
        vtype_map = {
            'process': VariableType.ProcessInput,
            'param': VariableType.ParamInput,
            'static': VariableType.Static,
            'ignored': VariableType.Ignored,
            'local': VariableType.Local,
        }
        vtype = vtype_map[category]
        return Variable(vtype, value)

    def get_category(self, tinput_id: str, value: Optional[str | list[str]]) -> str:
        if tinput_id in self.process_inputs:
            return 'process'
        elif tinput_id in self.param_inputs:
            return 'param'
        elif tinput_id in self.internal_inputs:
            if value is None:
                return 'ignored'
            else:
                return 'static'
        else:
            return 'local'
    
    def to_string(self) -> str:
        out: str = ''
        for tinput_id, history in self.data_structure.items():
            block = ''
            block += f'{tinput_id}:\n'
            for var in history.items:
                block += f'[{var.vtype}] {var.value}\n'
            block += '\n'
            out += block
        return out
    


    # # helper methods
    # def _apply_index(self, var: Variable, index: Optional[int]) -> Variable:
    #     """
    #     for same variables, the var.value is actually a list of varnames.
        
    #     this will happen in the case of a ToolInput, where the type is a secondary file eg BamBai. 
    #     in the nextflow process, this bambai may appear like the following:
    #         inputs:
    #         tuple path(bam), path(bai)
        
    #     in this case, the ToolInput is split into the var.value ['bam', 'bai']
    #     to get the current var.value for the bam, we want var.value[0].
    #     to get the current var.value for the bai, we want var.value[1].

    #     """
    #     if isinstance(var.value, list):
    #         if index is not None:
    #             var.value = var.value[index]
    #         else:
    #             var.value = var.value[0]
    #     return var

