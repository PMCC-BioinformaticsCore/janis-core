

from typing import Optional
from enum import Enum, auto
from dataclasses import dataclass, field

from janis_core import CommandTool, PythonTool

from ... import data_sources
from ...data_sources import DataSourceType
from ...scope import Scope



class VariableType(Enum): 
    PROCESS_INPUT   = auto()
    PARAM           = auto()
    STATIC          = auto()
    IGNORED         = auto()
    LOCAL           = auto()

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
        self.process_inputs = data_sources.task_inputs(self.scope)
        self.param_inputs = data_sources.param_inputs(self.scope)
        self.internal_inputs = data_sources.internal_inputs(self.scope)

    def get(self, tinput_id: str) -> VariableHistory:
        return self.data_structure[tinput_id]
        
    def update(self, tinput_id: str, vtype: VariableType, value: Optional[str | list[str]]) -> None:
        # create new variable history to track tinput_id if not exists
        if tinput_id not in self.data_structure:
            new_history = VariableHistory(tinput_id)
            self.data_structure[tinput_id] = new_history
        
        # gen new variable & update history
        history = self.data_structure[tinput_id]
        new_variable = Variable(vtype, value)
        history.items.append(new_variable)

    def update_for_tool(self, tool: CommandTool | PythonTool) -> None:
        for tinput in tool.inputs():
            ds = data_sources.get(self.scope, tinput)
            vtype = self.cast_to_variable_type(ds.dstype)
            self.update(tinput.id(), vtype, ds.value)

    def cast_to_variable_type(self, dstype: DataSourceType) -> VariableType:
        if dstype == DataSourceType.TASK_INPUT:
            return VariableType.PROCESS_INPUT
        elif dstype == DataSourceType.PARAM:
            return VariableType.PARAM
        elif dstype == DataSourceType.STATIC:
            return VariableType.STATIC
        elif dstype == DataSourceType.IGNORED:
            return VariableType.IGNORED
        else:
            raise RuntimeError
    
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
    
