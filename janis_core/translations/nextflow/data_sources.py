

from typing import Optional
from dataclasses import dataclass, field

from janis_core import (
    ToolInput,
    TInput,
)
from .scope import Scope
from enum import Enum, auto


class DataSourceType(Enum): 
    TASK_INPUT  = auto()
    PARAM       = auto()
    STATIC      = auto()
    IGNORED     = auto()


@dataclass 
class DataSource:
    tinput_id: str
    value: Optional[str | list[str]]
    dstype: DataSourceType


@dataclass
class TaskDSVariableRegister:
    
    data_structure: dict[str, dict[str, DataSource]] = field(default_factory=dict)
    
    """ 
    for each Workflow / CommandTool / PythonTool, stores a map of tinput_id: data_source.
    the data_source is the initial value which drives the ToolInput inside a task. 
    
    data_structure description:
    {
        Scope: {
            'my_file': str,
            'my_string': None,
            'my_secondary': list[str],
            'my_secondary_arr': str,
        }
    }
    """
    
    def add(
        self, 
        scope: Scope,
        dstype_str: str,
        tinput_id: str, 
        value: Optional[str | list[str]],

    ) -> None:
        """
        adds some data to our data_structure. 
        scope: tool scope
        
        """
        label = scope.to_string()

        # if scope not yet in data structure, add
        if label not in self.data_structure:
            self.data_structure[label] = {}

        # cast the data source subtype to enum value
        subtype_map = {
            'task_input': DataSourceType.TASK_INPUT,
            'param': DataSourceType.PARAM,
            'static': DataSourceType.STATIC,
            'ignored': DataSourceType.IGNORED,
        }

        # create new DataSource & add to data_structure
        dstype = subtype_map[dstype_str]
        ds = DataSource(tinput_id, value, dstype)
        self.data_structure[label][tinput_id] = ds

    def get(self, scope: Scope, tinput_id: str) -> DataSource:
        label = scope.to_string()
        return self.data_structure[label][tinput_id]
    
    def getall(self, scope: Scope) -> list[DataSource]:
        label = scope.to_string()
        return list(self.data_structure[label].values())
    
    def to_string(self) -> str:
        out: str = ''
        for scope_label, name_map in self.data_structure.items():
            out += f'\n{scope_label}\n'
            for tinput_id, ds in name_map.items():
                out += f'{tinput_id}: {ds.value}\n'
        return out

    

tvn_register = TaskDSVariableRegister()

def get(scope: Scope, inp: ToolInput | TInput) -> DataSource:
    return tvn_register.get(scope, inp.id())

def update(scope: Scope, dstype_str: str, tinput_id: str, value: Optional[str | list[str]]):
    tvn_register.add(scope, dstype_str, tinput_id, value)

def task_inputs(scope: Scope) -> set[str]:
    scoped_data_sources = tvn_register.getall(scope)
    return set([x.tinput_id for x in scoped_data_sources if x.dstype == DataSourceType.TASK_INPUT])

def param_inputs(scope: Scope) -> set[str]:
    scoped_data_sources = tvn_register.getall(scope)
    return set([x.tinput_id for x in scoped_data_sources if x.dstype == DataSourceType.PARAM])

def internal_inputs(scope: Scope) -> set[str]:
    scoped_data_sources = tvn_register.getall(scope)
    out: set[str] = set()
    for ds in scoped_data_sources:
        if ds.dstype in [DataSourceType.STATIC, DataSourceType.IGNORED]:
            out.add(ds.tinput_id)
    return out

def clear() -> None:
    tvn_register.data_structure = {}




    
# tds_register = TaskDSCategoryRegister()

# # categories
# def update_categories(scope: Scope, process_inputs: set[str], param_inputs: set[str], internal_inputs: set[str]):
#     tds_register.add(scope, 'process', process_inputs)
#     tds_register.add(scope, 'param', param_inputs)
#     tds_register.add(scope, 'internal', internal_inputs)


# @dataclass
# class TaskDSCategoryRegister:
    
#     data_structure: dict[str, dict[str, set[str]]] = field(default_factory=dict)
    
#     """
#     for each Workflow / CommandTool / PythonTool, for each ToolInput, stores data on whether the ToolInput is fed value via:
#         - a process input
#         - a global param
#         - not fed a value (internal input - static value)

#     data_structure description:
#     {
#         Scope: {
#             'process': set[str],
#             'param': set[str],
#             'internal': set[str]
#         }
#     }
#     """
    
#     def add(self, scope: Scope, subtype: str, tinput_ids: set[str]) -> None:
#         """
#         adds some data to our data_structure. 
#         scope: tool scope
#         subtype: one of 'process', 'param', 'internal'
#         tinput_ids: set of ToolInput identifiers belonging to that subtype
#         """
#         label = scope.to_string()
#         if label not in self.data_structure:
#             self.data_structure[label] = {}
#         self.data_structure[label][subtype] = tinput_ids

#     def get(self, scope: Scope, subtype: str) -> set[str]:
#         label = scope.to_string()
#         return self.data_structure[label][subtype]

#     def to_string(self) -> str:
#         out: str = ''
#         for scope_label, categories in self.data_structure.items():
#             out += f'\n{scope_label}\n'
#             for catname, tinput_ids in categories.items():
#                 out += f'{catname}: {" ".join(tinput_ids)}\n'
#         return out
    