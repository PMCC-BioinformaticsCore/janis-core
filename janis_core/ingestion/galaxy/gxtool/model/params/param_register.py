

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

from .param import XMLParam


class XMLParamRegister:
    def __init__(self):
        self.params: list[XMLParam] = []

    def list(self) -> list[XMLParam]:
        return self.params
    
    def add(self, param: XMLParam) -> None:
        """adds a param to register. enforces unique param var names"""
        self.params.append(param)

    def get(self, query: str, strategy: str='exact') -> Optional[XMLParam]:
        """performs search using the specified search strategy"""
        strategy_map = {
            'exact': ExactSearchStrategy(),
            'lca': LCASearchStrategy(),
            'filepath': FilepathSearchStrategy(),
        }

        search_strategy = strategy_map[strategy]
        return search_strategy.search(query, self.params)


class SearchStrategy(ABC):    
    @abstractmethod
    def search(self, query: str, params: list[XMLParam]) -> Optional[XMLParam]:
        """searches for a param using some concrete strategy"""
        ...

class ExactSearchStrategy(SearchStrategy):
    def search(self, query: str, params: list[XMLParam]) -> Optional[XMLParam]:
        """searches for a param using param name"""
        for param in params:
            if param.name == query:
                return param
        return None

@dataclass
class LCAParam:
    split_name: list[str]
    param: XMLParam

class LCASearchStrategy(SearchStrategy):
    def search(self, query: str, params: list[XMLParam]) -> Optional[XMLParam]:
        """searches for a param using LCA"""
        split_query = query.split('.')
        remaining_params = self.init_datastructure(params)

        for i in range(1, len(split_query) + 1):
            remaining_params = [p for p in remaining_params if len(p.split_name) >= i]
            remaining_params = [p for p in remaining_params if p.split_name[-i] == split_query[-i]]

        if len(remaining_params) > 0:
            # return shortest param name which matches the query
            remaining_params.sort(key=lambda x: len(x.split_name))
            return remaining_params[0].param

    def init_datastructure(self, params: list[XMLParam]):
        return [LCAParam(param.name.split('.'), param) for param in params]
        

class FilepathSearchStrategy(SearchStrategy):
    def search(self, query: str, params: list[XMLParam]) -> Optional[XMLParam]:
        """
        searches for a param by matching the specified 
        from_work_dir path to a given filepath
        """
        raise NotImplementedError
        # for param in params.values():
        #     if hasattr(param, 'from_work_dir') and param.from_work_dir == query:
        #         return param
