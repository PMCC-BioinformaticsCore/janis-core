

from typing import Any, Optional
import WDL
from collections import defaultdict
from dataclasses import dataclass
from copy import deepcopy


def get_entities_flat(workflow: WDL.Tree.Workflow) -> dict[str, list]:
    explorer = GraphExplorer(workflow)
    explorer.explore()
    return {
        'inputs': explorer.inputs,
        'calls': explorer.calls,
        'outputs': explorer.outputs
    }


@dataclass 
class CallContext:
    entity: WDL.Tree.Call
    scopedvars: list[Any]
    conditions: list[WDL.Tree.Conditional]
    scatter: Optional[WDL.Tree.Scatter]

@dataclass 
class OutputContext:
    entity: WDL.Tree.Decl 
    scopedvars: list[Any]

class GraphExplorer:
    def __init__(self, workflow: WDL.Tree.Workflow):
        self.workflow = workflow
        self.inputs: list[WDL.Tree.Decl] = []
        self.calls: list[CallContext] = []
        self.outputs: list[OutputContext] = []
        self.encountered_step = False

    def explore(self) -> None:
        self.explore_inputs()
        self.explore_calls()
        self.explore_outputs()

    def explore_inputs(self) -> None:
        if self.workflow.inputs is not None:
            for inp in self.workflow.inputs:
                self.inputs.append(inp)

    def explore_calls(self) -> None:
        if self.workflow.body is not None:
            self.explore_level(self.workflow.body, scopedvars=[], conditions=[])
        
    def explore_outputs(self) -> None:
        scopedvars = []
        if self.workflow.outputs is None:
            return None
        # tempvars
        for node in self.workflow.body:
            if isinstance(node, WDL.Tree.Decl):
                scopedvars.append(node)
        # actual outputs
        for out in self.workflow.outputs:
            flatout = OutputContext(out, scopedvars=scopedvars)
            self.outputs.append(flatout)
    
    def explore_level(
        self, 
        node: Any, 
        scopedvars: list[Any],
        conditions: list[WDL.Tree.Conditional], 
        scatter: Optional[WDL.Tree.Scatter]=None
        ) -> None:

        if not isinstance(node, list):
            node = [node]
        
        for item in node:
            if isinstance(item, WDL.Tree.Decl):
                scopedvars.append(item)
            elif isinstance(item, WDL.Tree.Conditional):
                this_vars = deepcopy(scopedvars)
                this_conds = deepcopy(conditions)
                this_conds.append(item)
                self.explore_level(item.body, this_vars, this_conds, scatter)
            elif isinstance(item, WDL.Tree.Scatter):
                this_vars = deepcopy(scopedvars)
                this_conds = deepcopy(conditions)
                self.explore_level(item.body, this_vars, this_conds, item)
            elif isinstance(item, WDL.Tree.Call):
                self.encountered_step = True
                call = CallContext(item, scopedvars, conditions, scatter)
                self.calls.append(call)