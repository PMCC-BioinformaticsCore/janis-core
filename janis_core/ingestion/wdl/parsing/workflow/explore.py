

from typing import Any, Optional
import WDL
from collections import defaultdict
from dataclasses import dataclass
from copy import deepcopy

def get_entities(workflow: WDL.Tree.Workflow) -> dict[str, list]:
    explorer = GraphExplorer(workflow)
    explorer.explore()
    return explorer.entities

@dataclass 
class InvertedCall:
    call: WDL.Tree.Call
    dependencies: list[WDL.Tree.Conditional]
    scatter: Optional[WDL.Tree.Scatter]

class GraphExplorer:
    def __init__(self, workflow: WDL.Tree.Workflow):
        self.workflow = workflow
        self.entities = defaultdict(list)

    def explore(self) -> None:
        if self.workflow.inputs is not None:
            for inp in self.workflow.inputs:
                self.entities['input'].append(inp)
        if self.workflow.body is not None:
            for node in self.workflow.body:
                self.explore_node(node, deps=[])
        if self.workflow.outputs is not None:
            for out in self.workflow.outputs:
                self.entities['output'].append(out)
    
    def explore_node(self, node: Any, deps: list[WDL.Tree.Conditional], scatter: Optional[WDL.Tree.Scatter]=None) -> None:
        if isinstance(node, list):
            # pass deps and scatter down
            this_deps = deepcopy(deps)
            for item in node:
                self.explore_node(item, this_deps, scatter=scatter)
        
        # ignoring dependencies for inputs
        elif isinstance(node, WDL.Tree.Decl):
            self.entities['input'].append(node)

        elif isinstance(node, WDL.Tree.Call):
            call = InvertedCall(call=node, dependencies=deps, scatter=scatter)
            self.entities['call'].append(call)
        
        elif isinstance(node, WDL.Tree.Scatter):
            # pass deps down with new scatter
            this_deps = deepcopy(deps)
            self.explore_node(node.body, this_deps, scatter=node)
        
        elif isinstance(node, WDL.Tree.Conditional):
            # pass scatter down with additional deps
            this_deps = deepcopy(deps)
            this_deps += [node]
            self.explore_node(node.body, this_deps, scatter=scatter)
        
        # don't explore children
        else:
            return None
      