
import WDL 
from collections import defaultdict
from typing import Any


def get_decl_refs(task: WDL.Tree.Task, decls: set[str]) -> dict[str, dict]:
    tracer = DeclReferenceTracer(task, decls)
    return tracer.trace()

class DeclReferenceTracer:
    """
    NOTE: this is incomplete. Does not check for indirect refs via scopedvar. 
    """
    def __init__(self, task: WDL.Tree.Task, decls: set[str]) -> None:
        self.task = task
        self.decls = decls
        self.references = defaultdict(dict)

    def trace(self) -> dict[str, dict]:
        self.init_references()
        self.trace_inputs()
        self.trace_scopedvars()
        self.trace_command()
        self.trace_outputs()
        self.trace_runtime()
        return self.references
    
    def init_references(self) -> None:
        for decl in self.decls:
            self.references[decl]['inputs'] = 0
            self.references[decl]['scopedvars'] = 0
            self.references[decl]['command'] = 0
            self.references[decl]['outputs'] = 0
            self.references[decl]['runtime'] = 0

    def trace_inputs(self) -> None:
        self.context = 'inputs'
        if self.task.inputs is None:
            return 
        for inp in self.task.inputs:
            if inp.expr is not None:
                self.explore(inp.expr)
    
    def trace_scopedvars(self) -> None:
        self.context = 'scopedvars'
        if self.task.postinputs is None:
            return 
        for node in self.task.postinputs:
            if isinstance(node, WDL.Tree.Decl):
                if node.expr is not None:
                    self.explore(node.expr)
    
    def trace_command(self) -> None:
        self.context = 'command'
        for child in self.task.command.children:
            if isinstance(child, WDL.Expr.Base):
                self.explore(child)
    
    def trace_outputs(self) -> None:
        self.context = 'outputs'
        if self.task.outputs is None:
            return 
        for out in self.task.outputs:
            if out.expr is not None:
                self.explore(out.expr)
    
    def trace_runtime(self) -> None:
        self.context = 'runtime'
        for value in self.task.runtime.values():
            if isinstance(value, WDL.Expr.Base):
                self.explore(value)

    def explore(self, expr: WDL.Expr.Base) -> None:
        if isinstance(expr, WDL.Expr.Ident) and str(expr.name) in self.decls:
            self.references[str(expr.name)][self.context] += 1
        
        for child in expr.children:
            if isinstance(child, WDL.Expr.Base):
                self.explore(child)