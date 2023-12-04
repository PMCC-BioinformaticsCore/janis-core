
import WDL 
from typing import Any, Optional
from copy import deepcopy

from janis_core import WorkflowBuilder, ScatterDescription
from janis_core.workflow.workflow import InputNode, OutputNode, StepNode
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory
from janis_core import ScatterMethod
from janis_core import ScatterDescription
from janis_core import AndOperator

from ..types import parse_type
from ..expressions import parse_expr
# from ..expressions import downstream_nodes
from ..expressions import expr_as_str
from ..EntityParser import WorkflowParser
from .explore import CallContext
from .explore import OutputContext





class WorkflowInputParser(WorkflowParser):

    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder, wdl_inp: WDL.Tree.Decl) -> None:
        super().__init__(wdl_wf, janis_wf)
        self.wdl_inp = wdl_inp

    def do_parse(self) -> InputNode:
        default = None
        if self.wdl_inp.expr:
            default, success = parse_expr(self.wdl_inp.expr, self.wdl_wf, self.janis_wf)
        dtype = parse_type(self.wdl_inp.type, self.wdl_wf, uuid=self.janis_wf.uuid)
        
        selector = self.janis_wf.input(
            self.wdl_inp.name, 
            dtype, 
            default=default
        )
        return selector.input_node
    
    def fallback(self) -> InputNode:
        # assumes issue is from parsing input default expression
        msg = f'Error parsing workflow input: {self.wdl_inp.name}'
        log_message(self.janis_wf.uuid, msg, category=ErrorCategory.FALLBACKS)
        dtype = parse_type(self.wdl_inp.type, self.wdl_wf, uuid=self.janis_wf.uuid)
        selector = self.janis_wf.input(
            self.wdl_inp.name, 
            dtype, 
            default=None
        )
        return selector.input_node




class WorkflowStepInputParser(WorkflowParser):

    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder, wdl_src: Any, flatcall: CallContext) -> None:
        super().__init__(wdl_wf, janis_wf)
        self.wdl_src = wdl_src
        self.flatcall = flatcall
    
    def do_parse(self) -> Any:
        src, success = parse_expr(self.wdl_src, self.wdl_wf, self.janis_wf, self.flatcall)
        return src
    
    def fallback(self) -> Any:
        raise NotImplementedError





def get_dependent_inputs(node: Any, call: WDL.Tree.Call) -> set[str]:
    tracer = StepScatterTracer(call)
    tracer.trace(node, [])
    return tracer.entities

class StepScatterTracer:
    def __init__(self, call: WDL.Tree.Call) -> None:
        self.call = call
        self.entities: set[str] = set()

    def trace(self, node: Any, idents: list[str]) -> None:
        if isinstance(node, WDL.Tree.Scatter):
            this_idents = deepcopy(idents)
            this_idents += [node.variable]
            for child in node.body:
                if isinstance(child, WDL.Tree.Decl):
                    for ident in this_idents:
                        if self.references_ident(child, ident):
                            this_idents.append(child.name)
                else:
                    self.trace(child, this_idents)

        elif isinstance(node, WDL.Tree.Conditional):
            this_idents = deepcopy(idents)
            if isinstance(child, WDL.Tree.Decl):
                for ident in this_idents:
                    if self.references_ident(child, ident):
                        this_idents.append(child.name)
            else:
                self.trace(child, this_idents)
        
        # found correct call
        elif isinstance(node, WDL.Tree.Call) and node.name == self.call.name:
            for inpname, source in node.inputs.items():
                for ident in idents:
                    if self.references_ident(source, ident):
                        self.entities.add(inpname)
        
        else:
            raise NotImplementedError

    def references_ident(self, query: Any, target: str) -> bool:
        if isinstance(query, WDL.Expr.Ident) and str(query.name) == target:
            return True 
        for child in query.children:
            if self.references_ident(child, target):
                return True
        return False

    

class WorkflowStepModifierParser(WorkflowParser):
    
    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder, flatcall: CallContext) -> None:
        super().__init__(wdl_wf, janis_wf)
        self.call = flatcall.entity
        self.conditions = flatcall.conditions
        self.scatter = flatcall.scatter
        self.jstep = self.janis_wf[self.call.name]

    def do_parse(self) -> None:
        self.scatter = self.parse_scatter()
        self.when = self.parse_conditional()

    def parse_scatter(self) -> Optional[ScatterDescription]:
        if not self.scatter:
            return None
        method = self.parse_scatter_method()
        fields = self.parse_scatter_fields()
        return ScatterDescription(fields, method=method)
        
    def parse_scatter_method(self) -> ScatterMethod:
        assert self.scatter is not None
        expr_str = expr_as_str(self.scatter.expr)
        if 'cross(' in expr_str:
            msg = f"Cross scatter detected: {{{expr_str}}}"
            log_message(self.jstep.uuid, msg, category=ErrorCategory.PLUMBING)
            return ScatterMethod.cross
        return ScatterMethod.dot

    def parse_scatter_fields(self) -> list[str]:
        assert self.scatter is not None
        fields = get_dependent_inputs(self.scatter, self.call)
        return list(fields)

    def parse_conditional(self) -> Any:
        if len(self.conditions) == 0:
            return None
        
        if len(self.conditions) == 1:
            res, success = parse_expr(self.conditions[0].expr, self.wdl_wf, self.janis_wf)
            return res 
        else:
            cond_res, success = parse_expr(self.conditions[0].expr, self.wdl_wf, self.janis_wf)
            for cond in self.conditions[1:]:
                res, success = parse_expr(cond.expr, self.wdl_wf, self.janis_wf)
                cond_res = AndOperator(cond_res, res)
            return cond_res

    def fallback(self) -> None:
        raise NotImplementedError
    



class WorkflowOutputParser(WorkflowParser):

    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder, flatout: OutputContext) -> None:
        super().__init__(wdl_wf, janis_wf)
        self.flatout = flatout

    def do_parse(self) -> OutputNode:
        wdl_out = self.flatout.entity
        if wdl_out.expr is None:
            raise Exception(f"Output {wdl_out.name} has no expression")
        # TODO UPDATE EXPRESSION PARSING FOR WORKFLOW SCOPE
        dtype = parse_type(wdl_out.type, self.wdl_wf, uuid=self.janis_wf.uuid)
        sel, success = parse_expr(wdl_out.expr, self.wdl_wf, self.janis_wf, self.flatout)

        # lmao janis doens't allow this but whatever
        identifier = wdl_out.name
        otp = OutputNode(
            self,
            identifier=identifier,
            datatype=dtype,
            source=sel,
            skip_typecheck=True
        )
        self.janis_wf.nodes[identifier] = otp
        self.janis_wf.output_nodes[identifier] = otp
        return otp
    
    def fallback(self) -> OutputNode:
        raise NotImplementedError



