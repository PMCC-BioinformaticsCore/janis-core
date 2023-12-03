
import WDL 
from typing import Any, Optional

from janis_core import WorkflowBuilder, ScatterDescription
from janis_core.workflow.workflow import InputNode, OutputNode, StepNode
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory

from ..types import parse_type
from ..expressions import parse_expr
from ..EntityParser import WorkflowParser
from .explore import InvertedCall


class WorkflowInputParser(WorkflowParser):

    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder, wdl_inp: WDL.Tree.Decl) -> None:
        super().__init__(wdl_wf, janis_wf)
        self.wdl_inp = wdl_inp

    def do_parse(self) -> InputNode:
        default = None
        if self.wdl_inp.expr:
            default = parse_expr(self.wdl_inp.expr, self.wdl_wf, self.janis_wf)
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

    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder, wdl_src: Any) -> None:
        super().__init__(wdl_wf, janis_wf)
        self.wdl_src = wdl_src
    
    def do_parse(self) -> Any:
       return parse_expr(self.wdl_src, self.wdl_wf, self.janis_wf)
    
    def fallback(self) -> Any:
        raise NotImplementedError




class WorkflowStepModifierParser(WorkflowParser):
    
    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder, inv_call: InvertedCall) -> None:
        super().__init__(wdl_wf, janis_wf)
        self.inv_call = inv_call
        self.call = inv_call.call
        self.dependencies = inv_call.dependencies
        self.scatter = inv_call.scatter
        self.jstep = self.janis_wf[self.call.name]

    def do_parse(self) -> None:
        # what about scatter & when? 
        self.scatter = self.parse_scatter()
        self.when = self.parse_dependencies()
    
    def parse_scatter(self) -> Optional[ScatterDescription]:
        if not self.scatter:
            return None
        
        raise NotImplementedError

        foreach = parse_expr(self.call.expr, self.wdl_wf, self.janis_wf)

        # TODO this uuid should be the janis step uuid. 
        scar_var_type = parse_type(self.call.expr.type, self.wdl_wf, uuid=self.janis_wf.uuid)
        if isinstance(scar_var_type, WDL.Type.Array):
            scar_var_type = scar_var_type.item_type

        # when we unwrap each step-input to the workflow, we want to replace 'self.call.variable' with
        #       lambda el: <operation with self.call.variable substituted for {el}>
        # if self.call.variable not in wf.input_nodes:
        #     wf.input(self.call.variable, scar_var_type)
        for inner_call in self.call.body:
            self.add_call_to_wf(
                wf, inner_call, foreach=foreach, expr_alias=self.call.variable
            )
    
    def parse_dependencies(self) -> Any:
        if not self.dependencies:
            return None
        
        raise NotImplementedError
        
        # TODO add j.ForEachSelector() to parse_expr in this case if list?
        for inner_call in self.call.body:

            # inner_call = call.body[0]
            self.add_call_to_wf(
                wf,
                inner_call,
                condition=self.translate_expr(
                    call.expr, input_selector_getter=selector_getter
                ),
                expr_alias=expr_alias,
                foreach=foreach,
            )

    def fallback(self) -> None:
        raise NotImplementedError



class WorkflowOutputParser(WorkflowParser):

    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder, wdl_out: WDL.Tree.Decl) -> None:
        super().__init__(wdl_wf, janis_wf)
        self.wdl_out = wdl_out

    def do_parse(self) -> OutputNode:
        if self.wdl_out.expr is None:
            raise Exception(f"Output {self.wdl_out.name} has no expression")
        # TODO UPDATE EXPRESSION PARSING FOR WORKFLOW SCOPE
        dtype = parse_type(self.wdl_out.type, self.wdl_wf, uuid=self.janis_wf.uuid)
        sel = parse_expr(self.wdl_out.expr, self.wdl_wf, self.janis_wf)
        return self.janis_wf.output(self.wdl_out.name, dtype, sel)
    
    def fallback(self) -> OutputNode:
        raise NotImplementedError



