#!/usr/bin/env python3

import functools
import os
import WDL
from typing import Any

from .parsing.workflow.explore import get_entities_flat
from janis_core import WorkflowBuilder, CommandToolBuilder
from janis_core.workflow.workflow import InputNode, StepNode, OutputNode

from janis_core.ingestion.common.graph import add_step_edges_to_graph
from .parsing import WorkflowInputParser
from .parsing import WorkflowStepModifierParser
from .parsing import WorkflowStepInputParser
from .parsing import WorkflowOutputParser

def error_boundary(return_value=None):
    def try_catch_translate_inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if not WdlParser.allow_errors:
                return func(*args, **kwargs)
            else:
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    j.Logger.log_ex(e)
                    return return_value

        return wrapper

    return try_catch_translate_inner


class WdlParser:

    allow_errors = False

    @staticmethod
    def from_doc(doc: str, base_uri=None):
        abs_path = os.path.relpath(doc)
        d = WDL.load(abs_path)

        parser = WdlParser()
        if d.workflow:
            return parser.ingest(d.workflow)
        tasks = []
        for t in d.tasks:
            tasks.append(parser.ingest(t))

        return tasks[0]

    def ingest(self, obj: WDL.SourceNode):
        """Main ingest entry point"""
        if isinstance(obj, WDL.Tree.Task):
            return self.ingest_task(obj)
        elif isinstance(obj, WDL.Tree.Workflow):
            return self.ingest_workflow(obj)
        else:
            raise RuntimeError(f"Unhandled WDL object type: {type(obj)}")
        
    def ingest_task(self, obj) -> CommandToolBuilder:
        """Tool ingest entry point"""
        from .parsing import parse_task
        return parse_task(obj) 
    
    def ingest_workflow(self, wdl_wf) -> WorkflowBuilder:
        """Workflow ingest entry point"""
        entities = get_entities_flat(wdl_wf)
        janis_wf = WorkflowBuilder(identifier=wdl_wf.name)
        self.ingest_workflow_inputs(wdl_wf, janis_wf, entities)
        self.ingest_workflow_steps(wdl_wf, janis_wf, entities)
        self.ingest_workflow_outputs(wdl_wf, janis_wf, entities)        
        return janis_wf

    def ingest_workflow_inputs(self, wdl_wf, janis_wf, entities) -> None:
        for inp in entities["inputs"]:
            self.ingest_workflow_input(wdl_wf, janis_wf, inp)

    def ingest_workflow_input(self, wdl_wf, janis_wf, wdl_inp) -> InputNode:
        parser = WorkflowInputParser(wdl_wf, janis_wf, wdl_inp) # type: ignore
        return parser.parse()
    
    def ingest_workflow_steps(self, wdl_wf, janis_wf, entities) -> None:
        for flatcall in entities['calls']:
            self.ingest_workflow_step(wdl_wf, janis_wf, flatcall)
        for flatcall in entities['calls']:
            self.ingest_workflow_step_inputs(wdl_wf, janis_wf, flatcall)
        for flatcall in entities['calls']:
            self.ingest_workflow_step_modifiers(wdl_wf, janis_wf, flatcall)

    def ingest_workflow_step(self, wdl_wf, janis_wf, flatcall) -> StepNode:
        call = flatcall.entity
        task = self.ingest(call.callee)
        return janis_wf.step(
            identifier=call.name,
            tool=task,
            ignore_missing=True
        )
            
    def ingest_workflow_step_inputs(self, wdl_wf, janis_wf, flatcall) -> StepNode:
        call = flatcall.entity
        task = call.callee
        inp_map = {}
        for k, v in call.inputs.items():
            if k in task.ignored_inputs:
                continue 
            parser = WorkflowStepInputParser(wdl_wf, janis_wf, v, flatcall)  
            inp_map[k] = parser.parse()
        inputs_dict = inp_map
        jstep = janis_wf[call.name]
        add_step_edges_to_graph(jstep, inputs_dict, janis_wf) 
        return jstep

    def ingest_workflow_step_modifiers(self, wdl_wf, janis_wf, flatcall) -> StepNode:
        parser = WorkflowStepModifierParser(wdl_wf, janis_wf, flatcall) 
        parser.parse()
        
        jstep = janis_wf[flatcall.entity.name]
        if parser.scatter is not None:
            jstep.scatter = parser.scatter
            janis_wf.has_scatter = True
        if parser.when is not None:
            jstep.when = parser.when
        
        return jstep

    def ingest_workflow_outputs(self, wdl_wf, janis_wf, entities) -> None:
        for flatout in entities["outputs"]:
            self.ingest_workflow_output(wdl_wf, janis_wf, flatout)

    def ingest_workflow_output(self, wdl_wf, janis_wf, flatout) -> OutputNode:
        parser = WorkflowOutputParser(wdl_wf, janis_wf, flatout) 
        return parser.parse()



    # @error_boundary()
    # def add_call_to_wf(
    #     self,
    #     wf: j.WorkflowBase,
    #     call: WDL.WorkflowNode,
    #     condition=None,
    #     foreach=None,
    #     expr_alias: str = None,
    # ):
       
    #     def selector_getter(exp):
    #         if exp == expr_alias:
    #             return j.ForEachSelector()

    #         return self.workflow_selector_getter(wf, exp)

    #     if isinstance(call, WDL.Call):
    #         task = self.ingest(call.callee)
    #         inp_map = {}
    #         for k, v in call.inputs.items():
    #             new_expr = self.translate_expr(v, input_selector_getter=selector_getter)
    #             inp_map[k] = new_expr
    #         return wf.step(call.name, task(**inp_map), when=condition, _foreach=foreach)

    #     elif isinstance(call, WDL.Conditional):
    #         # if len(call.body) > 1:
    #         #     raise NotImplementedError(
    #         #         f"Janis can't currently support more than one call inside the conditional: {', '.join(str(c) for c in call.body)}")
    #         for inner_call in call.body:
    #             # inner_call = call.body[0]
    #             self.add_call_to_wf(
    #                 wf,
    #                 inner_call,
    #                 condition=self.translate_expr(
    #                     call.expr, input_selector_getter=selector_getter
    #                 ),
    #                 expr_alias=expr_alias,
    #                 foreach=foreach,
    #             )
    #     elif isinstance(call, WDL.Scatter):
    #         # for scatter, we want to take the call.expr, and pass it to a step.foreach

    #         foreach = self.translate_expr(call.expr)

    #         scar_var_type = self.parse_wdl_type(call.expr.type)
    #         if isinstance(scar_var_type, WDL.Type.Array):
    #             scar_var_type = scar_var_type.item_type

    #         # when we unwrap each step-input to the workflow, we want to replace 'call.variable' with
    #         #       lambda el: <operation with call.variable substituted for {el}>
    #         # if call.variable not in wf.input_nodes:
    #         #     wf.input(call.variable, scar_var_type)
    #         for inner_call in call.body:
    #             self.add_call_to_wf(
    #                 wf, inner_call, foreach=foreach, expr_alias=call.variable
    #             )

    #     elif isinstance(call, WDL.Decl):
    #         self.add_decl_to_wf_input(wf, call)
    #     else:
    #         raise NotImplementedError(f"body type: {type(call)}")

    # def add_decl_to_wf_input(self, wf: j.WorkflowBase, inp: WDL.Decl):
    #     default = None
    #     if inp.expr:

    #         def selector_getter(exp):
    #             return self.workflow_selector_getter(wf, exp)

    #         default = self.translate_expr(
    #             inp.expr, input_selector_getter=selector_getter
    #         )

    #     return wf.input(inp.name, self.parse_wdl_type(inp.type, uuid=inp.name), default=default)






if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        raise Exception("Expected 1 argument, the name of a CWL tool.")

    toolname = sys.argv[1]

    try:
        tool = WdlParser.from_doc(toolname)
        tool.translate("janis")

    except WDL.Error.MultipleValidationErrors as err:
        for exc in err.exceptions:
            print(exc, file=sys.stderr)
            print(exc.pos, file=sys.stderr)
            print(exc.node, file=sys.stderr)
    except WDL.Error.ValidationError as exc:
        print(exc, file=sys.stderr)
        print(exc.pos, file=sys.stderr)
        print(exc.node, file=sys.stderr)








### DEPRECATED ####

    # ### OUTPUT ###
    # def parse_wf_output(self, wf: j.WorkflowBase, obj: WDL.Tree.Workflow, out: WDL.Tree.Decl):
    #     if out.expr is None:
    #         raise Exception(f"Output {out.name} has no expression")
    #     # TODO UPDATE EXPRESSION PARSING FOR WORKFLOW SCOPE
    #     dtype = parse_type(out.type, obj, uuid=wf.uuid)
        
    #     def selector_getter(exp):
    #         return self.workflow_selector_getter(wf, exp)
        
    #     sel = self.translate_expr(out.expr, input_selector_getter=selector_getter)
    #     wf.output(out.name, dtype, sel)


    # def translate_expr(
    #     self, expr: WDL.Expr.Base, input_selector_getter: Callable[[str], any] = None
    # ) -> Optional[Union[j.Selector, List[j.Selector], int, str, float, bool]]:
    #     if expr is None:
    #         return None

    #     tp = lambda exp: self.translate_expr(
    #         exp, input_selector_getter=input_selector_getter
    #     )

    #     if isinstance(expr, WDL.Expr.Array):
    #         # a literal array
    #         return [self.translate_expr(e) for e in expr.items]
    #     if isinstance(expr, WDL.Expr.String):
    #         return self.translate_wdl_string(expr)
    #     elif isinstance(expr, (WDL.Expr.Int, WDL.Expr.Boolean, WDL.Expr.Float)):
    #         return expr.literal.value
    #     if isinstance(expr, WDL.Expr.Placeholder):
    #         return self.translate_expr(expr.expr)
    #     if isinstance(expr, WDL.Expr.IfThenElse):
    #         return j.If(tp(expr.condition), tp(expr.consequent), tp(expr.alternative))
    #     elif isinstance(expr, WDL.Expr.Get):
    #         n = str(expr.expr)
    #         if input_selector_getter:
    #             return input_selector_getter(n)
    #         return j.InputSelector(n)
    #     elif isinstance(expr, WDL.Expr.Apply):
    #         return self.translate_apply(
    #             expr, input_selector_getter=input_selector_getter
    #         )

    #     raise Exception(f"Unsupported WDL expression type: {expr} ({type(expr)})")

    # ### EXPRESSIONS ###
    # def workflow_selector_getter(self, wf, exp: str):
    #     if "." in exp:
    #         node, *tag = exp.split(".")
    #         if len(tag) > 1:
    #             raise Exception(f"Couldn't parse source ID: {exp} - too many '.'")
    #         return wf[node][tag[0]]
    #     return wf[exp]
    
    # def translate_wdl_string(self, s: WDL.Expr.String):
    #     if s.literal is not None:
    #         return str(s.literal).lstrip('"').rstrip('"')

    #     elements = {}
    #     counter = 1
    #     _format = str(s).lstrip('"').rstrip('"')

    #     for placeholder in s.children:
    #         if isinstance(placeholder, (str, bool, int, float)):
    #             continue

    #         token = f"JANIS_WDL_TOKEN_{counter}"
    #         if str(placeholder) not in _format:
    #             # if the placeholder came up again
    #             continue

    #         _format = _format.replace(str(placeholder), f"{{{token}}}")
    #         elements[token] = self.translate_expr(placeholder)
    #         counter += 1

    #     if len(elements) == 0:
    #         return str(s)

    #     _format.replace("\\n", "\n")

    #     return j.StringFormatter(_format, **elements)

    # def file_size_operator(self, src, *args):
    #     multiplier = None
    #     if len(args) > 1:
    #         f = args[1].lower()
    #         multiplier_heirarchy = [
    #             ("ki" in f, 1024),
    #             ("k" in f, 1000),
    #             ("mi" in f, 1.024),
    #             ("gi" in f, 0.001024),
    #             ("g" in f, 0.001),
    #         ]
    #         if not any(m[0] for m in multiplier_heirarchy):
    #             j.Logger.warn(
    #                 f"Couldn't determine prefix {f} for FileSizeOperator, defaulting to MB"
    #             )
    #         else:
    #             multiplier = [m[1] for m in multiplier_heirarchy if m[0] is True][0]

    #     if isinstance(src, list):
    #         return multiplier * sum(j.FileSizeOperator(s) for s in src)

    #     base = j.FileSizeOperator(src, *args)
    #     if multiplier is not None and multiplier != 1:
    #         return multiplier * base
    #     return base

    # def basename_operator(self, src, *args):
    #     retval = j.BasenameOperator(src)
    #     if len(args) > 0:
    #         retval = retval.replace(args[0], "")

    #     return retval

    # def translate_apply(
    #     self, expr: WDL.Expr.Apply, **expr_kwargs
    # ) -> Union[j.Selector, List[j.Selector]]:

    #     # special case for select_first of array with one element
    #     if expr.function_name == "select_first" and len(expr.arguments) > 0:
    #         inner = expr.arguments[0]
    #         if isinstance(inner, WDL.Expr.Array) and len(inner.items) == 1:
    #             return self.translate_expr(inner.items[0]).assert_not_null()

    #     args = [self.translate_expr(e, **expr_kwargs) for e in expr.arguments]

    #     fn_map = {
    #         "_land": j.AndOperator,
    #         "defined": j.IsDefined,
    #         "select_first": j.FilterNullOperator,
    #         "basename": self.basename_operator,
    #         "length": j.LengthOperator,
    #         "_gt": j.GtOperator,
    #         "_gte": j.GteOperator,
    #         "_lt": j.LtOperator,
    #         "_lte": j.LteOperator,
    #         "sep": j.JoinOperator,
    #         "_add": j.AddOperator,
    #         "_interpolation_add": j.AddOperator,
    #         "stdout": j.Stdout,
    #         "_mul": j.MultiplyOperator,
    #         "_div": j.DivideOperator,
    #         "glob": j.WildcardSelector,
    #         "range": j.RangeOperator,
    #         "_at": j.IndexOperator,
    #         "_negate": j.NotOperator,
    #         "_sub": j.SubtractOperator,
    #         "size": self.file_size_operator,
    #         "ceil": j.CeilOperator,
    #         "select_all": j.FilterNullOperator,
    #         "sub": j.ReplaceOperator,
    #         "round": j.RoundOperator,
    #         "write_lines": lambda exp: f"JANIS: write_lines({exp})",
    #         "read_tsv": lambda exp: f"JANIS: j.read_tsv({exp})",
    #         "read_boolean": lambda exp: f"JANIS: j.read_boolean({exp})",
    #         "read_lines": lambda exp: f"JANIS: j.read_lines({exp})",
    #     }
    #     fn = fn_map.get(expr.function_name)
    #     if fn is None:
    #         raise Exception(f"Unhandled WDL apply function_name: {expr.function_name}")
    #     if isinstance(fn, LambdaType):
    #         return fn(args)
    #     return fn(*args)

    # ### DATATYPES ###
    # def parse_wdl_type(self, t: WDL.Type.Base, uuid: Optional[str]=None):
    #     optional = t.optional
    #     if isinstance(t, WDL.Type.Int):
    #         return j.Int(optional=optional)
    #     elif isinstance(t, WDL.Type.String):
    #         return j.String(optional=optional)
    #     elif isinstance(t, WDL.Type.Float):
    #         return j.Float(optional=optional)
    #     elif isinstance(t, WDL.Type.Boolean):
    #         return j.Boolean(optional=optional)
    #     elif isinstance(t, WDL.Type.File):
    #         return j.File(optional=optional)
    #     elif isinstance(t, WDL.Type.Directory):
    #         return j.Directory(optional=optional)
    #     elif isinstance(t, WDL.Type.Array):
    #         return j.Array(self.parse_wdl_type(t.item_type, uuid), optional=optional)
    #     elif isinstance(t, WDL.Type.StructInstance):
    #         if uuid:
    #             log_message(uuid, 'WDL Struct type unsupported. Has been cast to File type.', ErrorCategory.DATATYPES)
    #         return j.File(optional=optional)

    #     raise Exception(f"Didn't handle WDL type conversion for '{t}' ({type(t)})")





### DEPRECATED ###

    # cmdtool = j.CommandToolBuilder(
    #     tool=obj.name,
    #     version='DEV',
    #     container='ubuntu:latest',
    #     base_command=None,
    #     inputs=[],
    #     outputs=[]
    # )
    # # TODO something like this:
    # # parser = CommandParser()
    # # cmdtool._inputs, cmdtool.base_command, cmdtool._outputs = parser.parse()

    # cmdtool._base_command = self.parse_command_tool_base_command(obj) # type: ignore
    # cmdtool._container = self.container_from_runtime(obj.runtime, inputs=obj.inputs) # type: ignore
    
    # # inputs
    # for wdl_inp in obj.inputs:
    #     if wdl_inp.name.startswith("runtime_"):
    #         continue
    #     j_inp = self.parse_command_tool_input(wdl_inp)
    #     if j_inp is not None:
    #         cmdtool._inputs.append(j_inp)
    
    # # outputs
    # for wdl_out in obj.outputs:
    #     j_out = self.parse_command_tool_output(wdl_out)
    #     if j_out is not None:
    #         cmdtool._outputs.append(j_out)
    
    # # files to create
    # cmdtool._files_to_create = {"script.sh": self.translate_expr(obj.command)}
    # cmdtool._memory = self.parse_memory_requirement(obj.runtime.get("memory"))
    # cpus = self.translate_expr(obj.runtime.get("cpu"))
    # if isinstance(cpus, str):
    #     cpus = int(cpus)
    # cmdtool._cpus = cpus
    # cmdtool._disk = self.parse_disk_requirement(obj.runtime.get("disks"))

    # return cmdtool


            
    # def parse_command_tool_base_command(self, obj: WDL.Task) -> list[str]:
    #     return ["sh", "script.sh"]


    
    # def parse_command_tool_input(self, inp: WDL.Decl):
    #     default = None
    #     if inp.expr:
    #         default = self.translate_expr(inp.expr)

    #     # explicitly skip "runtime_*" inputs because they're from janis
    #     if inp.name.startswith("runtime_"):
    #         return None

    #     return j.ToolInput(inp.name, self.parse_wdl_type(inp.type, uuid=inp.name), default=default)

    # def parse_command_tool_output(self, outp: WDL.Decl):
    #     sel = self.translate_expr(outp.expr)

    #     return j.ToolOutput(outp.name, self.parse_wdl_type(outp.type, uuid=outp.name), selector=sel)