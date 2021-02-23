import os
from types import LambdaType
from typing import List, Union, Optional, Callable
import WDL

import janis_core as j


class WdlParser:

    @staticmethod
    def from_doc(doc: str, base_uri=None):
        abs_path = os.path.relpath(doc)
        d = WDL.load(abs_path)

        parser = WdlParser()

        if d.workflow:
            return parser.from_loaded_object(d.workflow)

        tasks = []
        for t in d.tasks:
            tasks.append(parser.from_loaded_object(t))

        return tasks[0]

    def from_loaded_object(self, obj: WDL.SourceNode):
        if isinstance(obj, WDL.Task):
            return self.from_loaded_task(obj)
        elif isinstance(obj, WDL.Workflow):
            return self.from_loaded_workflow(obj)


    def from_loaded_workflow(self, obj: WDL.Workflow):
        wf = j.WorkflowBuilder(identifier=obj.name)

        for inp in obj.inputs:
            self.add_decl_to_wf_input(wf, inp)

        for call in obj.body:
            self.add_call_to_wf(wf, call)

        return wf

    def workflow_selector_getter(self, wf, exp: str):
        if "." in exp:
            node, *tag = exp.split(".")
            if len(tag) > 1:
                raise Exception(f"Couldn't parse source ID: {exp} - too many '.'")
            return wf[node][tag[0]]

        return wf[exp]

    def add_call_to_wf(self, wf: j.WorkflowBase, call: WDL.WorkflowNode, condition=None, scatter=None):
        selector_getter = lambda exp: self.workflow_selector_getter(wf, exp)

        if isinstance(call, WDL.Call):
            task = self.from_loaded_object(call.callee)
            inp_map = {k: self.translate_expr(v, input_selector_getter=selector_getter) for k, v in call.inputs.items()}
            return wf.step(call.name, task(**inp_map), when=condition, scatter=scatter)

        elif isinstance(call, WDL.Conditional):
            # if len(call.body) > 1:
            #     raise NotImplementedError(
            #         f"Janis can't currently support more than one call inside the conditional: {', '.join(str(c) for c in call.body)}")
            for inner_call in call.body:
                inner_call = call.body[0]
                self.add_call_to_wf(wf, inner_call, condition=self.translate_expr(call.expr, input_selector_getter=selector_getter))
        elif isinstance(call, WDL.Scatter):
            scatter = None # self.translate_expr(call.expr)
            scar_var_type = self.parse_wdl_type(call.expr.type)
            if isinstance(scar_var_type, WDL.Type.Array):
                scar_var_type = scar_var_type.item_type

            if call.variable not in wf.input_nodes:
                wf.input(call.variable, scar_var_type)
            for inner_call in call.body:
                self.add_call_to_wf(wf, inner_call, scatter=scatter)


        elif isinstance(call, WDL.Decl):
            self.add_decl_to_wf_input(wf, call)
        else:
            raise NotImplementedError(f"body type: {type(call)}")

    def add_decl_to_wf_input(self, wf: j.WorkflowBase, inp: WDL.Decl):
        default = None
        if inp.expr:
            default = self.translate_expr(inp.expr)

        return wf.input(inp.name, self.parse_wdl_type(inp.type), default=default)




    def from_loaded_task(self, obj: WDL.Task):
        rt = obj.runtime
        translated_script = self.translate_expr(obj.command)
        c = j.CommandToolBuilder(
            tool=obj.name,
            base_command=["sh", "script.sh"],
            container=rt.get("container", "ubuntu:latest"),
            version="DEV",
            inputs=[
                self.parse_command_tool_input(i)
                for i in obj.inputs
                if not i.name.startswith("runtime_")
            ],
            outputs=[self.parse_command_tool_output(o) for o in obj.outputs],
            files_to_create={"script.sh": translated_script},
        )

        return c

    def translate_expr(
        self, expr: WDL.Expr.Base, input_selector_getter: Callable[[str], any]=None
    ) -> Optional[Union[j.Selector, List[j.Selector], int, str, float, bool]]:
        if expr is None:
            return None

        tp = lambda exp: self.translate_expr(exp, input_selector_getter=input_selector_getter)

        if isinstance(expr, WDL.Expr.Array):
            # a literal array
            return [self.translate_expr(e) for e in expr.items]
        if isinstance(expr, WDL.Expr.String):
            return self.translate_wdl_string(expr)
        elif isinstance(expr, (WDL.Expr.Int, WDL.Expr.Boolean, WDL.Expr.Float)):
            return expr.literal.value
        if isinstance(expr, WDL.Expr.Placeholder):
            return self.translate_expr(expr.expr)
        if isinstance(expr, WDL.Expr.IfThenElse):
            return j.If(tp(expr.condition), tp(expr.consequent), tp(expr.alternative))
        elif isinstance(expr, WDL.Expr.Get):
            n = str(expr.expr)
            if input_selector_getter:
                return input_selector_getter(n)
            return j.InputSelector(n)
        elif isinstance(expr, WDL.Expr.Apply):
            return self.translate_apply(expr)

        raise Exception(f"Unsupported WDL expression type: {expr} ({type(expr)})")

    def translate_wdl_string(self, s: WDL.Expr.String):
        if not s.command:
            return str(s.literal).lstrip('"').rstrip('"')

        elements = {}
        counter = 1
        _format = str(s)

        for placeholder in s.children:
            if isinstance(placeholder, (str, bool, int, float)):
                continue

            token = f"JANIS_WDL_TOKEN_{counter}"
            _format = _format.replace(str(placeholder), f"{{{token}}}")
            elements[token] = self.translate_expr(placeholder)

        if len(elements) == 0:
            return str(s)

        _format.replace("\\n", "\n")

        return j.StringFormatter(_format, **elements)

    def translate_apply(
        self, expr: WDL.Expr.Apply
    ) -> Union[j.Selector, List[j.Selector]]:

        # special case for select_first of array with one element
        if expr.function_name == "select_first" and len(expr.arguments) > 0:
            inner = expr.arguments[0]
            if isinstance(inner, WDL.Expr.Array) and len(inner.items) == 1:
                return self.translate_expr(inner.items[0]).assert_not_null()

        args = [self.translate_expr(e) for e in expr.arguments]

        fn_map = {
            "_land": j.AndOperator,
            "defined": j.IsDefined,
            "select_first": j.FilterNullOperator,
            "basename": j.BasenameOperator,
            "length": j.LengthOperator,
            "_gt": j.GtOperator,
            "_gte": j.GteOperator,
            "sep": j.JoinOperator,
            "_add": j.AddOperator,
            "_interpolation_add": j.AddOperator,
            "stdout": j.Stdout,
            "_mul": j.MultiplyOperator,
            "_div": j.DivideOperator,
            "glob": j.WildcardSelector,
            "range": j.RangeOperator,
            "_at": j.IndexOperator,
            "_negate": j.NotOperator,
            "write_lines": lambda exp: f"write_lines({exp})"
        }
        fn = fn_map.get(expr.function_name)
        if fn is None:
            raise Exception(f"Unhandled WDL apply function_name: {expr.function_name}")
        if isinstance(fn, LambdaType):
            return fn(args)
        return fn(*args)

    def parse_wdl_type(self, t: WDL.Type.Base):
        optional = t.optional
        if isinstance(t, WDL.Type.Int):
            return j.Int(optional=optional)
        elif isinstance(t, WDL.Type.String):
            return j.String(optional=optional)
        elif isinstance(t, WDL.Type.Float):
            return j.Float(optional=optional)
        elif isinstance(t, WDL.Type.Boolean):
            return j.Boolean(optional=optional)
        elif isinstance(t, WDL.Type.File):
            return j.File(optional=optional)
        elif isinstance(t, WDL.Type.Directory):
            return j.Directory(optional=optional)
        elif isinstance(t, WDL.Type.Array):
            return j.Array(self.parse_wdl_type(t.item_type), optional=optional)

        raise Exception(f"Didn't handle WDL type conversion for '{t}' ({type(t)})")

    def parse_command_tool_input(self, inp: WDL.Decl):
        default = None
        if inp.expr:
            default = self.translate_expr(inp.expr)

        # explicitly skip "runtime_*" inputs because they're from janis
        if inp.name.startswith("runtime_"):
            return None

        return j.ToolInput(inp.name, self.parse_wdl_type(inp.type), default=default)

    def parse_command_tool_output(self, outp: WDL.Decl):
        sel = self.translate_expr(outp.expr)

        return j.ToolOutput(outp.name, self.parse_wdl_type(outp.type), selector=sel)


if __name__ == "__main__":
    doc = "path/to/doc.wdl"
    t = WdlParser.from_doc(doc)

    t.translate("janis")
