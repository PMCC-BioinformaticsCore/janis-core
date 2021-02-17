import os
from types import LambdaType
from typing import List, Union
import WDL

import janis_core as j


class WdlParser:
    @staticmethod
    def from_doc(doc: str, base_uri=None):
        abs_path = os.path.relpath(doc)
        d = WDL.load(abs_path)

        parser = WdlParser()

        tasks = []
        for t in d.tasks:
            tasks.append(parser.from_loaded_object(t))

        return tasks[0]

    def from_loaded_object(self, obj: WDL.Task):
        rt = obj.runtime
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
            outputs=[],
            files_to_create={"script.sh": self.translate_expr(obj.command)},
        )

        return c

    def translate_expr(
        self, expr: WDL.Expr.Base
    ) -> Union[j.Selector, List[j.Selector], int, str, float, bool]:

        tp = self.translate_expr

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
            return j.InputSelector(str(expr.expr))
        elif isinstance(expr, WDL.Expr.Apply):
            return self.translate_apply(expr)

        raise Exception(f"Unsupported WDL expression type: {expr} ({type(expr)})")

    def translate_wdl_string(self, s: WDL.Expr.String):
        if not s.command:
            return str(s.literal)

        elements = {}
        counter = 1
        _format = str(s)

        for placeholder in s.children:
            if isinstance(placeholder, (str, bool, int, float)):
                continue

            token = f"JANIS_WDL_TOKEN_{counter}"
            _format.replace(str(placeholder), token)
            elements[token] = self.translate_expr(placeholder)

        if len(elements) > 0:
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
            "_interpolation_add": j.AddOperator,
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
        if inp.expr:
            print(
                f"Input {inp.name} has expression, need to work out how to parse this: '{inp.expr}'"
            )

        # explicitly skip "runtime_*" inputs because they're from janis
        if inp.name.startswith("runtime_"):
            return None

        return j.ToolInput(inp.name, self.parse_wdl_type(inp.type))


if __name__ == "__main__":
    doc = "path/to/doc.wdl"
    t = WdlParser.from_doc(doc)

    t.translate("wdl")
