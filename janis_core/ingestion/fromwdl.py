#!/usr/bin/env python3
import functools
import os
import re
from types import LambdaType

from typing import List, Union, Optional, Callable
import WDL

import janis_core as j


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

    @error_boundary()
    def add_call_to_wf(
        self,
        wf: j.WorkflowBase,
        call: WDL.WorkflowNode,
        condition=None,
        foreach=None,
        expr_alias: str = None,
    ):
        def selector_getter(exp):
            if exp == expr_alias:
                return j.ForEachSelector()

            return self.workflow_selector_getter(wf, exp)

        if isinstance(call, WDL.Call):
            task = self.from_loaded_object(call.callee)
            inp_map = {}
            for k, v in call.inputs.items():
                new_expr = self.translate_expr(v, input_selector_getter=selector_getter)

                inp_map[k] = new_expr

            return wf.step(call.name, task(**inp_map), when=condition, _foreach=foreach)

        elif isinstance(call, WDL.Conditional):
            # if len(call.body) > 1:
            #     raise NotImplementedError(
            #         f"Janis can't currently support more than one call inside the conditional: {', '.join(str(c) for c in call.body)}")
            for inner_call in call.body:
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
        elif isinstance(call, WDL.Scatter):
            # for scatter, we want to take the call.expr, and pass it to a step.foreach

            foreach = self.translate_expr(call.expr)

            scar_var_type = self.parse_wdl_type(call.expr.type)
            if isinstance(scar_var_type, WDL.Type.Array):
                scar_var_type = scar_var_type.item_type

            # when we unwrap each step-input to the workflow, we want to replace 'call.variable' with
            #       lambda el: <operation with call.variable substituted for {el}>
            # if call.variable not in wf.input_nodes:
            #     wf.input(call.variable, scar_var_type)
            for inner_call in call.body:
                self.add_call_to_wf(
                    wf, inner_call, foreach=foreach, expr_alias=call.variable
                )

        elif isinstance(call, WDL.Decl):
            self.add_decl_to_wf_input(wf, call)
        else:
            raise NotImplementedError(f"body type: {type(call)}")

    def add_decl_to_wf_input(self, wf: j.WorkflowBase, inp: WDL.Decl):
        default = None
        if inp.expr:

            def selector_getter(exp):
                return self.workflow_selector_getter(wf, exp)

            default = self.translate_expr(
                inp.expr, input_selector_getter=selector_getter
            )

        return wf.input(inp.name, self.parse_wdl_type(inp.type), default=default)

    @classmethod
    def container_from_runtime(cls, runtime, inputs: List[WDL.Decl]):
        container = runtime.get("container", runtime.get("docker"))
        if isinstance(container, WDL.Expr.Get):
            # relevant input
            inp = [i.expr for i in inputs if i.name == str(container.expr)]
            if len(inp) > 0:
                container = inp[0]
            else:
                j.Logger.warn(
                    f"Expression for determining containers was '{container}' "
                    f"but couldn't find input called {str(container.expr)}"
                )
        if isinstance(container, WDL.Expr.String):
            container = container.literal
        if isinstance(container, WDL.Value.String):
            container = container.value
        if container is None:
            container = "ubuntu:latest"
        if not isinstance(container, str):
            j.Logger.warn(
                f"Expression for determining containers ({container}) are not supported in Janis, using ubuntu:latest"
            )
            container = "ubuntu:latest"
        return container

    def parse_memory_requirement(self, value):
        s = self.translate_expr(value)
        if s is None:
            return 1.074
        elif isinstance(s, str):
            if s.lower().endswith("g"):
                return float(s[:-1].strip())
            if s.lower().endswith("gb"):
                return float(s[:-2].strip())
            elif s.lower().endswith("gib"):
                return float(s[:-3].strip()) * 1.074
            elif s.lower().endswith("mb"):
                return float(s[:-2].strip()) / 1000
            elif s.lower().endswith("mib"):
                return float(s[:-3].strip()) / 1024
            raise Exception(f"Memory type {s}")
        elif isinstance(s, (float, int)):
            # in bytes?
            return s / (1024 ** 3)
        elif isinstance(s, j.Selector):
            return s
        raise Exception(f"Couldn't recognise memory requirement '{value}'")

    def parse_disk_requirement(self, value):
        s = self.translate_expr(value)
        if isinstance(s, str):
            try:
                return int(s)
            except ValueError:
                pass
            pattern_matcher = re.match(r"local-disk (\d+) .*", s)
            if not pattern_matcher:
                raise Exception(f"Couldn't recognise disk type '{value}'")
            s = pattern_matcher.groups()[0]
            try:
                return int(s)
            except ValueError:
                pass
            if s.lower().endswith("gb"):
                return float(s[:-2].strip())
            elif s.lower().endswith("gib"):
                return float(s[:-3].strip()) * 1.074
            elif s.lower().endswith("mb"):
                return float(s[:-2].strip()) / 1000
            elif s.lower().endswith("mib"):
                return float(s[:-3].strip()) / 1024
            raise Exception(f"Disk type type {s}")
        elif isinstance(s, (float, int)):
            # in GiB
            return s * 1.07374
        elif isinstance(s, j.Selector):
            return s
        elif s is None:
            return 2.14748  # 2 GiB
        raise Exception(f"Couldn't recognise memory requirement '{value}'")

    def from_loaded_task(self, obj: WDL.Task):
        rt = obj.runtime
        translated_script = self.translate_expr(obj.command)
        inputs = obj.inputs

        cpus = self.translate_expr(rt.get("cpu"))
        if not isinstance(cpus, j.Selector) and cpus is not None and not isinstance(cpus, (int, float)):
            cpus = int(cpus)

        c = j.CommandToolBuilder(
            tool=obj.name,
            base_command=["sh", "script.sh"],
            container=self.container_from_runtime(rt, inputs=inputs),
            version="DEV",
            inputs=[
                self.parse_command_tool_input(i)
                for i in obj.inputs
                if not i.name.startswith("runtime_")
            ],
            outputs=[self.parse_command_tool_output(o) for o in obj.outputs],
            files_to_create={"script.sh": translated_script},
            memory=self.parse_memory_requirement(rt.get("memory")),
            cpus=cpus,
            disk=self.parse_disk_requirement(rt.get("disks")),
        )

        return c

    def translate_expr(
        self, expr: WDL.Expr.Base, input_selector_getter: Callable[[str], any] = None
    ) -> Optional[Union[j.Selector, List[j.Selector], int, str, float, bool]]:
        if expr is None:
            return None

        tp = lambda exp: self.translate_expr(
            exp, input_selector_getter=input_selector_getter
        )

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
            return self.translate_apply(
                expr, input_selector_getter=input_selector_getter
            )

        raise Exception(f"Unsupported WDL expression type: {expr} ({type(expr)})")

    def translate_wdl_string(self, s: WDL.Expr.String):
        if s.literal is not None:
            return str(s.literal).lstrip('"').rstrip('"')

        elements = {}
        counter = 1
        _format = str(s).lstrip('"').rstrip('"')

        for placeholder in s.children:
            if isinstance(placeholder, (str, bool, int, float)):
                continue

            token = f"JANIS_WDL_TOKEN_{counter}"
            if str(placeholder) not in _format:
                # if the placeholder came up again
                continue

            _format = _format.replace(str(placeholder), f"{{{token}}}")
            elements[token] = self.translate_expr(placeholder)
            counter += 1

        if len(elements) == 0:
            return str(s)

        _format.replace("\\n", "\n")

        return j.StringFormatter(_format, **elements)

    def file_size_operator(self, src, *args):
        multiplier = None
        if len(args) > 1:
            f = args[1].lower()
            multiplier_heirarchy = [
                ("ki" in f, 1024),
                ("k" in f, 1000),
                ("mi" in f, 1.024),
                ("gi" in f, 0.001024),
                ("g" in f, 0.001),
            ]
            if not any(m[0] for m in multiplier_heirarchy):
                j.Logger.warn(
                    f"Couldn't determine prefix {f} for FileSizeOperator, defaulting to MB"
                )
            else:
                multiplier = [m[1] for m in multiplier_heirarchy if m[0] is True][0]

        if isinstance(src, list):
            return multiplier * sum(j.FileSizeOperator(s) for s in src)

        base = j.FileSizeOperator(src, *args)
        if multiplier is not None and multiplier != 1:
            return multiplier * base
        return base

    def basename_operator(self, src, *args):
        retval = j.BasenameOperator(src)
        if len(args) > 0:
            retval = retval.replace(args[0], "")

        return retval

    def translate_apply(
        self, expr: WDL.Expr.Apply, **expr_kwargs
    ) -> Union[j.Selector, List[j.Selector]]:

        # special case for select_first of array with one element
        if expr.function_name == "select_first" and len(expr.arguments) > 0:
            inner = expr.arguments[0]
            if isinstance(inner, WDL.Expr.Array) and len(inner.items) == 1:
                return self.translate_expr(inner.items[0]).assert_not_null()

        args = [self.translate_expr(e, **expr_kwargs) for e in expr.arguments]

        fn_map = {
            "_land": j.AndOperator,
            "defined": j.IsDefined,
            "select_first": j.FilterNullOperator,
            "basename": self.basename_operator,
            "length": j.LengthOperator,
            "_gt": j.GtOperator,
            "_gte": j.GteOperator,
            "_lt": j.LtOperator,
            "_lte": j.LteOperator,
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
            "_sub": j.SubtractOperator,
            "size": self.file_size_operator,
            "ceil": j.CeilOperator,
            "select_all": j.FilterNullOperator,
            "sub": j.ReplaceOperator,
            "round": j.RoundOperator,
            "write_lines": lambda exp: f"JANIS: write_lines({exp})",
            "read_tsv": lambda exp: f"JANIS: j.read_tsv({exp})",
            "read_boolean": lambda exp: f"JANIS: j.read_boolean({exp})",
            "read_lines": lambda exp: f"JANIS: j.read_lines({exp})",
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
