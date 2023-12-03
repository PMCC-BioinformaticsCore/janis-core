
import WDL 
from typing import Any
from types import LambdaType
import janis_core as j


def parse_expr(expr: WDL.Expr.Base, wdl_entity: WDL.Tree.Task | WDL.Tree.Workflow, j_entity: j.CommandToolBuilder | j.WorkflowBuilder) -> Any:
    parser = WDlExprParser(wdl_entity, j_entity)
    return parser.parse(expr)

class WDlExprParser:
    def __init__(self, wdl_entity: WDL.Tree.Task | WDL.Tree.Workflow, j_entity: j.CommandToolBuilder | j.WorkflowBuilder):
        self.wdl_entity = wdl_entity    # wdl task or workflow
        self.j_entity = j_entity        # respective janis cmdtoolbuilder or workflowbuilder

    def parse(self, expr: WDL.Expr.Base) -> Any:
        if expr is None:
            return None
        if isinstance(expr, WDL.Expr.Array):
            # a literal array
            return [self.parse(e) for e in expr.items]
        if isinstance(expr, WDL.Expr.String):
            return self.parse_string(expr)
        elif isinstance(expr, (WDL.Expr.Int, WDL.Expr.Boolean, WDL.Expr.Float)):
            return expr.literal.value
        if isinstance(expr, WDL.Expr.Placeholder):
            return self.parse(expr.expr)
        if isinstance(expr, WDL.Expr.IfThenElse):
            return j.If(self.parse(expr.condition), self.parse(expr.consequent), self.parse(expr.alternative))
        elif isinstance(expr, WDL.Expr.Get):
            return self.parse_get(expr)
        elif isinstance(expr, WDL.Expr.Apply):
            return self.parse_apply(expr)
        raise Exception(f"Unsupported WDL expression type: {expr} ({type(expr)})")

    def parse_get(self, expr: WDL.Expr.Get):
        if isinstance(expr.expr, WDL.Expr.Get):
            expr = expr.expr
        assert isinstance(expr.expr, WDL.Expr.Ident)

        # tool - input | temp var?
        if isinstance(self.wdl_entity, WDL.Tree.Task):
            assert isinstance(self.j_entity, j.CommandToolBuilder)
            return j.InputSelector(str(expr.expr.name))
        
        # workflow - input | step | step output | temp var?
        elif isinstance(self.wdl_entity, WDL.Tree.Workflow):
            assert isinstance(self.j_entity, j.WorkflowBuilder)
            expr_str = str(expr.expr.name)
            if "." in expr_str:
                node, *tag = expr_str.split(".")
                if len(tag) > 1:
                    raise Exception(f"Couldn't parse source ID: {expr_str} - too many '.'")
                return self.j_entity[node][tag[0]]
            return self.j_entity[expr_str]

        else:
            raise RuntimeError

    def parse_string(self, s: WDL.Expr.String):
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
            elements[token] = self.parse(placeholder)
            counter += 1

        if len(elements) == 0:
            return str(s)

        _format.replace("\\n", "\n")

        return j.StringFormatter(_format, **elements)

    def parse_file_size(self, src, *args):
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

    def parse_basename(self, src, *args):
        retval = j.BasenameOperator(src)
        if len(args) > 0:
            retval = retval.replace(args[0], "")
        return retval

    def parse_apply(self, expr: WDL.Expr.Apply) -> j.Selector | list[j.Selector]:

        # special case for select_first of array with one element
        if expr.function_name == "select_first" and len(expr.arguments) > 0:
            inner = expr.arguments[0]
            if isinstance(inner, WDL.Expr.Array) and len(inner.items) == 1:
                return self.parse(inner.items[0]).assert_not_null()

        args = [self.parse(e) for e in expr.arguments]

        fn_map = {
            "_land": j.AndOperator,
            "defined": j.IsDefined,
            "select_first": j.FilterNullOperator,
            "basename": self.parse_basename,
            "length": j.LengthOperator,
            "_gt": j.GtOperator,
            "_gte": j.GteOperator,
            "_lt": j.LtOperator,
            "_lte": j.LteOperator,
            "sep": j.JoinOperator,
            "_add": j.AddOperator,
            "_interpolation_add": j.AddOperator,
            "stdout": j.Stdout,
            "_add": j.AddOperator,
            "_mul": j.MultiplyOperator,
            "_div": j.DivideOperator,
            "glob": j.WildcardSelector,
            "range": j.RangeOperator,
            "_at": j.IndexOperator,
            "_negate": j.NotOperator,
            "_sub": j.SubtractOperator,
            "size": self.parse_file_size,
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

