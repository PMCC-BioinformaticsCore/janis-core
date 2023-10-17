
import WDL 
from typing import Optional, Any, Callable
from types import LambdaType
import janis_core as j

def parse_expr(
    expr: WDL.Expr.Base, 
    input_selector_getter: Optional[Callable[[str], Any]]=None
    ) -> Any:
    
    if expr is None:
        return None

    tp = lambda exp: parse_expr(
        exp, input_selector_getter=input_selector_getter
    )

    if isinstance(expr, WDL.Expr.Array):
        # a literal array
        return [parse_expr(e) for e in expr.items]
    if isinstance(expr, WDL.Expr.String):
        return parse_string(expr)
    elif isinstance(expr, (WDL.Expr.Int, WDL.Expr.Boolean, WDL.Expr.Float)):
        return expr.literal.value
    if isinstance(expr, WDL.Expr.Placeholder):
        return parse_expr(expr.expr)
    if isinstance(expr, WDL.Expr.IfThenElse):
        return j.If(tp(expr.condition), tp(expr.consequent), tp(expr.alternative))
    elif isinstance(expr, WDL.Expr.Get):
        n = str(expr.expr)
        if input_selector_getter:
            return input_selector_getter(n)
        return j.InputSelector(n)
    elif isinstance(expr, WDL.Expr.Apply):
        return parse_apply(
            expr, input_selector_getter=input_selector_getter
        )

    raise Exception(f"Unsupported WDL expression type: {expr} ({type(expr)})")

def parse_string(s: WDL.Expr.String):
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
        elements[token] = parse_expr(placeholder)
        counter += 1

    if len(elements) == 0:
        return str(s)

    _format.replace("\\n", "\n")

    return j.StringFormatter(_format, **elements)


def file_size_operator(src, *args):
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

def basename_operator(src, *args):
    retval = j.BasenameOperator(src)
    if len(args) > 0:
        retval = retval.replace(args[0], "")
    return retval

def parse_apply(expr: WDL.Expr.Apply, **expr_kwargs) -> j.Selector | list[j.Selector]:

    # special case for select_first of array with one element
    if expr.function_name == "select_first" and len(expr.arguments) > 0:
        inner = expr.arguments[0]
        if isinstance(inner, WDL.Expr.Array) and len(inner.items) == 1:
            return parse_expr(inner.items[0]).assert_not_null()

    args = [parse_expr(e, **expr_kwargs) for e in expr.arguments]

    fn_map = {
        "_land": j.AndOperator,
        "defined": j.IsDefined,
        "select_first": j.FilterNullOperator,
        "basename": basename_operator,
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
        "size": file_size_operator,
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