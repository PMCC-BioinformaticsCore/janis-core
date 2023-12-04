
import WDL 
from typing import Any, Optional, Tuple
from types import LambdaType
import regex as re 

import janis_core as j
from janis_core import settings
from janis_core.messages import load_loglines
from janis_core.messages import LogLine
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory
from .workflow.explore import CallContext
from .workflow.explore import OutputContext
from copy import deepcopy

### WDL EXPR -> STRING ###

def expr_as_str(expr: WDL.Expr.Base) -> str:
    with open(expr.pos.abspath, 'r') as fp:
        doclines = fp.readlines()
        doclines = [x.rstrip('\n') for x in doclines]
    lines = doclines[expr.pos.line - 1: expr.pos.end_line]
    lines[-1] = lines[-1][:expr.pos.end_column - 1]
    lines[0] = lines[0][expr.pos.column - 1:]
    text = '\n'.join(lines)
    return text



### WDL -> JANIS MAPPING ###

def parse_expr(
    expr: WDL.Expr.Base, 
    wdl_entity: WDL.Tree.Task | WDL.Tree.Workflow, 
    j_entity: j.CommandToolBuilder | j.WorkflowBuilder,
    node_context: Optional[CallContext | OutputContext]=None
    ) -> Tuple[Any, bool]:
    
    if settings.ingest.SAFE_MODE == True:
        try:
            parser = WDlExprParser(wdl_entity, j_entity, node_context)
            return parser.parse(expr), True
        
        except Exception as e:
            expr_str = str(expr)
            loglines = load_loglines(category=ErrorCategory.SCRIPTING, entity_uuids=set([j_entity.uuid]))
            
            # expr already has token?
            token = _get_token_for_expr(expr_str, loglines)
            if token:
                return token, False
            
            # expr needs new token
            token = f'__TOKEN{len(loglines) + 1}__'
            msg = f'{token} = "{expr}"'
            log_message(j_entity.uuid, msg, ErrorCategory.SCRIPTING)
            return token, False
    else:
        parser = WDlExprParser(wdl_entity, j_entity, node_context)
        return parser.parse(expr), True

def _get_token_for_expr(expr: str, loglines: list[LogLine]) -> Optional[str]:
    token_p = r'__TOKEN\d+__'
    for line in loglines:
        if re.match(token_p, line.message):
            token, expr = line.message.split(' = ', 1)
            if expr.strip('"') == expr:
                return token
    return None


class WDlExprParser:
    def __init__(
        self, 
        wdl_entity: WDL.Tree.Task | WDL.Tree.Workflow, 
        j_entity: j.CommandToolBuilder | j.WorkflowBuilder,
        node_context: Optional[CallContext | OutputContext]=None
        ) -> None:
        self.wdl_entity = wdl_entity     # wdl task or workflow
        self.j_entity = j_entity         # respective janis cmdtoolbuilder or workflowbuilder
        self.node_context = node_context # contextual information about scoped vars / scatter available in this scope

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

        # task
        if isinstance(self.wdl_entity, WDL.Tree.Task):
            return self.parse_get_task(expr)
        # workflow step
        elif isinstance(self.wdl_entity, WDL.Tree.Workflow):
            return self.parse_get_workflow(expr)
        else:
            raise RuntimeError
        
    def parse_get_task(self, expr: WDL.Expr.Get):
        # tool - input | temp var?
        assert isinstance(self.j_entity, j.CommandToolBuilder)
        return j.InputSelector(str(expr.expr.name))
    
    def parse_get_workflow(self, expr: WDL.Expr.Get):
        # workflow - input | step | step output | scatter target | temp var
        expr_str = str(expr.expr.name)

        if self.is_scatter_target(expr_str):
            return self.parse_get_scatter_target(expr_str)
        elif self.is_scoped_var(expr_str):
            return self.parse_get_scoped_var(expr_str)
        # condition (when)?
        elif "." in expr_str:
            return self.parse_get_stepout_ref(expr_str)
        else:
            return self.parse_get_input_ref(expr_str)

    def is_scatter_target(self, expr: str):
        assert isinstance(self.j_entity, j.WorkflowBuilder)
        if not isinstance(self.node_context, CallContext):
            return False 
        if self.node_context.scatter is None:
            return False
        if self.node_context.scatter.variable == expr:
            return True
        return False
    
    def is_scoped_var(self, expr: str):
        assert isinstance(self.j_entity, j.WorkflowBuilder)
        if self.node_context is None:
            return False
        for tempvar in self.node_context.scopedvars:
            if tempvar.name == expr:
                return True
        return False
    
    def parse_get_scatter_target(self, expr: str):
        # TODO maybe call parse_expr() rather than self.parse? want the largest traceback.
        assert isinstance(self.node_context, CallContext)
        assert self.node_context is not None
        assert self.node_context.scatter is not None
        return self.parse(self.node_context.scatter.expr)
    
    def parse_get_scoped_var(self, expr: str):
        assert self.node_context is not None
        tempvar = [x for x in self.node_context.scopedvars if x.name == expr][0]
        return self.parse(tempvar.expr)

    def parse_get_stepout_ref(self, expr: str):
        assert isinstance(self.j_entity, j.WorkflowBuilder)
        assert '.' in expr 
        stp_id, *tag = expr.split(".")
        if len(tag) > 1:
            raise Exception(f"Couldn't parse source ID: {expr} - too many '.'")
        assert stp_id in self.j_entity.step_nodes
        stp = self.j_entity.step_nodes[stp_id]
        sout = tag[0]
        selector = stp.get_item(sout)
        return selector
        
    def parse_get_input_ref(self, expr: str):
        assert isinstance(self.j_entity, j.WorkflowBuilder)
        if expr in self.j_entity.input_nodes:
            node = self.j_entity.input_nodes[expr]
            return j.InputNodeSelector(node)
        raise NotImplementedError
    
    def parse_string(self, s: WDL.Expr.String):
        if s.literal is not None:
            return str(s.literal).lstrip('"').rstrip('"')

        elements = {}
        counter = 1
        _format = str(s).lstrip('"').rstrip('"')

        for placeholder in s.children:
            if isinstance(placeholder, (str, bool, int, float)):
                continue

            token = f"TOKEN{counter}"
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
    
    def parse_prefix(self, src, *args):
        raise NotImplementedError
    
    def parse_read_lines(self, src, *args):
        if len(args) > 0:
            raise NotImplementedError
        return j.SplitOperator(j.ReadContents(src), '\n')
    
    def parse_read_int(self, src, *args):
        if len(args) > 0:
            raise NotImplementedError
        return j.AsIntOperator(j.ReadContents(src))
    
    def parse_read_float(self, src, *args):
        if len(args) > 0:
            raise NotImplementedError
        return j.AsFloatOperator(j.ReadContents(src))
    
    def parse_read_boolean(self, src, *args):
        if len(args) > 0:
            raise NotImplementedError
        return j.AsBoolOperator(j.ReadContents(src))
    
    def parse_apply(self, expr: WDL.Expr.Apply) -> j.Selector | list[j.Selector]:
        args = [self.parse(e) for e in expr.arguments]

        fn_map = {
            'read_lines': self.parse_read_lines,
            'read_tsv': None,
            'read_json': j.ReadJsonOperator,
            'read_map': None,
            'read_object': None,
            'read_objects': None,
            'read_string': j.ReadContents,
            'read_int': self.parse_read_int,
            'read_float': self.parse_read_float,
            'read_boolean': self.parse_read_boolean,
            'write_lines': None,
            'write_tsv': None,
            'write_json': None,
            'write_map': None,
            'write_object': None,
            'write_objects': None,
            "range": j.RangeOperator,
            'transpose': j.TransposeOperator,
            'zip': None,
            'cross': None,
            "length": j.LengthOperator,
            'flatten': j.FlattenOperator,
            'prefix': self.parse_prefix,
            "select_first": j.FirstOperator,
            'select_all': j.FilterNullOperator,
            'defined': j.IsDefined,
            "basename": self.parse_basename,
            'floor': j.FloorOperator,
            "sep": j.JoinOperator,
            "stdout": j.Stdout,
            "glob": j.WildcardSelector,
            "size": self.parse_file_size,
            "ceil": j.CeilOperator,
            "sub": j.ReplaceOperator,
            "round": j.RoundOperator,
            
            "_lor": j.OrOperator,
            "_eqeq": j.EqualityOperator,
            "_land": j.AndOperator,
            "_gt": j.GtOperator,
            "_gte": j.GteOperator,
            "_lt": j.LtOperator,
            "_lte": j.LteOperator,
            "_add": j.AddOperator,
            "_interpolation_add": j.AddOperator,
            "_add": j.AddOperator,
            "_mul": j.MultiplyOperator,
            "_div": j.DivideOperator,
            "_at": j.IndexOperator,
            "_negate": j.NotOperator,
            "_sub": j.SubtractOperator,
        }
        # TODO 
        # CWL flatten: https://github.com/common-workflow-library/cwl-patterns/blob/main/javascript_snippets/flatten-nestedarray.cwl

        # special case for select_first of array with one element
        if expr.function_name == "select_first" and len(expr.arguments) > 0:
            inner = expr.arguments[0]
            if isinstance(inner, WDL.Expr.Array) and len(inner.items) == 1:
                return self.parse(inner.items[0]).assert_not_null()
        
        # uncaught error.  
        # log message and bail
        if expr.function_name not in fn_map:
            raise Exception(f"Unhandled WDL apply function_name: {expr.function_name}")
        
        fn = fn_map[expr.function_name]
        
        # caught error.  we can't parse this func. ignore, log a message, keep going.
        if fn is None:
            arg = self.parse(expr.arguments[0])
            msg = f"Function {expr.function_name}({arg}) not supported. Ignored function."
            log_message(self.j_entity.uuid, msg, ErrorCategory.SCRIPTING)
            return arg

        if isinstance(fn, LambdaType):
            return fn(args)
        return fn(*args)

