

import os
from lark import Lark
from lark import Tree
from lark import Token
from typing import Tuple, Any
import janis_core as j
from janis_core.messages import log_error
from janis_core.messages import get_messages
from janis_core.messages import ErrorCategory

GRAMMAR_PATH = f'{os.path.dirname(os.path.abspath(__file__))}/grammar.ebnf'

def parse_expression(expr: str, entity: j.CommandToolBuilder | j.WorkflowBuilder) -> Any:
    result, success = ExpressionParser().parse(expr)
    if not success:
        messages = get_messages(entity.id(), level='ERROR')
        messages = [x for x in messages if x.category == ErrorCategory.SCRIPT]
        new_token = f'__TOKEN{len(messages) + 1}__'
        msg = f'{new_token} = "{expr}"'
        log_error(entity.id(), ErrorCategory.SCRIPT, msg)
        return new_token, False
    return result, success

class ExpressionParser:
    file_attr_map = {
        'attr_basename': j.BasenameOperator,
        'attr_dirname': j.DirnameOperator,
        'attr_nameroot': j.NamerootOperator,
        'attr_nameext': j.NameextOperator,
        'attr_size': j.FileSizeOperator,
        'attr_contents': j.ReadContents,
        'attr_length': j.LengthOperator,
    }

    math_map = {
        'floor': j.FloorOperator,
        'ceil': j.CeilOperator,
        'round': j.RoundOperator,
    }
    
    two_value_map = {
        'and': j.AndOperator,
        'or': j.OrOperator,
        'deep_eq': j.EqualityOperator,
        'eq': j.EqualityOperator,
        'deep_ineq': j.InequalityOperator,
        'ineq': j.InequalityOperator,
        'gteq': j.GteOperator,
        'gt': j.GtOperator,
        'lteq': j.LteOperator,
        'lt': j.LtOperator,
        'add': j.AddOperator,
        'sub': j.SubtractOperator,
        'mul': j.MultiplyOperator,
        'div': j.DivideOperator,
    }
    
    with open(GRAMMAR_PATH) as fp:
        grammar = fp.read()
    parser = Lark(grammar, start='the_text')
    
    def parse(self, expr: str) -> Tuple[Any, bool]:
        try:
            tree = self.parser.parse(expr)
            return self.parse_node(tree), True
        except Exception as e:
            return None, False

    def generate_string_formatter(self, node: Tree) -> j.StringFormatter:
        token_replacers = {}
        string_format = ""
        for child in node.children:
            if child.data == 'javascript':
                key = f"token{len(token_replacers)+1}"
                value = self.parse_node(child)
                # update token replacers
                token_replacers[key] = value
                # update string format
                string_format += f"{{{key}}}"
            elif child.data == 'text':
                string_format += str(child.children[0].value)
            else:
                raise RuntimeError
            
        if len(token_replacers) == 0:
            raise RuntimeError
        
        return j.StringFormatter(string_format, **token_replacers)
            
    def parse_node(self, node: Tree | Token) -> Any:
        if isinstance(node, Token):
            if node.type == 'TRUE':
                return True
            elif node.type == 'FALSE':
                return False
            elif node.type == 'NULL':
                return None
            elif node.type == 'SIGNED_NUMBER':
                if '.' in node.value:
                    return float(node.value)
                return int(node.value)
            else:
                return str(node.value)
        
        t_name = str(node.data)

        # entry points
        if t_name == 'the_text':
            return self.generate_string_formatter(node)
        if t_name == 'javascript':
            return self.parse_node(node.children[0])
        if t_name == 'text':
            return str(node.value)
        
        # objects
        elif t_name == 'input':
            i_name = self.parse_node(node.children[0])
            return j.InputSelector(i_name)
        elif t_name == 'rt_outdir':
            return '.'
        elif t_name == 'rt_tmpdir':
            return '.'
        elif t_name == 'rt_outdir_size':
            return 1024
        elif t_name == 'rt_tmpdir_size':
            return 1024
        elif t_name == 'rt_cores':
            return j.CpuSelector()
        elif t_name == 'rt_ram':
            return j.MemorySelector()
        
        # logial - misc
        elif t_name == 'ternary':
            cond_check = self.parse_node(node.children[0])
            cond_true = self.parse_node(node.children[1])
            cond_false = self.parse_node(node.children[2])
            return j.If(cond_check, cond_true, cond_false)
        elif t_name == 'return_ifelse':
            cond_check = self.parse_node(node.children[0])
            cond_true = self.parse_node(node.children[1])
            cond_false = self.parse_node(node.children[2])
            return j.If(cond_check, cond_true, cond_false)
        elif t_name == 'return_inline':
            return self.parse_node(node.children[0])
        elif t_name == 'group':
            inner = self.parse_node(node.children[0])
            return j.GroupOperator(inner)
        
        # logial - math
        elif t_name in self.math_map:
            inner = self.parse_node(node.children[0])
            return self.math_map[t_name](inner)

        # logial - two_value
        elif t_name in self.two_value_map:
            left = self.parse_node(node.children[0])
            right = self.parse_node(node.children[1])
            return self.two_value_map[t_name](left, right)

        # attributes
        elif t_name in self.file_attr_map:
            file_object = self.parse_node(node.children[0])
            return self.file_attr_map[t_name](file_object)

        # methods - array
        elif t_name == 'meth_join':
            arr_object = self.parse_node(node.children[0])
            separator = self.parse_node(node.children[1])
            return j.JoinOperator(arr_object, separator)
        elif t_name == 'meth_slice':
            arr_object = self.parse_node(node.children[0])
            if isinstance(node.children[1], Token):
                start = int(self.parse_node(node.children[1]))
                end = None
            elif len(node.children[1].children) == 3:
                start = int(self.parse_node(node.children[1].children[0]))
                end = int(self.parse_node(node.children[1].children[2]))
            else:
                raise RuntimeError
            return j.SliceOperator(arr_object, start, end)
        elif t_name == 'meth_flat':
            arr_object = self.parse_node(node.children[0])
            return j.FlattenOperator(arr_object)
        elif t_name == 'meth_index':
            arr_object = self.parse_node(node.children[0])
            index = int(self.parse_node(node.children[1]))
            return j.IndexOperator(arr_object, index)
        
        # methods - string
        elif t_name == 'meth_split':
            string_object = self.parse_node(node.children[0])
            separator = self.parse_node(node.children[1])
            return j.SplitOperator(string_object, separator)
        elif t_name == 'meth_replace':
            string_object = self.parse_node(node.children[0])
            pattern = self.parse_node(node.children[1].children[0])
            repl = self.parse_node(node.children[1].children[2])
            return j.ReplaceOperator(string_object, pattern, repl)
        
        # methods - casting
        elif t_name == 'meth_tostr':
            inner = self.parse_node(node.children[0])
            return j.AsStringOperator(inner)

        # functions
        elif t_name == 'func_parseint':
            inner = self.parse_node(node.children[0])
            return j.AsIntOperator(inner)

        else:
            raise NotImplementedError ## this should trigger warning message
        

    
            




