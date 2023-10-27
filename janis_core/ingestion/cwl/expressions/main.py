

import os
from lark import Lark
from lark import Tree
from lark import Token
from typing import Tuple, Any
import janis_core as j

GRAMMAR_PATH = f'{os.path.dirname(os.path.abspath(__file__))}/grammar.ebnf'


def parse_expression(expr: str):
    return ExpressionParser().parse(expr)

class ExpressionParser:
    with open(GRAMMAR_PATH) as fp:
        grammar = fp.read()
    parser = Lark(grammar, start='the_text')
    
    def parse(self, expr: str) -> Tuple[Any, bool]:
        try:
            tree = self.parser.parse(expr)
            return self.parse_node(tree), True
        except Exception as e:
            return None, False
    
    def parse_node(self, node: Tree) -> Any:
        t_name = str(node.data)

        if t_name == 'javascript':
            return self.parse_node(node.children[0])
        elif t_name == 'input':
            i_name = str(node.children[0])
            return j.InputSelector(i_name)
        elif t_name == 'expr':
            pass
        
        # attributes
        elif t_name == 'attr_basename':
            file_object = self.parse_node(node.children[0])
            return j.BasenameOperator(file_object)
        elif t_name == 'attr_dirname':
            file_object = self.parse_node(node.children[0])
            return j.DirnameOperator(file_object)
        elif t_name == 'attr_nameroot':
            file_object = self.parse_node(node.children[0])
            return j.NamerootOperator(file_object)
        elif t_name == 'attr_nameext':
            file_object = self.parse_node(node.children[0])
            return j.NameextOperator(file_object)
        elif t_name == 'attr_size':
            file_object = self.parse_node(node.children[0])
            return j.FileSizeOperator(file_object)
        elif t_name == 'attr_contents':
            file_object = self.parse_node(node.children[0])
            return j.ReadContents(file_object)
        elif t_name == 'attr_length':
            file_object = self.parse_node(node.children[0])
            return j.LengthOperator(file_object)
        
        # methods
        elif t_name == 'meth_join':
            arr_object = self.parse_node(node.children[0])
            separator = str(node.children[1])
            return j.JoinOperator(arr_object, separator)
        elif t_name == 'meth_slice':
            arr_object = self.parse_node(node.children[0])
            if isinstance(node.children[1], Token):
                start = int(str(node.children[1]))
                end = None
            else:
                start = int(str(node.children[1].children[0]))
                end = int(str(node.children[1].children[2]))
            return j.SliceOperator(arr_object, start, end)
        elif t_name == 'meth_flat':
            arr_object = self.parse_node(node.children[0])
            return j.FlattenOperator(arr_object)
        elif t_name == 'meth_index':
            arr_object = self.parse_node(node.children[0])
            index = int(str(node.children[1]))
            return j.IndexOperator(arr_object, index)
        elif t_name == 'meth_split':
            string_object = self.parse_node(node.children[0])
            separator = str(node.children[1])
            return j.SplitOperator(string_object, separator)
        elif t_name == 'meth_replace':
            string_object = self.parse_node(node.children[0])
            pattern = str(node.children[1].children[0])
            repl = str(node.children[1].children[2])
            return j.ReplaceOperator(string_object, pattern, repl)
        
        else:
            pass
            
        raise NotImplementedError




