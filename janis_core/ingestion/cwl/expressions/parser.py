
import os
from lark import Lark


GRAMMAR_PATH = f'{os.path.dirname(os.path.abspath(__file__))}/grammar.ebnf'

with open(GRAMMAR_PATH) as fp:
    grammar = fp.read()

data1 = '${ return 100}'
data2 = '${ return(parseInt(runtime.ram/runtime.cores-100).toString() + "M") }'
data3 = """${
  if( inputs.output_name == null ){
    return inputs.bedgraph.basename;
  }
  else{
    return inputs.output_name;
  }
}"""

json_parser = Lark(grammar, start='the_text')
tree = json_parser.parse(data1)
print(tree.pretty())
tree = json_parser.parse(data2)
print(tree.pretty())
tree = json_parser.parse(data3)
print(tree.pretty())
print()




