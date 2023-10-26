
import os
from lark import Lark

GRAMMAR_PATH = f'{os.path.dirname(os.path.abspath(__file__))}/grammar.ebnf'

with open(GRAMMAR_PATH) as fp:
    grammar = fp.read()

data01 = '$(self === "name" ? true : false).fastq'
data02 = '$(self === "name" ? true : false).fastq$(inputs.index_name.ext)'
data1 = '$(self === "name" ? true : false)'
data2 = '$((inputs.index_name !== null) ? inputs.index_name : inputs.sequences.nameroot)'
data3 = '$(inputs.reads_1.basename).trimmed$(inputs.reads_1.nameext)'


json_parser = Lark(grammar, start='the_text')
tree = json_parser.parse(data01)
print(tree.pretty())
tree = json_parser.parse(data02)
print(tree.pretty())
tree = json_parser.parse(data1)
print(tree.pretty())
tree = json_parser.parse(data2)
print(tree.pretty())
tree = json_parser.parse(data3)
print(tree.pretty())
print()





