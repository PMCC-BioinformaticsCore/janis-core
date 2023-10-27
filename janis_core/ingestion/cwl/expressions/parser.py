
import os
from lark import Lark


GRAMMAR_PATH = f'{os.path.dirname(os.path.abspath(__file__))}/grammar.ebnf'

with open(GRAMMAR_PATH) as fp:
    grammar = fp.read()

dataa = "$(inputs.sequences.split('/'))"
data0 = "$(self.location.split('/'))"
data1 = "$(self.location.split('/').slice(-1)[0])"
data2 = '$(self === "name" ? true : false)'
data3 = "$(self === 'name' ? true : false).fastq"
data4 = '$(self === "name" ? true : false).fastq$(inputs.index_name.ext)'
data5 = '$((inputs.index_name !== null) ? inputs.index_name : inputs.sequences.nameroot)'
data6 = '$(inputs.reads_1.basename).trimmed$(inputs.reads_1.nameext)'
data7 = "$(self.location.split('/').slice(-1)[0].split('.').slice(-1)[0])"

json_parser = Lark(grammar, start='the_text')
tree = json_parser.parse(dataa)
print(tree.pretty())
tree = json_parser.parse(data0)
print(tree.pretty())
tree = json_parser.parse(data1)
print(tree.pretty())
tree = json_parser.parse(data2)
print(tree.pretty())
tree = json_parser.parse(data3)
print(tree.pretty())
tree = json_parser.parse(data4)
print(tree.pretty())
tree = json_parser.parse(data5)
print(tree.pretty())
tree = json_parser.parse(data6)
print(tree.pretty())
tree = json_parser.parse(data7)
print(tree.pretty())
print()





