import unittest
from janis_core.types.operator import *


class TestOperators(unittest.TestCase):
    def test_not_operator(self):
        op = NotOperator(1)
        self.assertEqual("!(1)", str(op))


class TestAndOperator(unittest.TestCase):
    def test_add_operator(self):
        op = AndOperator("cond1", "cond2")
        self.assertEqual("(cond1 && cond2)", str(op))

    def test_nested_add_operator(self):
        op = AndOperator("cond1", AndOperator("cond2", "cond3"))
        self.assertEqual("(cond1 && (cond2 && cond3))", str(op))

    def test_and_to_operator(self):
        op = AndOperator("cond1", "cond2").op_and("cond3")
        self.assertEqual("((cond1 && cond2) && cond3)", str(op))


class TestAddOperator(unittest.TestCase):
    def test_add_operator(self):
        op = AddOperator(1, 2)
        self.assertEqual("(1 + 2)", str(op))

    def test_nested_add_operator(self):
        op = AddOperator(1, AddOperator(2, 3))
        self.assertEqual("(1 + (2 + 3))", str(op))

    def test_radd_to_number(self):
        op = 1 + AddOperator(2, 3)
        self.assertEqual("(1 + (2 + 3))", str(op))

    def test_add_to_number(self):
        op = AddOperator(1, 2) + 3
        self.assertEqual("((1 + 2) + 3)", str(op))
