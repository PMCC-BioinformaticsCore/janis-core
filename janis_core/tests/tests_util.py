from unittest import TestCase
from janis_core.utils import (
    get_value_for_hints_and_ordered_resource_tuple,
    recursive_2param_wrap,
)


class TestOrderedHints(TestCase):

    ordered_hints = [
        (
            "hinttype1",
            {
                "ht1h1": 1,
                # "ht1h2": 2,
                "ht1h3": 3,
            },
        ),
        ("hinttype2", {"ht2h1": 11, "ht2h2": 12, "ht2h3": 13}),
    ]

    def test_first_tuple(self):
        hints = {"hinttype1": "ht1h1"}
        val = get_value_for_hints_and_ordered_resource_tuple(hints, self.ordered_hints)
        self.assertEqual(1, val)

    def test_first_tuple_2(self):
        hints = {"hinttype1": "ht1h3"}
        val = get_value_for_hints_and_ordered_resource_tuple(hints, self.ordered_hints)
        self.assertEqual(3, val)

    def test_second_tuple(self):
        hints = {"hinttype2": "ht2h2"}
        val = get_value_for_hints_and_ordered_resource_tuple(hints, self.ordered_hints)
        self.assertEqual(12, val)

    def test_none(self):
        hints = {"hinttype1": "ht1h2", "hinttype2": "ht2hn"}
        val = get_value_for_hints_and_ordered_resource_tuple(hints, self.ordered_hints)
        self.assertIsNone(val)

    def test_ordered_1(self):
        hints = {"hinttype2": "h2h1", "hinttype1": "ht1h1"}
        val = get_value_for_hints_and_ordered_resource_tuple(hints, self.ordered_hints)
        self.assertEqual(1, val)

    def test_ordered_first_no_exist(self):
        hints = {"hinttype1": "ht1h2", "hinttype2": "ht2h3"}
        val = get_value_for_hints_and_ordered_resource_tuple(hints, self.ordered_hints)
        self.assertEqual(13, val)

    def test_hint_doesnt_exist(self):
        hints = {"hinttype3": "ht3h1"}
        val = get_value_for_hints_and_ordered_resource_tuple(hints, self.ordered_hints)
        self.assertIsNone(val)


class Test2ItemMethodWrap(TestCase):
    def test_0_items(self):
        items = []
        self.assertRaises(Exception, recursive_2param_wrap, "zip", items)

    def test_1_item(self):
        items = ["item"]
        self.assertRaises(Exception, recursive_2param_wrap, "zip", items)

    def test_2_items(self):
        items = ["item1", "item2"]
        self.assertEqual("zip(item1, item2)", recursive_2param_wrap("zip", items))

    def test_3_items(self):
        items = ["item1", "item2", "item3"]
        self.assertEqual(
            "zip(item1, zip(item2, item3))", recursive_2param_wrap("zip", items)
        )

    def test_4_items(self):
        items = ["item1", "item2", "item3", "item4"]
        self.assertEqual(
            "zip(item1, zip(item2, zip(item3, item4)))",
            recursive_2param_wrap("zip", items),
        )


# class testtranslate(TestCase):
#     def test1(self):
#         from janis_bioinformatics.tools.samtools.sort.sort import SamToolsSort_1_9
#
#         SamToolsSort_1_9().translate("cwl")


# from janis_pipelines import WGSGermlineGATK

# WGSGermlineGATK().translate("cwl")
