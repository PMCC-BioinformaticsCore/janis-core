import unittest

from janis_core import InputSelector
from janis_core.operators.stringformatter import StringFormatter
from janis_core.utils.bracketmatching import (
    get_keywords_between_braces,
    variable_name_validator,
)
from janis_core.utils.errors import (
    IncorrectArgsException,
    TooManyArgsException,
    ConflictingArgumentsException,
    InvalidByProductException,
)


class TestMatchValidator(unittest.TestCase):
    # Use cases from: https://stackoverflow.com/a/36331242/2860731

    def test_valid_1(self):
        self.assertTrue(variable_name_validator("x"))

    def test_valid_2(self):
        self.assertTrue(variable_name_validator("X"))

    def test_valid_3(self):
        self.assertTrue(variable_name_validator("X123"))

    def test_invalid_1(self):
        self.assertFalse(variable_name_validator("2"))

    def test_invalid_2(self):
        self.assertFalse(variable_name_validator("while"))


class TestBraceDetection(unittest.TestCase):
    def test_no_group(self):
        k = "there's no match here"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(0, len(res))

    def test_one_group(self):
        k = "there's {one} match here"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(1, len(res))
        self.assertSetEqual({"one"}, res)

    def test_two_groups(self):
        k = "there's {one} matches {here}"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(2, len(res))
        self.assertSetEqual({"one", "here"}, res)

    def test_nested_group(self):
        k = "theres {a nested {group} here} but not just {that}"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(1, len(res))
        self.assertSetEqual({"that"}, res)

    def test_unwrap_first_group(self):
        k = "theres } the wrong bracket {first}"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(1, len(res))
        self.assertSetEqual({"first"}, res)


class TestMatchDetectionAndValidation(unittest.TestCase):
    # Use cases from: https://stackoverflow.com/a/36331242/2860731

    def test_valid_1(self):
        k = "test {x}"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(1, len(res))
        self.assertTrue("x" in res)

    def test_valid_2(self):
        k = "test {X}"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(1, len(res))
        self.assertTrue("X" in res)

    def test_valid_3(self):
        k = "test {X123}"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(1, len(res))
        self.assertTrue("X123" in res)

    def test_invalid_1(self):
        k = "test {2}"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(0, len(res))
        self.assertFalse("2" in res)

    def test_invalid_2(self):
        k = "test {while}"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(0, len(res))
        self.assertFalse("while" in res)

    def test_combination(self):
        k = "The format {is} {valid} {while} some of the args are invalid"
        res, _ = get_keywords_between_braces(k)
        self.assertEqual(1, len(res))
        self.assertFalse("is" in res)
        self.assertTrue("valid" in res)
        self.assertFalse("while" in res)


class TestStringFormatterInit(unittest.TestCase):
    def test_no_format(self):
        b = StringFormatter("this _format has no arguments")
        self.assertDictEqual({}, b.kwargs)

    def test_one_correct_argument(self):
        b = StringFormatter("this {test} has one argument", test="My test _format")
        self.assertDictEqual({"test": "My test _format"}, b.kwargs)

    def test_one_incorrect_argument(self):
        self.assertRaises(
            IncorrectArgsException,
            StringFormatter,
            format="this {other} has one argument",
            test="this will fail",
        )

    def test_too_many_args(self):
        self.assertRaises(
            TooManyArgsException,
            StringFormatter,
            format="this {other} has one argument",
            other="this will fail",
            extra="because",
            it="will",
        )

    def test_too_many_args2(self):
        self.assertRaises(
            TooManyArgsException,
            StringFormatter,
            format="this should fail",
            extra="will",
        )


class TestStringFormatterConcat(unittest.TestCase):
    def test_formatter_and_string(self):
        b = StringFormatter("no args") + " here"
        self.assertEqual("no args here", b._format)
        self.assertEqual(0, len(b.kwargs))

    def test_two_formatters_no_format(self):
        b = StringFormatter("no args") + StringFormatter(" in this")
        self.assertEqual("no args in this", b._format)
        self.assertEqual(0, len(b.kwargs))

    def test_two_formatter_no_overlap(self):
        b = StringFormatter("one {arg} ", arg="hi") + StringFormatter(
            "a {different} arg", different="diff"
        )
        self.assertEqual("one {arg} a {different} arg", b._format)
        self.assertEqual(2, len(b.kwargs))
        self.assertEqual("hi", b.kwargs["arg"])
        self.assertEqual("diff", b.kwargs["different"])

    def test_two_formatters_equal_overlap(self):
        b = StringFormatter("one {arg} ", arg="hi") + StringFormatter(
            "and the same {arg}", arg="hi"
        )
        self.assertEqual("one {arg} and the same {arg}", b._format)
        self.assertEqual(1, len(b.kwargs))
        self.assertEqual("hi", b.kwargs["arg"])

    def test_two_formatters_nonequal_overlap(self):
        try:
            b = StringFormatter("one {arg} ", arg="hi") + StringFormatter(
                "and the same {arg}", arg="hi2"
            )
            self.fail(
                "Test 'test_two_formatters_nonequal_overlap' should fail as the two args have different values"
            )
        except ConflictingArgumentsException:
            self.assertTrue(True)

    def test_byproduct_by_concat_string(self):
        try:
            b = StringFormatter("one {") + "test}"
            self.fail(
                "Test 'test_byproduct_by_concat_string' should fail as concatting the string creates a new placeholder"
            )

        except InvalidByProductException:
            self.assertTrue(True)

    def test_byproduct_by_concat_string_covered(self):
        b = StringFormatter("one {arg} + another {", arg="q") + "arg}"
        self.assertEqual("one {arg} + another {arg}", b._format)
        self.assertEqual(1, len(b.kwargs))

    def test_reverse_add(self):
        b = "Hello, " + StringFormatter("world")
        self.assertEqual("Hello, world", b._format)


class TestStringFormatterResolve(unittest.TestCase):
    def test_formatter_no_replacement(self):
        b = StringFormatter("no args")
        resolved = b.resolve_with_resolved_values()
        self.assertEqual("no args", resolved)

    def test_formatter_one_replacement(self):
        b = StringFormatter("one {arg} was resolved", arg=None)
        resolved = b.resolve_with_resolved_values(arg="value")
        self.assertEqual("one value was resolved", resolved)

    def test_formatter_two_replacement(self):
        b = StringFormatter(
            "{howmany} argument(s) {wasORwere} resolved", howmany=None, wasORwere=None
        )
        resolved = b.resolve_with_resolved_values(howmany=1, wasORwere="was")
        self.assertEqual("1 argument(s) was resolved", resolved)

    def test_formatter_duplicate_replacement(self):
        b = StringFormatter("{arg} is the same as {arg}", arg=None)
        resolved = b.resolve_with_resolved_values(
            arg="S07E25"
        )  # ;) https://www.youtube.com/watch?v=7WCfTREZSdQ
        self.assertEqual("S07E25 is the same as S07E25", resolved)

    def test_evaluate_from_input_selector_resolved(self):
        namesel = InputSelector("name")
        adjsel = InputSelector("adjsel")
        b = StringFormatter("{name} is {adjective}", name=namesel, adjective=adjsel)
        inp = {"name": "Janis", "adjsel": "pretty good"}
        self.assertEqual("Janis is pretty good", b.evaluate(inp))


class TestInputSelectorConversion(unittest.TestCase):
    def test_input_selector_conversion(self):
        inpsel = InputSelector("test")
        formatter = inpsel.to_string_formatter()
        self.assertEqual(1, len(formatter.kwargs))
        self.assertEqual(inpsel, formatter.kwargs.get("test"))


class TestDictionaryCombinationBuilder(unittest.TestCase):
    def test_double(self):
        d = {"a": [1, 2], "b": ["a", "b", "c"]}
        combinations = StringFormatter.generate_combinations_of_input_dicts(
            list(d.items())
        )
        self.assertEqual(6, len(combinations))
        self.assertDictEqual({"a": 1, "b": "a"}, combinations[0])
        self.assertDictEqual({"a": 1, "b": "b"}, combinations[1])
        self.assertDictEqual({"a": 1, "b": "c"}, combinations[2])
        self.assertDictEqual({"a": 2, "b": "a"}, combinations[3])
        self.assertDictEqual({"a": 2, "b": "b"}, combinations[4])
        self.assertDictEqual({"a": 2, "b": "c"}, combinations[5])

    def test_multi_string_formatter(self):
        sf = StringFormatter("iteration_{i}", i=InputSelector("iters"))
        sfs = sf.evaluate({"iters": [1, 2, 3]})

        self.assertEqual(3, len(sfs))
        self.assertEqual("iteration_1", sfs[0])
        self.assertEqual("iteration_2", sfs[1])
        self.assertEqual("iteration_3", sfs[2])

    def test_dotproduct_string_formatter(self):
        sf = StringFormatter(
            "iteration_{i}_{j}", i=InputSelector("it1"), j=InputSelector("it2")
        )
        sfs = sf.evaluate({"it1": [1, 2, 3], "it2": ["a", "b", "c"]})

        self.assertEqual(3, len(sfs))
        self.assertEqual("iteration_1_a", sfs[0])
        self.assertEqual("iteration_2_b", sfs[1])
        self.assertEqual("iteration_3_c", sfs[2])

    def test_crossproduct_string_formatter(self):
        sf = StringFormatter(
            "iteration_{i}_{j}", i=InputSelector("it1"), j=InputSelector("it2")
        )
        sfs = sf.evaluate({"it1": [1, 2, 3], "it2": ["a", "b"]})

        self.assertEqual(6, len(sfs))
        self.assertEqual("iteration_1_a", sfs[0])
        self.assertEqual("iteration_1_b", sfs[1])
        self.assertEqual("iteration_2_a", sfs[2])
        self.assertEqual("iteration_2_b", sfs[3])
        self.assertEqual("iteration_3_a", sfs[4])
        self.assertEqual("iteration_3_b", sfs[5])
