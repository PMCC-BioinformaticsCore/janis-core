from unittest import TestCase
from typing import List, Union, Optional

from janis_core.utils.generics_util import (
    is_generic,
    is_base_generic,
    is_qualified_generic,
)


class TestTypingGenerics(TestCase):
    def test_unqualified_generic_list(self):
        self.assertTrue(is_base_generic(List))
        self.assertFalse(is_qualified_generic(List))
        self.assertTrue(is_generic(List))

    def test_unqualified_generic_union(self):
        self.assertTrue(is_base_generic(Union))
        self.assertFalse(is_qualified_generic(Union))
        self.assertTrue(is_generic(Union))

    def test_unqualified_generic_optional(self):
        self.assertTrue(is_base_generic(Optional))
        self.assertFalse(is_qualified_generic(Optional))
        self.assertTrue(is_generic(Optional))

    def test_qualified_generic_list(self):
        self.assertFalse(is_base_generic(List[str]))
        self.assertTrue(is_qualified_generic(List[str]))
        self.assertTrue(is_generic(List[str]))

    def test_qualified_generic_union(self):
        self.assertFalse(is_base_generic(Union[str, int]))
        self.assertTrue(is_qualified_generic(Union[str, int]))
        self.assertTrue(is_generic(Union[str, int]))

    def test_qualified_generic_optional(self):
        self.assertFalse(is_base_generic(Optional[str]))
        self.assertTrue(is_qualified_generic(Optional[str]))
        self.assertTrue(is_generic(Optional[str]))
