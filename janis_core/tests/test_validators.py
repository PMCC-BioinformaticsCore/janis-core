import unittest

from janis_core.utils.validators import Validators
from janis_core import settings

class TestValidators(unittest.TestCase):

    def setUp(self) -> None:
        settings.validation.STRICT_IDENTIFIERS = True
        settings.validation.VALIDATE_STRINGFORMATTERS = True

    def test_valid_identifiers(self):
        self.assertTrue(Validators.validate_identifier("test_workflow"))

    def test_invalid_identifiers(self):
        self.assertFalse(Validators.validate_identifier("test-workflow"))

    def test_invalid_sample_name(self):
        self.assertFalse(Validators.validate_identifier("fastqs_CDG-025-156R_PDX"))

    def test_invalid_sample_name_error(self):
        error = Validators.reason_for_failure("fastqs_CDG-025-156R_PDX")
        self.assertNotEqual("Undefined", error)

    def test_transform_sample_name(self):
        self.assertEqual(
            "fastqs_CDG025156R_PDX",
            Validators.transform_identifier_to_be_valid("fastqs_CDG-025-156R_PDX"),
        )
