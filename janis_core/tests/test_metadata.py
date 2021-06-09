import unittest
from datetime import date

from janis_core.utils.metadata import Metadata


class TestMetadata(unittest.TestCase):
    def test_update(self):

        m = Metadata()
        m2 = m.update(contributors=["c2"])
        self.assertListEqual(["c2"], m.contributors)
        self.assertListEqual(["c2"], m2.contributors)

    def test_serialize(self):
        m = Metadata(contributors=["Michael Franklin"])
        d = m.get_dict({"calculateChecksum": "ofThisDictionary"})
        self.assertIn("contributors", d)
        self.assertEqual("86a02269745ef99ac7e4999f8f66491a20a2d4ba26123f70e47cf0f4", d["checksum"])

        # Check for all attributes that are none, and check they're not in the output
        for k, v in vars(m).items():
            if v is not None:
                continue
            self.assertNotIn(k, d)

    def test_date_generated(self):
        m = Metadata()
        d = m.get_dict({})
        self.assertIn("dateGenerated", d)
        self.assertEqual(date.today().strftime("%Y-%m-%d"), d["dateGenerated"])
