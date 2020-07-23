import unittest

import janis_core.translations.nfgen as nfgen
from janis_core.translations.nfgen import ProcessScriptType


class NfProcessTests(unittest.TestCase):
    def test_process_basic(self):
        script = "echo 'Hello, Nextflow!'"
        process = nfgen.Process("doMoreThings", script=script)

        expected = """\
process doMoreThings {

  \"""
  echo 'Hello, Nextflow!'
  \"""

}
"""
        self.assertEqual(expected, process.get_string())

    def test_process_basic_shell(self):
        script = "echo 'Hello, Nextflow!'"
        process = nfgen.Process(
            "doMoreThings", script=script, script_type=ProcessScriptType.shell
        )

        expected = """\
process doMoreThings {

  shell:
    \"""
    echo 'Hello, Nextflow!'
    \"""

}
"""
        self.assertEqual(expected, process.get_string())
