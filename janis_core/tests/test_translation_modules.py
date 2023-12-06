import unittest
from os import getcwd
import os.path

from janis_core.translation_deps.exportpath import ExportPathKeywords
from janis_core.ingestion import ingest
from janis_core.translations import CwlTranslator
from janis_core.translations import NextflowTranslator
from janis_core.translations import WdlTranslator


class TestExportPath(unittest.TestCase):
    """
    These tests are for components that apply to every single translation.
    Perhaps the export path function, or anything in translationbase.py.
    """

    def test_user_path(self):
        self.assertEqual(
            os.path.expanduser("~"), ExportPathKeywords.resolve("~", None, None)
        )

    def test_workflow_spec(self):
        self.assertEqual(
            "/my/path/to/cwl",
            ExportPathKeywords.resolve("/my/path/to/{language}", "cwl", None),
        )

    def test_workflow_name(self):
        self.assertEqual(
            "/my/workflow_name/path",
            ExportPathKeywords.resolve("/my/{name}/path", None, "workflow_name"),
        )

    def test_multi_replace(self):
        self.assertEqual(
            "/test_multi_replace/test_multi_replace/test_multi_replace",
            ExportPathKeywords.resolve(
                "/{name}/{name}/{name}", None, "test_multi_replace"
            ),
        )

    def test_combo_replace(self):
        self.assertEqual(
            os.path.expanduser("~") + "/Desktop/workflowname/wdl/",
            ExportPathKeywords.resolve(
                "~/Desktop/{name}/{language}/", "wdl", "workflowname"
            ),
        )

    def test_replace_only_cwd_dot(self):
        self.assertEqual(getcwd(), ExportPathKeywords.resolve(".", None, None))

    def test_replace_none(self):
        self.assertEqual(getcwd(), ExportPathKeywords.resolve(None, None, None))

    def test_replace_empty(self):
        self.assertEqual(getcwd(), ExportPathKeywords.resolve("", None, None))

    def test_replace_falsey(self):
        self.assertEqual(getcwd(), ExportPathKeywords.resolve(False, None, None))

    def test_replace_cwd_in_scope(self):
        self.assertEqual(
            os.path.join(getcwd(), "test"),
            ExportPathKeywords.resolve("./test", None, None),
        )

    def test_replace_dontreplace(self):
        path = "/.myfile/starting/with/dot"
        self.assertEqual(path, ExportPathKeywords.resolve(path, None, None))

    def test_random_dot(self):
        path = "/mypath/./starting/with/dot"
        self.assertEqual(path, ExportPathKeywords.resolve(path, None, None))

    def test_no_preceding_slash(self):
        path = "test/"
        self.assertEqual(
            os.path.join(getcwd(), path), ExportPathKeywords.resolve(path, None, None)
        )

    def test_no_spec_except(self):
        self.assertRaises(
            Exception,
            ExportPathKeywords.resolve,
            path="{language}",
            workflow_spec=None,
            workflow_name="name",
        )

    def test_no_name_except(self):
        self.assertRaises(
            Exception,
            ExportPathKeywords.resolve,
            path="{name}",
            workflow_spec="spec",
            workflow_name=None,
        )

