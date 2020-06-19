# import unittest
# import tempfile
# import logging
#
# from janis_core.translations import CwlTranslator
# from janis_core.types import Stdout
#
# from janis_core import CommandToolBuilder, ToolInput, ToolOutput, Tool
#
#
# @unittest.skipUnless(
#     False, reason="MiniWDL + cwltool dependency haven't been properly resolved"
# )
# class TestJanisFunctions(unittest.TestCase):
#     @classmethod
#     def setUpClass(cls):
#         import WDL
#
#         logging.basicConfig(
#             level=logging.DEBUG, format="%(name)s %(levelname)s %(message)s"
#         )
#         logger = logging.getLogger(cls.__name__)
#         cfg = WDL.runtime.config.Loader(logger, [])
#         WDL.runtime.task.SwarmContainer.global_init(cfg, logger)
#
#     def setUp(self):
#         self._dir = tempfile.mkdtemp(prefix="miniwdl_test_stdlib_")
#
#     def _test_wdl_task(
#         self, wdl: str, inputs=None, expected_exception: Exception = None, cfg=None
#     ):
#         """
#          Source: MiniWDL project (MIT License)
#
#          https://github.com/chanzuckerberg/miniwdl
#          (Based on: https://github.com/chanzuckerberg/miniwdl/blob/481ead80cac4d765979998f0e7959c889bd5dc75/tests/test_5stdlib.py#L21-L43)
#         """
#         import WDL
#
#         cfg = cfg or WDL.runtime.config.Loader(logging.getLogger(self.id()), [])
#         try:
#             doc = WDL.parse_document(wdl)
#             assert len(doc.tasks) == 1
#             doc.typecheck()
#             assert (
#                 len(
#                     doc.tasks[0].required_inputs.subtract(doc.tasks[0].available_inputs)
#                 )
#                 == 0
#             )
#             if isinstance(inputs, dict):
#                 inputs = WDL.values_from_json(
#                     inputs, doc.tasks[0].available_inputs, doc.tasks[0].required_inputs
#                 )
#             rundir, outputs = WDL.runtime.run(
#                 cfg, doc.tasks[0], (inputs or WDL.Env.Bindings()), run_dir=self._dir
#             )
#         except WDL.runtime.RunFailed as exn:
#             if expected_exception is not None:
#                 self.assertIsInstance(exn.__context__, expected_exception)
#                 return exn.__context__
#             raise exn.__context__
#         except Exception as exn:
#             if expected_exception is not None:
#                 self.assertIsInstance(exn, expected_exception)
#                 return exn.__context__
#             raise
#         if expected_exception is not None:
#             self.assertFalse(str(expected_exception) + " not raised")
#         return WDL.values_to_json(outputs)
#
#     def _test_cwl(self, tool: Tool):
#         from os.path import join
#         import cwltool.factory
#
#         fac = cwltool.factory.Factory()
#         outdir = tempfile.gettempdir()
#         tooldict = CwlTranslator().translate(
#             tool, export_path=outdir, to_console=True, to_disk=True
#         )
#         outname = join(outdir, CwlTranslator.tool_filename(tool))
#         tool = fac.make(outname)
#         results = tool()
#         print(results)
#
#     def test_simple(self):
#         EchoTestTool = CommandToolBuilder(
#             tool="testtool_",
#             base_command=["echo"],
#             inputs=[ToolInput("inp", str, default="Hello, World", position=0)],
#             outputs=[ToolOutput("out", Stdout)],
#             container="ubuntu",
#             version="v0.1.0",
#         )
#
#         task = EchoTestTool()
#
#         # self._test_cwl(task)
#         outputs = self._test_wdl_task(task.translate("wdl"))
#
#         with open(outputs["out"]) as f:
#             out = f.readline().strip()
#         print("OUTPUT: \n\t" + str(out))
#         self.assertEqual("Hello, World", out)
