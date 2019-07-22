import unittest

from janis_core import Logger


class TestCommandTool(unittest.TestCase):
    def setUp(self):
        Logger.mute()

        from janis_core.unix.tools.tar import Tar

        self.tarTool = Tar()

        Logger.unmute()
