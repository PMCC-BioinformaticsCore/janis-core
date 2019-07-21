import unittest

from core import Logger


class TestCommandTool(unittest.TestCase):
    def setUp(self):
        Logger.mute()

        from core.unix.tools.tar import Tar

        self.tarTool = Tar()

        Logger.unmute()
