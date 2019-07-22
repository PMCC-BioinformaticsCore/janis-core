from abc import ABC

from janis_core import CommandTool


class UnixTool(CommandTool, ABC):
    @staticmethod
    def tool_module():
        return "unix"

    @staticmethod
    def docker():
        return "ubuntu:latest"
