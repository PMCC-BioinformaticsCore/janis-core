

from abc import ABC
from janis_core import CommandTool, PythonTool


BIOINFORMATICS_MODULE = "bioinformatics"

class BioinformaticsTool(CommandTool, ABC):
    def tool_module(self):
        return BIOINFORMATICS_MODULE

class BioinformaticsPythonTool(PythonTool, ABC):
    def tool_module(self):
        return BIOINFORMATICS_MODULE
    
class UnixTool(CommandTool, ABC):
    def tool_module(self):
        return "unix"

    def container(self):
        return "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"

    def version(self):
        return "v1.0.0"

    

