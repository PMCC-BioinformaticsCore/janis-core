
from abc import ABC
from janis_core import Workflow

 
BIOINFORMATICS_MODULE = "bioinformatics"

class BioinformaticsWorkflow(Workflow, ABC):
    def tool_module(self):
        return BIOINFORMATICS_MODULE