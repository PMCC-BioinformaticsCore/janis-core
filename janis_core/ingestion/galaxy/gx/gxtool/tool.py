

from dataclasses import dataclass

from ..configfiles.Configfile import Configfile
from .metadata import ToolXMLMetadata
from .param.ParamRegister import ParamRegister
from .TestRegister import TestRegister
from .requirements import ContainerRequirement, CondaRequirement

Requirement = ContainerRequirement | CondaRequirement


@dataclass
class XMLToolDefinition:
    """
    High-level component representing a tool. 
    Does not depend on lower level representations or parsing.
    Permits storing and retreiving data about the tool.
    """
    metadata: ToolXMLMetadata
    raw_command: str
    configfiles: list[Configfile]
    inputs: ParamRegister
    outputs: ParamRegister
    tests: TestRegister

    # def get_input(self, query: str, strategy: str='exact') -> Optional[Param]:
    #     return self.inputs.get(query.lstrip('$'), strategy=strategy)
    
    # def list_inputs(self) -> list[Param]:
    #     return self.inputs.list()

    # def get_output(self, query: str, strategy: str='exact') -> Optional[Param]:
    #     return self.outputs.get(query.lstrip('$'), strategy=strategy)
    
    # def list_outputs(self) -> list[OutputParam]:
    #     return self.outputs.list()

    # def list_tests(self) -> list[TTestCase]:
    #     return self.tests.list()

    # def get_requirements(self) -> list[Requirement]:
    #     return self.metadata.requirements
    
    # def get_main_requirement(self) -> Requirement:
    #     return self.metadata.get_main_requirement()













