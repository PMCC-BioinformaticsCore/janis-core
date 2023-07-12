

from dataclasses import dataclass

from .configfile import XMLConfigfile
from .script import XMLScript
from .metadata import XMLMetadata
from .params.param_register import XMLParamRegister
from .tests import XMLTestRegister


@dataclass
class XMLTool:
    """
    High-level component representing a tool. 
    Does not depend on lower level representations or parsing.
    Permits storing and retreiving data about the tool.
    """
    metadata: XMLMetadata
    raw_command: str
    configfiles: list[XMLConfigfile]
    scripts: list[XMLScript]
    inputs: XMLParamRegister
    outputs: XMLParamRegister
    tests: XMLTestRegister

    # def get_input(self, query: str, strategy: str='exact') -> Optional[Param]:
    #     return self.inputs.get(query.lstrip('$'), strategy=strategy)
    
    # def list_inputs(self) -> list[Param]:
    #     return self.inputs.list()

    # def get_output(self, query: str, strategy: str='exact') -> Optional[Param]:
    #     return self.outputs.get(query.lstrip('$'), strategy=strategy)
    
    # def list_outputs(self) -> list[XMLOutputParam]:
    #     return self.outputs.list()

    # def list_tests(self) -> list[TTestCase]:
    #     return self.tests.list()

    # def get_requirements(self) -> list[Requirement]:
    #     return self.metadata.requirements
    
    # def get_main_requirement(self) -> Requirement:
    #     return self.metadata.get_main_requirement()













