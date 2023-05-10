
import os 
import sys
import subprocess
from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition
from janis_core.ingestion.galaxy.gx.gxtool.requirements.model import Requirement


def create_mulled_container(xmltool: XMLToolDefinition, base_container_uri: str, other_reqs: list[Requirement]) -> str:
    manager = MulledContainerGenerator(xmltool, base_container_uri, other_reqs)
    manager.run_tasks()
    return manager.mulled_uri 

class MulledContainerGenerator:
    def __init__(self, xmltool: XMLToolDefinition, base_container_uri: str, other_reqs: list[Requirement]):
        self.xmltool = xmltool
        self.base_container_uri = base_container_uri
        self.other_reqs = other_reqs

    @property 
    def mulled_uri(self) -> str:
        name = self.xmltool.metadata.id
        name = name.lower()
        name = name.replace('_', '-')
        return f'ppp-janis-translate-{name}-{self.xmltool.metadata.version}'

    def run_tasks(self) -> None:
        print(f'building container for {self.xmltool.metadata.id}_{self.xmltool.metadata.version}')
        print('the following conda packages will be mulled together:')
        for req in self.xmltool.metadata.requirements:
            print(f'  - {req.name}')
        print('(this may take some time)')
        self.do_build()
        print(f'built container {self.mulled_uri}')
    
    def do_build(self) -> None:
        os.environ["DEST_BASE_IMAGE"] = f"'{self.base_container_uri}'"
        command: list[str] = []
        command.append('mulled-build')
        command.append('--verbose')
        command.append('build-and-test')
        reqs_joined = ','.join([req.name for req in self.other_reqs])
        command.append(f"'{reqs_joined}'")
        result = subprocess.run(command, stderr=sys.stderr, stdout=sys.stdout)
        print(result)
    

