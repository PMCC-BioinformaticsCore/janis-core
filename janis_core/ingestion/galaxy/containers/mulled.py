
import os 
import sys
import subprocess
from typing import Optional, Tuple

from janis_core.ingestion.galaxy.gx.gxtool import XMLToolDefinition
from janis_core.ingestion.galaxy.gx.gxtool.requirements.model import Requirement

from .Container import Container
from .fetch import fetch_online
from .fetch import DEFAULT_CONTAINER


def gen_mulled_image(xmltool: XMLToolDefinition) -> Optional[Container]:
    _check_docker_available()
    base_image_uri, other_reqs = get_base_image_uri(xmltool)
    mulled_image = create_mulled_container(xmltool, base_image_uri, other_reqs)
    return mulled_image

def get_base_image_uri(xmltool: XMLToolDefinition) -> Tuple[str, list[Requirement]]:
    # start with main requirement as base container if we can find one
    main_req = xmltool.metadata.main_requirement
    assert(main_req)
    ordered_reqs = [main_req] + [x for x in xmltool.metadata.requirements if x.name != main_req.name]

    for req in ordered_reqs:
        container = fetch_online(req)
        if container is not None:
            other_reqs = [x for x in ordered_reqs if x.name != req.name]
            return (container.uri, other_reqs)  # type: ignore
    
    return (DEFAULT_CONTAINER, xmltool.metadata.requirements)  # type: ignore

def create_mulled_container(xmltool: XMLToolDefinition, base_container_uri: str, other_reqs: list[Requirement]) -> Container:
    manager = MulledContainerGenerator(xmltool, base_container_uri, other_reqs)
    manager.run_tasks()
    uri = manager.mulled_uri
    return Container({
        'image_type': 'docker',
        'repo': uri.split('/')[-1].split(':')[0],
        'tag': uri.split('/')[-1].split(':')[1],
        'uri': uri,
        '_timestamp': '',
    })

def _check_docker_available() -> None:
    try:
        subprocess.run(['docker', 'version'], stderr=sys.stderr, stdout=sys.stdout)
    except FileNotFoundError:
        print('you have supplied "--build-galaxy-tool-images", but docker was not found.')
        print('please install docker / run the app if not on linux.')
        print('alternatively, re-run janis translate without "--build-galaxy-tool-images"')
        sys.exit(1)

# quay.io/biocontainers/bioconductor-limma:3.34.9--r3.4.1_0

class MulledContainerGenerator:
    def __init__(self, xmltool: XMLToolDefinition, base_container_uri: str, other_reqs: list[Requirement]):
        self.xmltool = xmltool
        self.base_container_uri = base_container_uri
        self.other_reqs = other_reqs
        self.success: bool = False

    @property 
    def mulled_uri(self) -> str:
        name = self.xmltool.metadata.id
        name = name.lower()
        name = name.replace('_', '-')
        return f'ppp-janis-translate-{name}-{self.xmltool.metadata.version}'

    def run_tasks(self) -> None:
        print(f'building container for {self.xmltool.metadata.id}_{self.xmltool.metadata.version}')
        print('(this may take some time)')
        print('the following conda packages will be mulled together:')
        for req in self.xmltool.metadata.requirements:
            print(f'  - {req.name}')
        self.do_build()
        print(f'built container {self.mulled_uri}')
    
    def do_build(self) -> None:
        # create command
        command: list[str] = []
        command.append('mulled-build')
        command.append('--name-override')
        command.append(f'{self.mulled_uri}')
        command.append('--verbose')
        command.append('build-and-test')
        reqs_joined = ','.join([req.name for req in self.other_reqs])
        command.append(f'{reqs_joined}')
        
        # create env
        os.environ["DEST_BASE_IMAGE"] = f"'{self.base_container_uri}'"
        my_env = os.environ.copy()
        
        # run subprocess
        subprocess.run(command, stderr=sys.stderr, stdout=sys.stdout, env=my_env)
    

