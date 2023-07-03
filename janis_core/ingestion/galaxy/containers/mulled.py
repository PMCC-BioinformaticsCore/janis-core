
import os
import sys
import subprocess
from typing import Optional, Tuple

from janis_core.ingestion.galaxy.gxtool.model import XMLRequirement
from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from .Container import Container
from .fetch import fetch_online
from .fetch import DEFAULT_CONTAINER


def gen_mulled_image(xmltool: XMLTool) -> Optional[Container]:
    _check_docker_available()
    base_image_uri, other_reqs = get_base_image_uri(xmltool)
    mulled_image = create_mulled_container(xmltool, base_image_uri, other_reqs)
    return mulled_image

def get_base_image_uri(xmltool: XMLTool) -> Tuple[str, list[XMLRequirement]]:
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

def create_mulled_container(xmltool: XMLTool, base_container_uri: str, other_reqs: list[XMLRequirement]) -> Container:
    # user messages
    print(f'building container for {xmltool.metadata.id}_{xmltool.metadata.version}')
    print('the following conda packages will be mulled together:')
    for req in xmltool.metadata.requirements:
        print(f'  - {req.name}')
    print('(this may take some time)')

    # doing mulled-build
    manager = MulledContainerGenerator(xmltool, base_container_uri, other_reqs)
    manager.do_build()

    # managing result
    if manager.success:
        uri = manager.mulled_uri
        print(f'built container {manager.mulled_uri}')
    else:
        uri = base_container_uri
        print('mulled container generation failed. falling back to base container.')

    # returning container
    return Container({
        'image_type': 'docker',
        'repo': uri.split('/')[-1].split(':')[0],
        'tag': uri.split('/')[-1].split(':')[1],
        'uri': uri,
        '_timestamp': '',
    })

def _check_docker_available() -> None:
    try:
        completed_process = subprocess.run(['docker', 'version'], shell=True, capture_output=True)
        if completed_process.returncode == 0:
            print('docker found. proceeding with mulled container generation.')
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        print('you have supplied "--build-galaxy-tool-images", but docker was not found.')
        print('please install docker / run the app if not on linux.')
        print('alternatively, re-run janis translate without "--build-galaxy-tool-images"')
        sys.exit(1)

# quay.io/biocontainers/bioconductor-limma:3.34.9--r3.4.1_0

class MulledContainerGenerator:
    def __init__(self, xmltool: XMLTool, base_container_uri: str, other_reqs: list[XMLRequirement]):
        self.xmltool = xmltool
        self.base_container_uri = base_container_uri
        self.other_reqs = other_reqs
        self.success: bool = False

    @property
    def mulled_name(self) -> str:
        version = self.xmltool.metadata.version
        name = self.xmltool.metadata.id
        name = name.lower()
        name = name.replace('_', '-')
        return f'{name}:{version}'

    @property
    def mulled_uri(self) -> str:
        return f'quay.io/ppp-janis-translate/{self.mulled_name}'

    def do_build(self) -> None:
        # gen command
        command: list[str] = []
        command.append(f"DEST_BASE_IMAGE='{self.base_container_uri}'")
        command.append("mulled-build")
        command.append("-c")
        command.append("iuc,conda-forge,bioconda")
        command.append("--namespace")
        command.append("'ppp-janis-translate'")
        command.append("--name-override")
        command.append(f"'{self.mulled_name}'")
        command.append("--verbose")
        command.append("build")
        extra_reqs = ','.join([req.name for req in self.other_reqs])
        command.append(f"'{extra_reqs}'")
        cmdstr = ' '.join(command)

        # run subprocess
        completed_process = subprocess.run(cmdstr, shell=True, capture_output=True)
        if completed_process.returncode == 0:
            self.success = True
        else:
            self.success = False
        
        # stderr = completed_process.stderr.decode('utf-8')
        # stdout = completed_process.stdout.decode('utf-8')
        # print('\n--- STDERR ---')
        # print(stderr)
        # print('\n--- STDOUT ---')
        # print(stdout)
        # print()


