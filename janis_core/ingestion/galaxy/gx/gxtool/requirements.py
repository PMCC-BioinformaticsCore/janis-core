

from abc import ABC, abstractmethod
from dataclasses import dataclass


"""
EXAMPLE
<requirements>
    <requirement type="package" version="1.2.0">bowtie</requirement>
    <container type="docker">google/deepvariant:@TOOL_VERSION@</container>
</requirements>
"""


@dataclass
class Requirement(ABC):
    """models a tool XML requirement"""

    @property
    @abstractmethod
    def name(self) -> str:
        ...
    
    @property
    @abstractmethod
    def version(self) -> str:
        ...
    
    # @abstractmethod
    # def to_dict(self) -> dict[str, str]:
    #     ...



@dataclass
class CondaRequirement(Requirement):
    _name: str
    _version: str
    subtype: str = 'conda'

    @property
    def name(self) -> str:
        return self._name

    @property
    def version(self) -> str:
        return self._version
    
    # def to_dict(self) -> dict[str, str]:
    #     return {
    #         'type': self.subtype,
    #         'name': self.name,
    #         'version': self.version
    #     }


@dataclass
class ContainerRequirement(Requirement):
    url: str
    flavour: str
    registry_host: str
    image_type: str
    subtype: str = 'container'

    @property
    def name(self) -> str:
        return self.url.rsplit('/', 1)[-1].split(':', 1)[0]
    
    @property
    def version(self) -> str:
        return self.url.rsplit('/', 1)[-1].split(':', 1)[-1]

    # def to_dict(self) -> dict[str, str]:
    #     return {
    #         'type': self.subtype,
    #         'url': self._url,
    #         'flavour': self._flavour,
    #         'registry_host': self._registry_host
    #     }



