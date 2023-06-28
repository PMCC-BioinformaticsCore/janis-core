


from abc import ABC, abstractmethod
from typing import Any

from ..Container import Container



class Fetcher(ABC):
    """
    searched for the most recent container build which matches the name and version
    can just format a container requirement, or use an API to search online
    individual builds not supported, except in the case of of 'container' requirements as it just directly grabs the url from the requirement text 
    """

    @abstractmethod
    def fetch(self, requirement: Any) -> list[Container]:
        """finds relevant biocontainers for a given tool requirement"""
        ...


