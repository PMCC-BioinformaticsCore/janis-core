
from janis_core.ingestion.galaxy.logs import logging
import tempfile
from typing import Optional

from janis_core.ingestion.galaxy.gx.gxtool.requirements import Requirement, CondaRequirement, ContainerRequirement
from janis_core.ingestion.galaxy.paths import USER_DATA_DIR, CONTAINER_CACHE

from .ContainerCache import ContainerCache
from .Container import Container

from .fetching.Fetcher import Fetcher
from .fetching.ContainerReqFetcher import ContainerReqFetcher
#from .fetching.GA4GHFetcher import GA4GHFetcher
from .fetching.QuayIOFetcher import QuayIOFetcher
from .fetching.presets import get_images_preset

from .selection.selection import select_best_container_match




DISABLE_CACHE = False

def fetch_container(requirement: Requirement) -> Optional[str]:
    containers: list[Container] = []
    if not containers:
        containers = _fetch_cache(requirement)
    if not containers:
        containers = _fetch_presets(requirement)
    if not containers:
        containers = _fetch_online(requirement)
    if not containers:
        logging.no_container()
        return None
    else:
        container = select_best_container_match(containers, requirement)
        return container.url
    
def _fetch_cache(requirement: Requirement) -> list[Container]:
    cache = _load_cache()
    return cache.get(requirement.name)

def _load_cache() -> ContainerCache:
    if DISABLE_CACHE:
        temp = tempfile.TemporaryFile()
        cache_path = f'{tempfile.gettempdir()}/{temp.name}'
    else:
        cache_path = f'{USER_DATA_DIR}/{CONTAINER_CACHE}'
    return ContainerCache(cache_path)

def _fetch_presets(requirement: Requirement) -> list[Container]:
    return get_images_preset(requirement)

def _fetch_online(requirement: Requirement) -> list[Container]:
    strategy = _select_strategy(requirement)
    containers = strategy.fetch(requirement)
    if containers:
        _update_cache(containers)
    return containers
    
def _select_strategy(requirement: Requirement) -> Fetcher:
    match requirement:
        case CondaRequirement():
            return QuayIOFetcher()
        #case CondaRequirement():
        #    return GA4GHFetcher()
        case ContainerRequirement():
            return ContainerReqFetcher()
        case _:
            raise RuntimeError()
    
def _update_cache(containers: list[Container]) -> None:
    cache = _load_cache()
    for container in containers: 
        cache.add(container)



