


from janis_core.ingestion.galaxy.gx.gxtool.requirements.model import Requirement
from ..Container import Container


def select_best_container_match(containers: list[Container], requirement: Requirement) -> Container:
    containers = select_version_matches(containers, requirement)
    containers = sort_by_create_date(containers)
    return containers[0]

def select_version_matches(containers: list[Container], requirement: Requirement) -> list[Container]:
    matches = [c for c in containers if c.tool_version == requirement.version]
    if matches:
        return matches
    # TODO logging.version_mismatch
    return containers

def sort_by_create_date(containers: list[Container]) -> list[Container]: 
    return sorted(containers, key=lambda x: x.timestamp, reverse=True)


# @dataclass
# class SimilarityScore:
#     score: float
#     obj: Any

# def most_similar_version(tool_data: dict[str, Any]) -> dict[str, Any]:
#     vm = VersionMatcher()
#     target_version = requirement.version
#     selected: Optional[dict[str, str]] = None

#     selected = vm.get_version_exact(tool_data, target_version)
#     if not selected:
#         selected = vm.get_version_trimmed(tool_data, target_version)
#     if not selected:
#         selected = vm.get_most_recent(tool_data)
#     return selected
#     # if not selected:
#     #     selected = vm.get_version_trimmed_inexact(query_versions, target_version)

# def select_from_available(images: list[dict[str, Any]]) -> Optional[dict[str, Any]]:
#     for image in images:
#         if image['registry_host'] == 'quay.io/':
#             return image
#     for image in images:
#         if image['registry_host'] == 'depot.galaxyproject.org/singularity/':
#             return image
#     return None



        

