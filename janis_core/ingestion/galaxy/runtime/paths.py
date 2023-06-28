
import os

package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

USER_DATA_DIR = f'{package_dir}/data'
GALAXY_CONFIG = f'{USER_DATA_DIR}/galaxy_config.yaml'
GALAXY_DATATYPES_YAML = f'{USER_DATA_DIR}/galaxy_janis_types.yaml'
CONTAINER_CACHE = f'{USER_DATA_DIR}/container_url_cache.json'
WRAPPER_CACHE = f'{USER_DATA_DIR}/wrappers.json'
DOWNLOADED_WRAPPERS_DIR = f'.janis/galaxy_wrappers'

