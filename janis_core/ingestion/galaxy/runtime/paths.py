
import os

package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGE_DATA_DIR = f'{package_dir}/data'
JANIS_DATA_DIR = '.janis'

GALAXY_CONFIG = f'{PACKAGE_DATA_DIR}/galaxy_config.yaml'
GALAXY_DATATYPES_YAML = f'{PACKAGE_DATA_DIR}/janis_types.yaml'
CONTAINER_CACHE = f'{JANIS_DATA_DIR}/galaxy_containers/cache.json'
WRAPPER_CACHE = f'{JANIS_DATA_DIR}/galaxy_wrappers/cache.json'   
DOWNLOADED_WRAPPERS_DIR = f'{JANIS_DATA_DIR}/galaxy_wrappers'

