
import os
import pathlib


GEN_IMAGES = False
DISABLE_CONTAINER_CACHE = False

JANIS_DATA_DIR = os.path.join(os.getcwd(), '.janis')
JANIS_INSTALL_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
PACKAGE_DATA_DIR = os.path.join(JANIS_INSTALL_DIR, 'ingestion', 'galaxy', 'data')

GALAXY_CONFIG = f'{PACKAGE_DATA_DIR}/galaxy_config.yaml'
GALAXY_DATATYPES_YAML = f'{PACKAGE_DATA_DIR}/janis_types.yaml'
CONTAINER_CACHE = f'{JANIS_DATA_DIR}/galaxy_containers/cache.json'
WRAPPER_CACHE = f'{JANIS_DATA_DIR}/galaxy_wrappers/cache.json'   
DOWNLOADED_WRAPPERS_DIR = f'{JANIS_DATA_DIR}/galaxy_wrappers'
