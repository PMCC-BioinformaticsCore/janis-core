
import os
import pathlib

# private

_JANIS_DATA_DIR = os.path.join(os.getcwd(), '.janis')
_JANIS_INSTALL_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
_INGEST_DATA_DIR = os.path.join(_JANIS_INSTALL_DIR, 'ingestion', 'data')
_GALAXY_DATA_DIR = os.path.join(_JANIS_INSTALL_DIR, 'ingestion', 'data', 'galaxy')
_TEST_DATA_DIR = os.path.join(_JANIS_INSTALL_DIR, 'tests', 'data')

# public

GEN_IMAGES = False
DISABLE_CONTAINER_CACHE = False
GALAXY_CONFIG = f'{_GALAXY_DATA_DIR}/galaxy_config.yaml'
DATATYPES_YAML = f'{_INGEST_DATA_DIR}/janis_types.yaml'
CONTAINER_CACHE = f'{_JANIS_DATA_DIR}/galaxy_containers/cache.json'
WRAPPER_CACHE = f'{_JANIS_DATA_DIR}/galaxy_wrappers/cache.json'   
DEFAULT_WRAPPERS_DIR = f'{_JANIS_DATA_DIR}/galaxy_wrappers'
TESTING_WRAPPERS_DIR = f'{_TEST_DATA_DIR}/galaxy/wrappers'
