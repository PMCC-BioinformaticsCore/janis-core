


import yaml
from typing import Optional
from janis_core.ingestion.galaxy.settings.paths import GALAXY_DATATYPES_YAML

from .JanisDatatype import JanisDatatype


class DatatypeRegister:
    def __init__(self):
        self.dtype_map: dict[str, JanisDatatype] = {}

    def get(self, datatype: str) -> Optional[JanisDatatype]:
        if datatype in self.dtype_map:
            return self.dtype_map[datatype]

    def populate(self) -> None:
        """
        func loads the combined datatype yaml then converts it to dict with format as keys
        provides structue where we can search all the galaxy and janis types given what we see
        in galaxy 'format' attributes.
        """
        path = GALAXY_DATATYPES_YAML
        with open(path, 'r') as fp:
            datatypes = yaml.safe_load(fp)
        for type_data in datatypes['types']:
            janistype = self._init_type(type_data)
            # multiple keys per datatype
            self.dtype_map[type_data['format']] = janistype
            self.dtype_map[type_data['classname']] = janistype 
            if type_data['extensions']:
                for ext in type_data['extensions']:
                    self.dtype_map[ext] = janistype 

    def _init_type(self, dtype: dict[str, str]) -> JanisDatatype:
        return JanisDatatype(
            format=dtype['format'],
            source=dtype['source'],
            classname=dtype['classname'],
            extensions=dtype['extensions'],
            import_path=dtype['import_path']
        )


# SINGLETON
register = DatatypeRegister()

