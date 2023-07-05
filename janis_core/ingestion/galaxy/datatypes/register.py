


import yaml
from typing import Optional
from janis_core.settings.ingest.galaxy import GALAXY_DATATYPES_YAML

from .JanisDatatype import JanisDatatype


class DatatypeRegister:
    def __init__(self):
        self.format_map: dict[str, JanisDatatype] = {}
        self.extension_map: dict[str, JanisDatatype] = {}

    def get_from_extension(self, extension: str) -> Optional[JanisDatatype]:
        if extension in self.extension_map:
            return self.extension_map[extension]

    def get_from_format(self, format: str) -> Optional[JanisDatatype]:
        if format in self.format_map:
            return self.format_map[format]

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
            self.format_map[janistype.format] = janistype
            if janistype.extensions is not None:
                for ext in janistype.extensions.split(','):
                    self.extension_map[ext] = janistype

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

