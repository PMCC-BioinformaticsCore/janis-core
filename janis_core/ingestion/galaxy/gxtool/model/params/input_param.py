

from dataclasses import dataclass
from typing import Any, Optional
from .param import XMLParam


class XMLInputParam(XMLParam):
    def __init__(self, name: str):
        self.name: str = name
        self.formats: list[str] = []
        self.label: str = ''
        self.argument: Optional[str] = None
        self.helptext: str = ''
        self._optional: bool = False

    def set_optionality(self, value: bool) -> None:
        self._optional = value

    @property
    def default(self) -> Any:
        return None

    @property
    def docstring(self) -> str:
        return self.generic_get_docstring()

    def generic_get_docstring(self) -> str:
        if self.helptext and self.label:
            return self.label + '. ' + self.helptext
        elif self.helptext:
            return self.helptext
        elif self.label:
            return self.label
        return ''
    
    @property
    def optional(self) -> bool:
        return self._optional
    
    @property
    def array(self) -> bool:
        return False



class XMLTextParam(XMLInputParam):
    def __init__(self, name: str):
        super().__init__(name)
        self.value: Optional[str] = None
        self.formats: list[str] = ['string']

    @property
    def default(self) -> Any:
        if self.value:
            return self.value
        return None
    

class XMLIntegerParam(XMLInputParam):
    def __init__(self, name: str):
        super().__init__(name)
        self.value: Optional[str] = None
        self.min: Optional[str] = None
        self.max: Optional[str] = None
        self.formats: list[str] = ['integer']

    @property
    def default(self) -> Any:
        if self.value:
            return self.value
        # elif self.min:
        #     return self.min
        # elif self.max:
        #     return self.max
        return None


class XMLFloatParam(XMLInputParam):
    def __init__(self, name: str):
        super().__init__(name)
        self.value: Optional[str] = None
        self.min: Optional[str] = None
        self.max: Optional[str] = None
        self.formats: list[str] = ['float']

    @property
    def default(self) -> Any:
        if self.value:
            return self.value
        elif self.min:
            return self.min
        elif self.max:
            return self.max
        return None


class XMLBoolParam(XMLInputParam):
    def __init__(self, name: str):
        super().__init__(name)
        self.checked: bool = True
        self.truevalue: str = ''
        self.falsevalue: str = ''
        self.formats: list[str] = ['boolean']

    @property
    def default(self) -> str:
        if self.checked:
            return self.truevalue
        return self.falsevalue

    @property
    def docstring(self) -> str:
        docstring = self.generic_get_docstring()
        bool_values = [self.truevalue, self.falsevalue]
        bool_values = [v for v in bool_values if v != '']
        bool_values.sort()
        bool_str = ', '.join(bool_values)
        return f'{docstring}. possible values: {bool_str}'
    
    @property
    def optional(self) -> bool:
        return True

    @property
    def non_null_values(self) -> list[str]:
        null_values = ['', 'none', 'None', 'null', None]
        values = [self.truevalue, self.falsevalue]
        values = [v for v in values if v not in null_values]
        return values


@dataclass
class XMLSelectOption:
    value: str
    selected: bool
    ui_text: str

class XMLSelectParam(XMLInputParam):
    def __init__(self, name: str):
        super().__init__(name)
        self.options: list[XMLSelectOption] = []
        self.multiple: bool = False
        self.formats: list[str] = []

    @property
    def default(self) -> Any:
        for option in self.options:
            if option.selected:
                return option.value
        if len(self.options) > 0:
            return self.options[0].value
        return ''        

    @property
    def docstring(self) -> str:
        docstring = self.generic_get_docstring()
        option_values = [v.value for v in self.options]
        option_values.sort()
        option_str = ', '.join(option_values[:5])
        return f'{docstring}. possible values: {option_str}'
    
    @property
    def optional(self) -> bool:
        if self.multiple or self._optional:
            return True
        return False
    
    @property
    def array(self) -> bool:
        if self.multiple:
            return True
        return False

    @property
    def non_null_values(self) -> list[str]:
        null_values = ['', 'none', 'None', 'null', None]
        values = [opt.value for opt in self.options]
        values = [v for v in values if v not in null_values]
        return values
    

class XMLDataParam(XMLInputParam):
    def __init__(self, name: str):
        super().__init__(name)
        self.formats: list[str] = []
        self.multiple: bool = False
    
    @property
    def array(self) -> bool:
        return self.multiple


class XMLDataCollectionParam(XMLInputParam):
    def __init__(self, name: str, collection_type: str):
        super().__init__(name)
        self.collection_type = collection_type
        self.formats: list[str] = []
    
    @property
    def array(self) -> bool:
        return True

