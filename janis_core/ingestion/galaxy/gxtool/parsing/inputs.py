

from typing import Any
from janis_core.ingestion.galaxy.runtime.exceptions import ParamNotSupportedError

from ..model import (
    XMLParam,
    XMLInputParam,
    XMLSelectOption,
    XMLTextParam,
    XMLIntegerParam,
    XMLFloatParam,
    XMLBoolParam,
    XMLSelectParam,
    XMLDataParam,
    XMLDataCollectionParam
)


def parse_input_param(gxparam: Any) -> XMLParam:
    factory = XMLInputParamFactory()
    return factory.produce(gxparam)

class XMLInputParamFactory:
    def produce(self, gxparam: Any) -> XMLParam:
        match gxparam.type: # type: ignore
            case 'text':
                param = self.parse_text_param(gxparam)
            case 'integer':
                param = self.parse_int_param(gxparam)
            case 'float':
                param = self.parse_float_param(gxparam)
            case 'boolean':
                param = self.parse_bool_param(gxparam)
            case 'select':
                param = self.parse_select_param(gxparam)
            case 'data':
                param = self.parse_data_param(gxparam)
            case 'data_collection':
                param = self.parse_data_collection_param(gxparam)
            case 'data_column':
                param = self.parse_int_param(gxparam)
            case 'color':
                param = self.parse_text_param(gxparam)
                #logging.color_param_ignored()
            case _:
                raise ParamNotSupportedError(f'unknown param type: {str(gxparam.type)}')
        
        param = self.map_common_fields(gxparam, param)
        return param 
        
    def map_common_fields(self, gxparam: Any, param: XMLInputParam) -> XMLInputParam:
        param.label = str(gxparam.label)
        param.helptext = str(gxparam.help)
        param.argument = gxparam.argument
        param.set_optionality(bool(gxparam.optional))
        return param

    def parse_text_param(self, gxparam: Any) -> XMLTextParam:
        param = XMLTextParam(str(gxparam.flat_name))
        param.value = gxparam.value
        return param

    def parse_int_param(self, gxparam: Any) -> XMLIntegerParam:
        param = XMLIntegerParam(str(gxparam.flat_name))
        if hasattr(gxparam, 'value'):
            param.value = gxparam.value
        if hasattr(gxparam, 'min'):
            param.min = gxparam.min
        if hasattr(gxparam, 'max'):
            param.max = gxparam.max
        return param

    def parse_float_param(self, gxparam: Any) -> XMLFloatParam:
        param = XMLFloatParam(str(gxparam.flat_name))
        param.value = gxparam.value
        param.min = gxparam.min
        param.max = gxparam.max
        return param

    def parse_bool_param(self, gxparam: Any) -> XMLBoolParam:
        param = XMLBoolParam(str(gxparam.flat_name))
        param.checked = bool(gxparam.checked)
        param.truevalue = str(gxparam.truevalue)
        param.falsevalue = str(gxparam.falsevalue)
        return param

    def parse_select_param(self, gxparam: Any) -> XMLSelectParam:
        # TODO this could be dynamic options!
        param = XMLSelectParam(str(gxparam.flat_name))
        param.multiple = bool(gxparam.multiple)
        if hasattr(gxparam, 'static_options') and len(gxparam.static_options) > 0:
            for opt in gxparam.static_options:
                option = XMLSelectOption(value=opt[1], selected=opt[2], ui_text=opt[0])
                param.options.append(option)
        return param

    def parse_data_param(self, gxparam: Any) -> XMLDataParam:
        param = XMLDataParam(str(gxparam.flat_name))
        param.formats = gxparam.extensions
        param.multiple = bool(gxparam.multiple)
        return param

    def parse_data_collection_param(self, gxparam: Any) -> XMLDataCollectionParam:
        name = str(gxparam.flat_name)
        collection_type = gxparam.collection_types[0]
        param = XMLDataCollectionParam(name, collection_type)
        param.formats = gxparam.extensions
        return param

