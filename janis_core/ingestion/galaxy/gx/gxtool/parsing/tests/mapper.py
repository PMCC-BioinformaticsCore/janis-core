


import operator
from typing import Callable

from janis_core.ingestion.galaxy.runtime.exceptions import AttributeNotSupportedError
from janis_core.tool.test_classes import (
    TTestExpectedOutput,
    TTestPreprocessor
)

from .checks import ValidCheck
from . import custom_ops


def unsupported_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    raise AttributeNotSupportedError()

def not_implemented_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    raise NotImplementedError(f'no test strategy for {vcheck.ctype}')

def lines_diff_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    return TTestExpectedOutput(
        tag=f'{vcheck.outname}_{vcheck.ctype}',
        preprocessor=TTestPreprocessor.LinesDiff,
        operator=operator.eq,
        file_diff_source=vcheck.reffile,
        expected_value=int(vcheck.value),
    )

def has_size_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    if vcheck.reffile:
        return TTestExpectedOutput(
            tag=f'{vcheck.outname}_{vcheck.ctype}',
            preprocessor=TTestPreprocessor.FileSize,
            operator=operator.eq,
            expected_file=str(vcheck.reffile)
        )
    elif vcheck.value:
        return TTestExpectedOutput(
            tag=f'{vcheck.outname}_{vcheck.ctype}',
            preprocessor=TTestPreprocessor.FileSize,
            operator=operator.eq,
            expected_value=int(vcheck.value)
        )

def has_text_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    if vcheck.reffile:
        return TTestExpectedOutput(
            tag=f'{vcheck.outname}_{vcheck.ctype}',
            preprocessor=TTestPreprocessor.FileContent,
            operator=custom_ops.has_text,
            expected_file=str(vcheck.reffile)
        )
    elif vcheck.value:
        return TTestExpectedOutput(
            tag=f'{vcheck.outname}_{vcheck.ctype}',
            preprocessor=TTestPreprocessor.FileContent,
            operator=custom_ops.has_text,
            expected_value=str(vcheck.value['text'])
        )
    raise AttributeError('either value or reffile needs to be supplied to each ValidCheck')
    
def has_n_lines_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    return TTestExpectedOutput(
        tag=f'{vcheck.outname}_{vcheck.ctype}',
        preprocessor=TTestPreprocessor.LineCount,
        operator=operator.eq,
        expected_value=int(vcheck.value['n'])
    )

def has_n_columns_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    return TTestExpectedOutput(
        tag=f'{vcheck.outname}_{vcheck.ctype}',
        preprocessor=TTestPreprocessor.FileContent,
        operator=custom_ops.has_n_columns,
        expected_value=int(vcheck.value['n'])
    )

def has_text_matching_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    return TTestExpectedOutput(
        tag=f'{vcheck.outname}_{vcheck.ctype}',
        preprocessor=TTestPreprocessor.FileContent,
        operator=custom_ops.has_text_matching,
        expected_value=str(vcheck.value['expression'])
    )

def ftype_strategy(vcheck: ValidCheck) -> TTestExpectedOutput:
    # return TTestExpectedOutput(
    #     tag=f'{vcheck.outname}_{vcheck.ctype}',
    #     preprocessor=custom_ftype_preprocessor,
    #     operator=in_list_op,
    #     #operator=operator.eq,
    #     file_diff_source=vcheck.reffile,
    #     expected_value=str(vcheck.value),
    # )
    raise AttributeNotSupportedError('cannot check filetype')




strategy_map: dict[str, Callable[..., TTestExpectedOutput]] = {
    'lines_diff': lines_diff_strategy,
    'ftype': ftype_strategy,
    'has_size': has_size_strategy,
    'has_text': has_text_strategy,
    'has_n_lines': has_n_lines_strategy,
    'has_n_columns': has_n_columns_strategy,
    'has_text_matching': has_text_matching_strategy,

    # not Implemented
    'not_has_text': not_implemented_strategy,
    'has_line': not_implemented_strategy,
    'has_line_matching': not_implemented_strategy,

    # not v0.1
    're_match': unsupported_strategy,
    're_match_multiline': unsupported_strategy,
    'has_archive_member': unsupported_strategy,
    'is_valid_xml': unsupported_strategy,
    'xml_element': unsupported_strategy,
    'has_element_with_path': unsupported_strategy,
    'has_n_elements_with_path': unsupported_strategy,
    'element_text_matches': unsupported_strategy,
    'element_text_is': unsupported_strategy,
    'attribute_matches': unsupported_strategy,
    'attribute_is': unsupported_strategy,
    'element_text': unsupported_strategy,
    'has_h5_keys': unsupported_strategy,
    'has_h5_attribute': unsupported_strategy
}


def map_ttestout(vcheck: ValidCheck) -> TTestExpectedOutput:
    map_func = strategy_map[vcheck.ctype]
    return map_func(vcheck)