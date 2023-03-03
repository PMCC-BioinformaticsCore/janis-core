

from typing import Any

import re
import janis_core as j


function_token_matcher = re.compile("^\$\{\s*?return\s+?(.+?);*?\s*?\}$")

single_token_matcher = re.compile("^\$\((.+)\)$")  
# "$(some value here)"

inline_expression_matcher = re.compile("\$\((.+?)\)")  
# valueFrom: "Hello, my name is $(name)

input_selector_matcher = re.compile("^inputs\.([A-z0-9_]+)$")  
# valueFrom: $(inputs.myName) -> InputSelector("myName")

string_matcher = re.compile('^".+?"$')  
# literal strings

int_matcher = re.compile("^\s*\d+\s*$")

float_matcher = re.compile("^\s*(\d*\.\d+)|(\d+\.\d*)\s*$")



def parse_number_from_string(text: str) -> float | int:
    """
    Parse a string that is expected to contain a number.
    :param text: str. the number in string.
    :return: float or int. Parsed num.
    """
    # if not isinstance(text, str):  # optional - check type
    #     raise TypeError("num should be a str. Got {}.".format(type(text)))
    if int_matcher.search(text):
        return int(text)
    if float_matcher.search(text):
        return float(text)
    raise ValueError("num is not a number. Got {}.".format(text))  # optional


def parse_basic_expression(expr: Any) -> Any:
    if expr is None:
        return None
    if not isinstance(expr, str):
        return expr
    match = single_token_matcher.match(expr)
    if match:
        return convert_javascript_token(match.groups()[0])

    bigger_match = function_token_matcher.match(expr)
    if bigger_match:
        return convert_javascript_token(bigger_match.groups()[0])

    tokens = set(inline_expression_matcher.findall(expr))

    string_format = f"{expr}"
    token_replacers = {}

    for token, idx in zip(tokens, range(len(tokens))):
        key = f"JANIS_CWL_TOKEN_{idx+1}"
        string_format = string_format.replace(f"$({token})", f"{{{key}}}")
        token_replacers[key] = convert_javascript_token(token)

    if len(token_replacers) == 0:
        return string_format

    return j.StringFormatter(string_format, **token_replacers)


def convert_javascript_token(token: str) -> Any:
    input_selector_match = input_selector_matcher.match(token)
    if input_selector_match:
        return j.InputSelector(input_selector_match.groups()[0])

    if token.endswith(".size"):
        return j.FileSizeOperator(convert_javascript_token(token[:-5]))
    if token.endswith(".basename"):
        return j.BasenameOperator(convert_javascript_token(token[:-9]))
    if token.endswith(".path"):
        # Ignore it because Janis will automatically put this back in where relevant
        return convert_javascript_token(token[:-5])
    if token.endswith(".contents"):
        return j.ReadContents(convert_javascript_token(token[:-9]))

    is_string = string_matcher.match(token)
    if is_string:
        return token[1:-1]

    try:
        return parse_number_from_string(token)
    except ValueError:
        pass

    j.Logger.warn(
        f"Couldn't translate javascript token, will use the placeholder '<expr>{token}</expr>'"
    )
    return f"<expr>{token}</expr>"