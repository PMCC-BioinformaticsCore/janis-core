


class ParamNotSupportedError(Exception):
    """Raise when the tool XML possesses a non-supported param"""

class TagNotSupportedError(Exception):
    """Raise when the tool XML possesses a non-supported XML tag"""

class AttributeNotSupportedError(Exception):
    """Raise when the tool XML possesses a non-supported XML tag"""

class DuplicateParamError(Exception):
    """Raise when two param names clash - another param with the same name was already known"""

class InputError(Exception):
    """Raise when user inputs are not complete / correct"""

