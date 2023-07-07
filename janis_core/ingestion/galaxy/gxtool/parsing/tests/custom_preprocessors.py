

import os
from typing import Any


def custom_ftype_preprocessor(output_value: Any, **kwargs: dict[str, Any]):
    return os.path.splitext(output_value)

