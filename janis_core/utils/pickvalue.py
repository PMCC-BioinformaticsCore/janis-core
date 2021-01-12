from enum import Enum


class PickValue(Enum):
    """
    Based on: https://github.com/common-workflow-language/cwl-v1.2/blob/conditionals/design-documents/conditionals-2019.md

    first_non_null: Picks the first value in a list that's not null
    single_non_null: Picks the ONLY element in a list that's non null (and crashes if there are multiple non-null)
    all_non_null: Filters out the null elements
    """

    first_non_null = "first_non_null"
    single_non_null = "single_non_null"
    all_non_null = "all_non_null"
