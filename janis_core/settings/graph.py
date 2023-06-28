


# step scatter
ALLOW_NON_ARRAY_SCATTER_INPUT: bool = True  # whether to allow a non-array source to feed a scattered step input
ALLOW_UNKNOWN_SCATTER_FIELDS: bool = True  # whether to allow scatter for fields which are not tool inputs

# StepTagInput
ALLOW_INCORRECT_NUMBER_OF_SOURCES: bool = True  # whether to allow the number of sources on a StepTagInput to be incorrect for that tool Input

# Edge
ALLOW_INCOMPATIBLE_TYPES: bool = True # whether to allow data links between incompatibile types: eg String -> File
ALLOW_UNKNOWN_SOURCE: bool = True  # whether to allow edges where the source isn't registered on the workflow
