



ALLOW_INCORRECT_NUMBER_OF_SOURCES: bool = False  # whether to allow the number of sources on a StepTagInput to be incorrect for that tool Input
ALLOW_NON_ARRAY_SCATTER_INPUT: bool = False  # whether to allow a non-array source to feed a scattered step input
ALLOW_INCOMPATIBLE_TYPES: bool = True # whether to allow data links between incompatibile types: eg String -> File
ALLOW_UNKNOWN_SOURCE: bool = False  # whether to allow edges where the source isn't registered on the workflow
ALLOW_UNKNOWN_SCATTER_FIELDS: bool = False  # whether to allow scatter for fields which are not tool inputs
