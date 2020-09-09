class ToolEvaluation:
    def __init__(
        self,
        can_translate_to_cwl: bool,
        can_translate_to_wdl: bool,
        has_friendly_name: bool,
        has_contributors: bool,
    ):
        self.can_translate_to_cwl = can_translate_to_cwl
        self.can_translate_to_wdl = can_translate_to_wdl
        self.has_friendly_name = has_friendly_name
        self.has_contributors = has_contributors
