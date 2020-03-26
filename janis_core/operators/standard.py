from janis_core.types.common_data_types import String
from janis_core.operators.operator import FunctionOperator


class BasenameOperator(FunctionOperator):
    def to_wdl(self, unwrap_operator, *args):
        arg = args[0]
        return f"basename({unwrap_operator(arg)})"

    def to_cwl(self, unwrap_operator, *args):
        return unwrap_operator(args[0]) + ".basename"

    def argtypes(self):
        return [String]

    def returntype(self):
        return String
