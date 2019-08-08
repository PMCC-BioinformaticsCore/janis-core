from janis_core import Array, CommandTool, ToolInput, String, ToolOutput


class SingleTestTool(CommandTool):
    @staticmethod
    def tool():
        return "TestStepTool"

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [ToolInput("inputs", String())]

    def friendly_name(self):
        return None

    def outputs(self):
        return [ToolOutput("out", String())]

    @staticmethod
    def container():
        return None

    @staticmethod
    def version():
        return None


class ArrayTestTool(CommandTool):
    @staticmethod
    def tool():
        return "ArrayStepTool"

    def friendly_name(self):
        return None

    @staticmethod
    def base_command():
        return "echo"

    def inputs(self):
        return [ToolInput("inputs", Array(String()))]

    def outputs(self):
        return [ToolOutput("outs", Array(String()))]

    @staticmethod
    def container():
        return None

    @staticmethod
    def version():
        return None
