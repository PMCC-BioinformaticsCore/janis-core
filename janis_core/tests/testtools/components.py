


from typing import Optional

from janis_core import (
    ToolOutput,
    ToolInput,
    CommandTool,
)

from janis_core.types import (
    String,
    File,
    Array,
    Int,
    Boolean,
    Stdout
)



### TOOLS ### -------------------------------------------------------

class ComponentsMandatoryTestTool(CommandTool):
    def tool(self) -> str:
        return "ComponentsMandatoryTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("flag_true", Boolean, position=4, prefix="--flag-true", default=True),
            ToolInput("flag_false", Boolean, position=5, prefix="--flag-false", default=False),
            
            ToolInput("pos_basic", File, position=1),
            ToolInput("pos_default", Int, default=95, position=2),

            ToolInput("opt_basic", String, position=6, prefix="--opt-basic"),
            ToolInput("opt_default", Int, position=7, default=5, prefix="--opt-default"),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class ComponentsOptionalTestTool(CommandTool):
    def tool(self) -> str:
        return "ComponentsOptionalTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            ToolInput("flag_true", Boolean(optional=True), position=4, prefix="--flag-true", default=True),
            ToolInput("flag_false", Boolean(optional=True), position=5, prefix="--flag-false", default=False),
            
            ToolInput("pos_optional", String(optional=True), position=3),
            ToolInput("opt_optional", String(optional=True), position=8, prefix="--opt-optional"),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class ComponentsMandatoryArrayTestTool(CommandTool):
    def tool(self) -> str:
        return "ComponentsMandatoryArrayTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            # ToolInput("flag_true", Array(Boolean), position=4, prefix="--flag-true", default=True),
            # ToolInput("flag_false", Array(Boolean), position=5, prefix="--flag-false", default=False),
            
            ToolInput("pos_basic_arr", Array(File), position=1),
            ToolInput("pos_default_arr", Array(Int), default=[1, 2], position=2),

            ToolInput("opt_basic_arr", Array(String), position=3, prefix="--opt-basic"),
            ToolInput("opt_default_arr", Array(Int), position=4, default=[100, 200, 300], prefix="--opt-default"),
            
            ToolInput(
                "opt_basic_arr_prefixeach", 
                Array(String), 
                position=5, 
                prefix="--opt-basic-prefixeach", 
                prefix_applies_to_all_elements=True
            ),
            ToolInput(
                "opt_default_arr_prefixeach", 
                Array(String), 
                position=6, 
                default=['hi', 'there'], 
                prefix="--opt-default-prefixeach", 
                prefix_applies_to_all_elements=True
            ),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"


class ComponentsOptionalArrayTestTool(CommandTool):
    def tool(self) -> str:
        return "ComponentsOptionalArrayTestTool"

    def base_command(self) -> Optional[str | list[str]]:
        return "echo"

    def inputs(self) -> list[ToolInput]:
        return [
            # ToolInput("flag_true", Array(Boolean), position=4, prefix="--flag-true", default=True),
            # ToolInput("flag_false", Array(Boolean), position=5, prefix="--flag-false", default=False),
            
            ToolInput("pos_optional_arr", Array(File, optional=True), position=1),
            ToolInput("opt_optional_arr", Array(String, optional=True), position=2, prefix="--opt-optional-arr"),
            ToolInput(
                "opt_optional_arr_prefixeach", 
                Array(String, optional=True), 
                position=3, 
                prefix="--opt-optional-arr-prefixeach", 
                prefix_applies_to_all_elements=True
            ),
        ]

    def outputs(self):
        return [ToolOutput("out", Stdout)]

    def container(self) -> str:
        return "ubuntu:latest"

    def version(self) -> str:
        return "TEST"

