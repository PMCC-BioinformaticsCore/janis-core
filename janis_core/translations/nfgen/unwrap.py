
from typing import Any, Optional

from janis_core import Logger
from janis_core import (
    CommandTool, 
    ToolInput, 
    Operator, 
    Selector,
    AliasSelector, 
    InputSelector, 
    WildcardSelector, 
    StringFormatter
)
from janis_core.workflow.workflow import StepNode
from janis_core.types import (
    Filename,
    File,
    Directory
)


def unwrap_expression(
    value: Any,
    input_in_selectors: dict[str, Any],
    tool: Optional[CommandTool]=None,
    inputs_dict: Optional[dict[str, ToolInput]]=None,
    quote_string: bool=True,
    for_output: bool=False,
    skip_inputs_lookup: bool=False,
    in_shell_script: bool=False,
    var_indicator: Optional[str]=None,
    step_indicator: Optional[str]=None,
    **debugkwargs: Any,
) -> Any:
    """
    The main logic to unwrap a janis expression and represent it in Nextflow translation
    """
    assert(tool)

    if value is None:
        if quote_string:
            return "null"
        return None

    if isinstance(value, StepNode):
        raise Exception(
            f"The Step node '{value.id()}' was found when unwrapping an expression, "
            f"you might not have selected an output."
        )

    if isinstance(value, list):
        toolid = str(debugkwargs.get("tool_id", "unwrap_list_expression"))
        elements: list[Any] = []
        for i in range(len(value)):
            el = unwrap_expression(
                value=value[i],
                input_in_selectors=input_in_selectors,
                quote_string=quote_string,
                tool=tool,
                for_output=for_output,
                inputs_dict=inputs_dict,
                skip_inputs_lookup=skip_inputs_lookup,
                in_shell_script=in_shell_script,
                var_indicator=var_indicator,
                step_indicator=step_indicator,
                tool_id=toolid + "." + str(i),
            )
            elements.append(el)
        list_representation = f"[{', '.join(elements)}]"
        return list_representation

    elif isinstance(value, str):
        if quote_string:
            return f"'{value}'"
        return value

    elif isinstance(value, bool):
        if quote_string:
            return f"'{value}'"
        return value

    elif isinstance(value, int) or isinstance(value, float):
        return str(value)

    elif isinstance(value, Filename):
        formatted = value.generated_filename()
        return formatted

    elif isinstance(value, StringFormatter):
        return translate_string_formatter(
            selector=value,
            tool=tool,
            input_in_selectors=input_in_selectors,
            in_shell_script=in_shell_script,
            inputs_dict=inputs_dict,
            skip_inputs_lookup=skip_inputs_lookup,
            **debugkwargs,
        )

    elif isinstance(value, InputSelector):
        if for_output:
            el = prepare_filename_replacements_for(
                inp=value, 
                inputs_dict=inputs_dict, 
                input_in_selectors=input_in_selectors
            )
            return el
        return translate_input_selector(
            selector=value,
            inputs_dict=inputs_dict,
            tool=tool,
            input_in_selectors=input_in_selectors,
            skip_inputs_lookup=False,
            in_shell_script=in_shell_script,
        )

    elif isinstance(value, AliasSelector):
        return unwrap_expression(
            value=value.inner_selector,
            input_in_selectors=input_in_selectors,
            quote_string=quote_string,
            tool=tool,
            inputs_dict=inputs_dict,
            skip_inputs_lookup=skip_inputs_lookup,
            for_output=for_output,
            in_shell_script=in_shell_script,
            var_indicator=var_indicator,
            step_indicator=step_indicator,
        )

    elif isinstance(value, WildcardSelector):
        # raise Exception(
        #     f"A wildcard selector cannot be used as an argument value for '{debugkwargs}' {tool.id()}"
        # )
        return f"'{value.wildcard}'"
    elif isinstance(value, Operator):
        unwrap_expression_wrap = lambda x: unwrap_expression(
            value=x,
            input_in_selectors=input_in_selectors,
            quote_string=quote_string,
            tool=tool,
            for_output=for_output,
            inputs_dict=inputs_dict,
            skip_inputs_lookup=skip_inputs_lookup,
            in_shell_script=in_shell_script,
            **debugkwargs,
        )
        return value.to_nextflow(unwrap_expression_wrap, *value.args)

    elif callable(getattr(value, "nextflow", None)):
        if var_indicator is not None and step_indicator is not None:
            return value.nextflow(
                var_indicator=var_indicator, step_indicator=step_indicator
            )
        else:
            return value.nextflow()

    raise Exception(
        "Could not detect type %s to convert to input value" % type(value)
    )

def translate_string_formatter(
    selector: StringFormatter,
    tool: CommandTool,
    input_in_selectors: dict[str, Any],
    in_shell_script: bool=False,
    inputs_dict: Optional[dict[str, ToolInput]]=None,
    skip_inputs_lookup: bool=False,
    **debugkwargs: Any,
) -> str:
    """
    Translate Janis StringFormatter data type to Nextflow
    """
    if len(selector.kwargs) == 0:
        return str(selector)

    kwargreplacements = {
        k: f"{unwrap_expression(value=v, tool=tool, input_in_selectors=input_in_selectors, inputs_dict=inputs_dict, skip_inputs_lookup=skip_inputs_lookup, **debugkwargs)}"
        for k, v in selector.kwargs.items()
    }

    arg_val = selector._format
    for k in selector.kwargs:
        arg_val = arg_val.replace(f"{{{k}}}", f"${{{str(kwargreplacements[k])}}}")

    if in_shell_script:
        arg_val = arg_val.replace("\\", "\\\\")

    return arg_val

def prepare_filename_replacements_for(
    inp: Optional[Selector], 
    inputs_dict: Optional[dict[str, ToolInput]],
    input_in_selectors: dict[str, Any]
) -> Optional[str]:
    """
    Generate a string expression to represent a filename in Nextflow
    """
    if inp is None or not isinstance(inp, InputSelector):
        return None

    if not inputs_dict:
        return f"${inp.input_to_select}.name"
        # raise Exception(
        #     f"Couldn't generate filename as an internal error occurred (inputs_dict did not contain {inp.input_to_select})"
        # )

    if isinstance(inp, InputSelector):
        if inp.input_to_select not in inputs_dict:
            raise Exception(
                f"The InputSelector '{inp.input_to_select}' did not select a valid input"
            )

        tinp = inputs_dict.get(inp.input_to_select)
        intype = tinp.intype

        if intype.is_base_type((File, Directory)):
            potential_extensions = (
                intype.get_extensions() if intype.is_base_type(File) else None
            )

            base = f"{tinp.id()}"
            if intype.has_secondary_files():
                base = f"{tinp.id()}[0]"

            if inp.remove_file_extension and potential_extensions:
                base = f"{base}.simpleName"
            elif hasattr(tinp, "localise_file") and tinp.localise_file:
                base = f"{base}.name"

        elif isinstance(intype, Filename):
            base = str(
                unwrap_expression(
                    value=intype.generated_filename(),
                    input_in_selectors=input_in_selectors,
                    inputs_dict=inputs_dict,
                    for_output=True,
                )
            )
        
        elif (
            intype.is_array()
            and isinstance(intype.fundamental_type(), (File, Directory))
            and tinp.localise_file
        ):
            base = f"{tinp.id()}.map{{ el.name }}"
        
        else:
            base = f"{tinp.id()}"

        if intype.optional:
            default = "'generated'"
            if isinstance(intype, Filename):
                default = base
            replacement = f"({inp.input_to_select} && {inp.input_to_select} != 'None' && {inp.input_to_select} != '' ? {inp.input_to_select} : {default})"
        else:
            replacement = f"{base}"

        # return f"\"${{{replacement}}}\""
        return replacement

def translate_input_selector(
    selector: InputSelector,
    inputs_dict: Optional[dict[str, ToolInput]],
    tool: CommandTool,
    input_in_selectors: dict[str, Any],
    skip_inputs_lookup: bool=False,
    in_shell_script: bool=False,
):
    """
    Translate Janis InputSelector data type into Nextflow expressions
    """
    if tool.versioned_id() not in input_in_selectors:
        input_in_selectors[tool.versioned_id()] = set()

    input_in_selectors[tool.versioned_id()].add(selector.input_to_select)

    sel: str = selector.input_to_select
    if not sel:
        raise Exception(
            "No input was selected for input selector: " + str(selector)
        )

    skip_lookup = skip_inputs_lookup or sel.startswith("runtime_")

    if not skip_lookup:

        if inputs_dict is None:
            raise Exception(
                f"An internal error occurred when translating input selector '{sel}': the inputs dictionary was None"
            )
        if selector.input_to_select not in inputs_dict:
            raise Exception(
                f"Couldn't find the input '{sel}' for the InputSelector(\"{sel}\")"
            )

        tinp: ToolInput = inputs_dict[selector.input_to_select]

        intype = tinp.intype

        if intype.is_base_type((File, Directory)):
            if intype.has_secondary_files():
                sel = f"{sel}[0]"

        if selector.remove_file_extension:
            if intype.is_base_type((File, Directory)):

                potential_extensions = (
                    intype.get_extensions() if intype.is_base_type(File) else None
                )
                if selector.remove_file_extension and potential_extensions:
                    sel = f"{sel}.simpleName"

            elif intype.is_array() and isinstance(
                intype.fundamental_type(), (File, Directory)
            ):
                inner_type = intype.fundamental_type()
                extensions = (
                    inner_type.get_extensions()
                    if isinstance(inner_type, File)
                    else None
                )

                inner_sel = f"el.basename"
                if extensions:
                    for ext in extensions:
                        inner_sel += f'.replace(/{ext}$/, "")'
                sel = f"{sel}.map(function(el) {{ return {inner_sel}; }})"
            else:
                Logger.warn(
                    f"InputSelector {sel} is requesting to remove_file_extension but it has type {tinp.input_type.id()}"
                )

    if in_shell_script:
        sel = f"${{{sel}}}"

    return sel