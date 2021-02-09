import json
import re
from typing import Dict, Tuple, List, Optional

from janis_core import Logger
from janis_core.tool.commandtool import ToolArgument, ToolInput, Tool, ToolOutput
from janis_core.workflow.workflow import StepNode, InputNode, OutputNode, WorkflowBase

from janis_core.tool.tool import ToolType
from janis_core.translations import TranslatorBase

from janis_core.types import (
    InputSelector,
    WildcardSelector,
    CpuSelector,
    String,
    Selector,
    Directory,
    Stdout,
    Stderr,
    Array,
    Boolean,
    Filename,
    File,
)

from janis_core.operators import Operator, StringFormatter

class BashTranslator(TranslatorBase):
    def __init__(self):
        super().__init__(name="bash")

    SHEBANG = "#!/usr/bin/env bash"
    LIB_FILENAME = "lib"

    @classmethod
    def translate_workflow(
        cls,
        workflow,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> Tuple[any, Dict[str, any]]:

        inputsdict = workflow.inputs_map()
        toolinputs_dict = {k: ToolInput(k, v.intype) for k, v in inputsdict.items()}

        step_keys = list(workflow.step_nodes.keys())
        tool_scripts = {}
        tools = {}

        tool_scripts[cls.LIB_FILENAME] = cls.generate_generic_functions()
        for step_id in workflow.step_nodes:
            tool = workflow.step_nodes[step_id].tool
            tools[tool.versioned_id()] = tool
            input_file_prefix = f"{tool.versioned_id()}.input"
            tool_scripts[input_file_prefix] = cls.generate_wf_step_input_vars(tool, step_keys)
            tool_scripts[tool.versioned_id()] = cls.tool_script(tool, step_id, input_file_prefix)

        return cls.workflow_script(workflow, tools), tool_scripts


    @classmethod
    def workflow_script(cls, workflow: WorkflowBase, tools: Dict[str, Tool]):
        lib = cls.generate_generic_functions()
        esc = '\\"'
        outputs = cls.generate_wf_tool_outputs(workflow)
        command = ""
        for tool_id in tools:
            command += f"source $TOOLDIR/{cls.tool_filename(tool_id)}\n"

        return f"""{cls.SHEBANG}
export DIR=$(pwd)
export STDOUTPATH=$(pwd)/stdout
export STDERRPATH=$(pwd)/stderr
export TOOLDIR=$2

# Load functions
source $TOOLDIR/{cls.tool_filename(cls.LIB_FILENAME)}

# source twice so we do not need to worry about order of variables (no self referenced variables here)
source $1
source $1

{command}

outputs="{json.dumps(outputs).replace('"', esc)}"
outputs="${{outputs//STDOUT/$STDOUTPATH}}"
outputs="${{outputs//STDERR/$STDERRPATH}}"
export outputs="${{outputs//DIR/$DIR}}"
echo $outputs

# END
"""

    @classmethod
    def tool_script(cls, tool: Tool, step_id: str, input_file_prefix: str):
        command = cls.generate_tool_command(tool)
        output_vars = cls.generate_wf_step_output_vars(tool, step_id)
        input_filename = cls.tool_filename(input_file_prefix)

        return f"""{cls.SHEBANG}
# Inputs
source $TOOLDIR/{input_filename}

# Print out full command for debugging purpose
echo "{command}"

# Executing command
# Redirect stdout to a file
# Redirect stderr to another file and stderr
{command} > $DIR/{step_id}_stdout 2> >(tee -a $DIR/{step_id}_stderr >&2) 

#Outputs
{output_vars}
"""

    @classmethod
    def translate_tool_internal(
        cls,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        params_to_include = None
        if tool.connections:
            params_to_include = set(tool.connections.keys())

        doc = f"# {tool.id()} bash wrapper"
        meta = tool.bind_metadata() or tool.metadata
        if params_to_include:
            doc += "\n\tNB: this wrapper only contains a subset of the available parameters"
        if meta and meta.documentation:
            doc += "".join(
                "\n# " + l for l in meta.documentation.splitlines(keepends=False)
            )

        outputs = cls.generate_outputs(tool)
        esc = '\\"'

        command = cls.generate_tool_command(tool)
        lib = cls.generate_generic_functions()

        return f"""{cls.SHEBANG}
export DIR=$(pwd)
export STDOUTPATH=$(pwd)/stdout
export STDERRPATH=$(pwd)/stderr
export TOOLDIR=$2

{lib}

# source twice so we do not need to worry about order of variables
source $1
source $1

{doc}

# Print out full command for debugging purpose
echo "{command}"

# Executing command
# Redirect stdout to a file
# Redirect stderr to another file and stderr
{command} > $STDOUTPATH 2> >(tee -a $STDERRPATH >&2)

outputs="{json.dumps(outputs).replace('"', esc)}"
outputs="${{outputs//STDOUT/$STDOUTPATH}}"
outputs="${{outputs//STDERR/$STDERRPATH}}"
outputs="${{outputs//DIR/$DIR}}"
echo $outputs
"""

    @classmethod
    def generate_tool_command(cls, tool):
        args = []
        for a in tool.arguments() or []:
            if a.prefix is not None or a.position is not None:
                args.append(a)

        for a in tool.inputs() or []:
            if a.prefix is not None or a.position is not None:
                args.append(a)

        args = sorted(args, key=lambda a: (a.prefix is None))
        args = sorted(args, key=lambda a: (a.position or 0))

        params_to_include = None
        if tool.connections:
            params_to_include = set(tool.connections.keys())

        bc = tool.base_command()
        if bc is None:
            bc = []
        elif not isinstance(bc, list):
            bc = [bc]

        output_args = []
        inputsdict = tool.inputs_map()
        for a in args:
            if isinstance(a, ToolInput):
                # if params_to_include and a.id() not in params_to_include:
                #     # skip if we're limiting to specific commands
                #     continue
                arg = translate_command_input(tool_input=a, inputsdict=inputsdict)
                if not arg:
                    Logger.warn(f"Parameter {a.id()} was skipped")
                    continue
                output_args.append(arg)
            else:
                arg = cls.translate_command_argument(tool_arg=a, inputsdict=inputsdict)
                if not arg:
                    Logger.warn(f"Argument {a} was skipped")
                    continue
                output_args.append(arg)

        str_bc = " ".join(f"{c}" for c in bc)
        command = " \\\n".join([str_bc, *[a for a in output_args]])

        return command

    @classmethod
    def generate_wf_step_input_vars(cls, tool: Tool, step_keys: List[str]):
        provided_inputs = cls.generate_wf_tool_inputs(tool, step_keys)
        tool_inputs = cls.build_inputs_file(
            tool, merge_resources=True,
            additional_inputs=provided_inputs
        )

        input_vars = ""
        for key in tool_inputs:
            val = tool_inputs[key] if not None else ""

            if key.endswith("WithPrefix"):
                var_no_prefix = key.split("WithPrefix")[0]

                input_vars += f"""{cls.SHEBANG}

if [[ ! -z "${var_no_prefix}" ]]
then
    {key}="{val}"
else
    {key}=""
fi
"""
            else:
                input_vars += f"export {key}=\"{val}\"\n"

        return input_vars

    @classmethod
    def generate_wf_step_output_vars(cls, tool: Tool, step_id: str):
        tool_outputs = cls.generate_outputs(tool)

        output_vars = ""
        for key in tool_outputs:
            var_name = f"{step_id}{key}"
            output_val = tool_outputs[key]
            if isinstance(output_val, list):
                output_val = " ".join(output_val)

            output_vars += f"""
{var_name}="{output_val}"
{var_name}="${{{var_name}//STDOUT/$STDOUTPATH}}"
{var_name}="${{{var_name}//STDERR/$STDERRPATH}}"
export {var_name}="${{{var_name}//DIR/$DIR}}"
echo "OUTPUT {var_name}: ${var_name}"
"""

        return output_vars

    @classmethod
    def generate_generic_functions(cls):
        return """
first()
{   
    for var in $@
    do  
        if [[ ! -z $var ]];
        then
            echo $var;
            break;
        fi
    done

    exit;
}

list() { 
    echo "$1"; 
}

join() { 
    local IFS="$1"; 
    shift; 
    echo "$*"; 
}
"""

    @classmethod
    def generate_outputs(cls, tool: Tool):

        inputsdict = tool.inputs_map()

        outputs = {}
        for out in tool.outputs():
            if isinstance(out.output_type, Array):
                val = []
                for sel in out.selector:
                    val.append("DIR/" + cls.unwrap_expression(sel, inputs_dict=inputsdict))

            elif isinstance(out.output_type, Stdout):
                val = "STDOUT"
            elif isinstance(out.output_type, Stderr):
                val = "STDERR"
            elif isinstance(out.output_type, File):
                sel = out.selector
                if sel is None:
                    sel = out.glob

                val = "DIR/" + cls.unwrap_expression(sel, inputs_dict=inputsdict)

            else:
                val = "XXX"

            outputs[out.tag] = val

        return outputs

    @classmethod
    def generate_wf_tool_inputs(cls, tool: Tool, step_keys: List[str]):
        inputs = {}
        for key in tool.connections:
            val = str(tool.connections[key])

            if val == "None":
                val = None
            else:
                if "inputs." in val:
                    # e.g. replace inputs.fastq to $fastq
                    val = val.replace("inputs.", "$")

                # e.g. replace bwamem.out to $bwamemout
                for tool_id in step_keys:
                    keyword = f"{tool_id}."
                    if keyword in val:
                        val = val.replace(keyword, f"${tool_id}")

            inputs[key] = val

        return inputs


    @classmethod
    def generate_wf_tool_outputs(cls, wf: WorkflowBase):
        step_keys = wf.step_nodes.keys()

        outputs = {}
        for o in wf.output_nodes:
            val = str(wf.output_nodes[o].source)

            if "inputs." in val:
                # e.g. replace inputs.fastq to $fastq
                val = val.replace("inputs.", "$")

            # e.g. replace bwamem.out to $bwamemout
            for tool_id in step_keys:
                keyword = f"{tool_id}."
                if keyword in val:
                    val = val.replace(keyword, f"${tool_id}")

            outputs[o] = val

        return outputs

    @staticmethod
    def snake_to_camel_case(string: str):
        parts = string.split("_")
        return parts[0] + "".join(x.title() for x in parts[1:])

    @classmethod
    def unwrap_expression(
            cls,
            value,
            code_environment=True,
            selector_override=None,
            tool=None,
            for_output=False,
            inputs_dict=None,
            skip_inputs_lookup=False,
            **debugkwargs,
    ):
        if value is None:
            if code_environment:
                return "null"
            return None

        if isinstance(value, StepNode):
            raise Exception(
                f"The Step node '{value.id()}' was found when unwrapping an expression, "
                f"you might not have selected an output."
            )

        if isinstance(value, list):
            toolid = debugkwargs.get("tool_id", "unwrap_list_expression")
            inner = " ".join(
                cls.unwrap_expression(
                    value[i],
                    code_environment=True,
                    selector_override=selector_override,
                    tool=tool,
                    tool_id=toolid + "." + str(i),
                    inputs_dict=inputs_dict,
                    skip_inputs_lookup=skip_inputs_lookup
                )
                for i in range(len(value))
            )
            # return cls.wrap_in_codeblock_if_required(
            #     f"$(list \"{inner}\")"
            #     f"", is_code_environment=code_environment
            # )

            return f"$(list \"{inner}\")"

        if isinstance(value, str):
            return value
            # if not code_environment:
            #     return value
            # return cls.quote_values_if_code_environment(
            #     cls.prepare_escaped_string(value), code_environment
            # )
        elif isinstance(value, int) or isinstance(value, float):
            return str(value)
        elif isinstance(value, Filename):
            # value.generated_filenamecwl() if code_environment else f"$({value.generated_filenamecwl()})"
            return cls.quote_values_if_code_environment(
                value.generated_filename(), code_environment
            )
        # elif isinstance(value, AliasSelector):
        #     return cls.unwrap_expression(
        #         value.inner_selector,
        #         code_environment=code_environment,
        #         selector_override=selector_override,
        #         inputs_dict=inputs_dict,
        #         for_output=for_output,
        #         tool=tool,
        #         **debugkwargs,
        #     )
        #
        elif isinstance(value, StringFormatter):
            return cls.translate_string_formatter(
                value,
                selector_override=selector_override,
                code_environment=code_environment,
                tool=tool,
                inputs_dict=inputs_dict,
                skip_inputs_lookup=skip_inputs_lookup,
                **debugkwargs,
            )
        # elif isinstance(value, InputNodeSelector):
        #     return translate_input_selector(
        #         InputSelector(value.id()),
        #         code_environment=code_environment,
        #         selector_override=selector_override,
        #         inputs_dict=inputs_dict,
        #         skip_inputs_lookup=True,
        #     )
        # elif isinstance(value, StepOutputSelector):
        #     sel = f"{value.node.id()}/{value.tag}"
        #     if sel in selector_override:
        #         return selector_override[sel]
        #     raise Exception(
        #         "An internal error occurred when unwrapping an operator, found StepOutputSelector with no alias"
        #     )
        # elif isinstance(value, ResourceSelector):
        #     if not tool:
        #         raise Exception(
        #             f"Tool must be provided when unwrapping ResourceSelector: {type(value).__name__}"
        #         )
        #     operation = value.get_operation(tool, hints={})
        #     return cls.unwrap_expression(
        #         operation,
        #         code_environment=code_environment,
        #         tool=tool,
        #         inputs_dict=inputs_dict,
        #         **debugkwargs,
        #     )
        #
        # elif for_output and isinstance(value, (Stderr, Stdout)):
        #     # next few ones we rely on the globs being
        #     if isinstance(value, Stdout):
        #         return "self[0]"
        #     elif isinstance(value, Stderr):
        #         return "self[1]"

        elif isinstance(value, InputSelector):
            if for_output:
                el = cls.prepare_filename_replacements_for(value, inputsdict=inputs_dict)
                return cls.wrap_in_codeblock_if_required(
                    el, is_code_environment=code_environment
                )
            return cls.translate_input_selector(
                selector=value,
                code_environment=code_environment,
                selector_override=selector_override,
                inputs_dict=inputs_dict,
                skip_inputs_lookup=skip_inputs_lookup
            )
        elif isinstance(value, WildcardSelector):
            raise Exception(
                f"A wildcard selector cannot be used as an argument value for '{debugkwargs}'"
            )
        elif isinstance(value, Operator):
            unwrap_expression_wrap = lambda exp: cls.unwrap_expression(
                exp,
                code_environment=True,
                selector_override=selector_override,
                tool=tool,
                for_output=for_output,
                inputs_dict=inputs_dict,
                skip_inputs_lookup=skip_inputs_lookup,
                **debugkwargs,
            )
            return cls.wrap_in_codeblock_if_required(
                value.to_shell(unwrap_expression_wrap, *value.args),
                is_code_environment=code_environment,
            )
        elif callable(getattr(value, "cwl", None)):
            return value.cwl()
        # elif isinstance(value, Operator):

        raise Exception(
            "Could not detect type %s to convert to input value" % type(value)
        )

    @classmethod
    def translate_input_selector(cls,
            selector: InputSelector,
            code_environment,
            inputs_dict,
            selector_override=None,
            skip_inputs_lookup=False,
    ):
        # TODO: Consider grabbing "path" of File

        sel: str = selector.input_to_select
        if not sel:
            raise Exception("No input was selected for input selector: " + str(selector))

        skip_lookup = skip_inputs_lookup or sel.startswith("runtime_")

        if selector_override and sel in selector_override:
            sel = selector_override[sel]

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
            if selector.remove_file_extension:
                if intype.is_base_type((File, Directory)):
                    potential_extensions = (
                        intype.get_extensions() if intype.is_base_type(File) else None
                    )
                    if selector.remove_file_extension and potential_extensions:
                        # sel = f"{sel}.basename"
                        # for ext in potential_extensions:
                        #     sel += f'.replace(/{ext}$/, "")'
                        for ext in potential_extensions:
                            sel = f"{{{sel}%{ext}}}"

                        sel = f"(basename \"${sel}\")"

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
            # elif tinp.localise_file:
            #     if intype.is_base_type((File, Directory)):
            #         sel += ".basename"
            #     elif intype.is_array() and isinstance(
            #             intype.fundamental_type(), (File, Directory)
            #     ):
            #         sel = f"{sel}.map(function(el) {{ return el.basename; }})"


        sel = f"${sel}"
        return sel if code_environment else f"$({sel})"

    @classmethod
    def translate_string_formatter(
            cls,
            selector: StringFormatter,
            selector_override,
            tool,
            code_environment=True,
            inputs_dict=None,
            skip_inputs_lookup=False,
            **debugkwargs,
    ):
        if len(selector.kwargs) == 0:
            return str(selector)

        kwargreplacements = {
            k: f"{cls.unwrap_expression(v, selector_override=selector_override, code_environment=True, tool=tool, inputs_dict=inputs_dict, skip_inputs_lookup=skip_inputs_lookup, **debugkwargs)}"
            for k, v in selector.kwargs.items()
        }

        arg_val = selector._format
        for k in selector.kwargs:
            arg_val = arg_val.replace(f"{{{k}}}", f"{str(kwargreplacements[k])}")

        return arg_val

    @classmethod
    def prepare_filename_replacements_for(cls,
            inp: Optional[Selector], inputsdict: Optional[Dict[str, ToolInput]]
    ) -> Optional[str]:
        if inp is None or not isinstance(inp, InputSelector):
            return None

        if not inputsdict:
            return "inputs." + inp.input_to_select + ".basename"
            # raise Exception(
            #     f"Couldn't generate filename as an internal error occurred (inputsdict did not contain {inp.input_to_select})"
            # )

        if isinstance(inp, InputSelector):
            if inp.input_to_select not in inputsdict:
                raise Exception(
                    f"The InputSelector '{inp.input_to_select}' did not select a valid input"
                )

            tinp = inputsdict.get(inp.input_to_select)
            intype = tinp.intype

            if intype.is_base_type((File, Directory)):
                potential_extensions = (
                    intype.get_extensions() if intype.is_base_type(File) else None
                )
                if inp.remove_file_extension and potential_extensions:
                    base = f"inputs.{tinp.id()}.basename"
                    for ext in potential_extensions:
                        base += f'.replace(/{ext}$/, "")'
                elif tinp.localise_file:
                    base = f"inputs.{tinp.id()}.basename"
                else:
                    base = f"inputs.{tinp.id()}"
            elif (
                    intype.is_array()
                    and isinstance(intype.fundamental_type(), (File, Directory))
                    and tinp.localise_file
            ):
                base = f"inputs.{tinp.id()}.map(function(el) {{ return el.basename; }})"
            else:
                base = "inputs." + tinp.id()

            if intype.optional:
                replacement = f'inputs.{tinp.id()} ? {base} : "generated"'
            else:
                replacement = f"{base}"

            return replacement

    @classmethod
    def wrap_in_codeblock_if_required(cls, value, is_code_environment):
        return value if is_code_environment else f"$({value})"

    @classmethod
    def quote_values_if_code_environment(cls, value, is_code_environment):
        return f'"{value}"' if is_code_environment else value

    @classmethod
    def prepare_escaped_string(cls, value: str):
        return json.dumps(value)[1:-1]

    @classmethod
    def translate_code_tool_internal(
        cls,
        tool,
        with_docker=True,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        raise Exception("CodeTool is not currently supported in bash translation")

    @classmethod
    def build_inputs_file(
            cls,
            tool,
            recursive=False,
            merge_resources=False,
            hints=None,
            additional_inputs: Dict = None,
            max_cores=None,
            max_mem=None,
            max_duration=None,
    ) -> Dict[str, any]:

        ad = additional_inputs or {}
        values_provided_from_tool = {}
        if tool.type() == ToolType.Workflow:
            values_provided_from_tool = {
                i.id(): i.value or i.default
                for i in tool.input_nodes.values()
                if i.value or (i.default and not isinstance(i.default, Selector))
            }
            # TODO: add mapping from inputdict

        elif tool.type() == ToolType.CommandTool:
            values_provided_from_tool = {
                i.id(): i.default
                for i in tool.tool_inputs()
            }

        # inp = {
        #     i.id(): ad.get(i.id(), values_provided_from_tool.get(i.id()))
        #     for i in tool.tool_inputs()
        #     if i.default is not None
        #        or not i.intype.optional
        #        or i.id() in ad
        #        or i.id() in values_provided_from_tool
        # }

        # Build input variables and another copy of each input variable with its prefix attached
        inp = {}
        #TODO: workflow input
        if tool.type() == ToolType.Workflow:
            for i in tool.tool_inputs():
                if i.default is not None \
                        or not i.intype.optional \
                        or i.id() in ad \
                        or i.id() in values_provided_from_tool:
                    val = ad.get(i.id(), values_provided_from_tool.get(i.id()))

                    if val == "" or val is None:
                        val = []

                    if not isinstance(val, list):
                        val = [val]

                    inp[i.tag] = " ".join(str(v) for v in val)
                    inp[i.tag + "WithPrefix"] = inp[i.tag]

        elif tool.type() == ToolType.CommandTool:
            for i in tool.inputs():
                if i.default is not None \
                   or not i.input_type.optional \
                   or i.tag in ad \
                   or i.tag in values_provided_from_tool:

                    prefix = i.prefix if i.prefix else ""
                    tprefix = prefix

                    separate_value_from_prefix = i.separate_value_from_prefix is not False
                    if prefix and separate_value_from_prefix:
                        tprefix += " "

                    ad.get(i.tag)
                    values_provided_from_tool.get(i.tag)

                    inputsdict = tool.inputs_map()

                    if isinstance(i.input_type, Filename):
                        val = cls.unwrap_expression(i.input_type.generated_filename(), inputs_dict=inputsdict)
                    elif isinstance(i.input_type, Boolean):
                        val = ad.get(i.tag, values_provided_from_tool.get(i.tag)) or ""
                        if val == "True":
                            val = True
                        if val == "False":
                            val = False
                    else:
                        val = ad.get(i.tag, values_provided_from_tool.get(i.tag)) or ""

                    if val == "" or val is None:
                        val = []

                    if not isinstance(val, list):
                        val = [val]

                    inp[i.tag] = " ".join(str(v) for v in val)

                    if len(val) > 0 and (i.prefix or i.position):
                        for v in val:
                            if isinstance(v, bool):
                                inp[i.tag + "WithPrefix"] = tprefix
                            else:
                                inp[i.tag + "WithPrefix"] = " ".join(tprefix + str(v) for v in val)
                    else:
                        inp[i.tag + "WithPrefix"] = ""

        if merge_resources:
            for k, v in cls.build_resources_input(
                    tool,
                    hints,
                    max_cores=max_cores,
                    max_mem=max_mem,
                    max_duration=max_duration,
            ).items():
                inp[k] = ad.get(k, v)

        return inp

    @staticmethod
    def stringify_translated_workflow(wf):
        return wf

    @staticmethod
    def stringify_translated_tool(tool):
        return tool

    @staticmethod
    def stringify_translated_inputs(inputs):
        lines = []
        Logger.debug(f"inputs: {inputs}")
        for key in inputs:
            val = str(inputs[key])
            lines.append(f"{key}=\"{val}\"")

        return "\n".join(lines)


    @staticmethod
    def workflow_filename(workflow):
        return workflow.versioned_id() + ".sh"

    @staticmethod
    def tool_filename(tool):
        prefix = tool
        if isinstance(tool, Tool):
            prefix = tool.versioned_id()

        return prefix + ".sh"

    @staticmethod
    def inputs_filename(workflow):
        return workflow.id() + ".input.sh"

    @staticmethod
    def resources_filename(workflow):
        return workflow.id() + "-resources.json"

    @staticmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        return None


    @classmethod
    def stdout_output_name(cls, tool):
        stdout_outputs = []
        for o in tool.outputs():
            if isinstance(o.output_type, Stdout):
                stdout_outputs.append(o.tag)

        if not stdout_outputs:
            return None

        if len(stdout_outputs) != 1:
            raise Exception("There is more than out one output with type Stdout")

        return stdout_outputs[0]


    @classmethod
    def translate_command_argument(cls, tool_arg: ToolArgument, inputsdict=None, **debugkwargs):
        # make sure it has some essence of a command line binding, else we'll skip it
        if not (tool_arg.position is not None or tool_arg.prefix):
            return None

        separate_value_from_prefix = tool_arg.separate_value_from_prefix is not False
        prefix = tool_arg.prefix if tool_arg.prefix else ""
        tprefix = prefix

        if prefix and separate_value_from_prefix:
            tprefix += " "

        arg_val = cls.unwrap_expression(tool_arg.value, inputsdict=inputsdict, skip_inputs_lookup=True)

        if tool_arg.shell_quote is not False:
            arg_val = f"\"{arg_val}\""

        arg_val = f"{tprefix}{arg_val}" if tprefix else f"{arg_val}"

        return arg_val


def translate_command_input(tool_input: ToolInput, inputsdict=None, **debugkwargs):
    # make sure it has some essence of a command line binding, else we'll skip it
    if not (tool_input.position is not None or tool_input.prefix):
        return None

    # separate_value_from_prefix = tool_input.separate_value_from_prefix is not False
    # prefix = tool_input.prefix if tool_input.prefix else ""
    # tprefix = prefix
    #
    # if prefix and separate_value_from_prefix:
    #     tprefix += " "
    #
    # name = tool_input.id()
    #
    # if tool_input.shell_quote is not False:
    #     name = f"'{name}'"
    #
    # # Replace all {inputs.VAR} with $VAR
    # while "{inputs." in name:
    #     regex = r'(.*)\{inputs\.(\w+)\}(.*)'
    #     name = re.sub(regex, r'\1$\2\3', name)
    #
    # arg_val = f"{tprefix}{name}" if tprefix else f"{name}"

    name = tool_input.id()
    return f"${name}WithPrefix"


if __name__ == "__main__":
    from janis_unix.tools import Echo

    Echo().translate("shell")
