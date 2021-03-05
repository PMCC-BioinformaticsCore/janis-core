from textwrap import indent
from typing import Tuple, Dict, Union

from janis_core import (
    Selector,
    ToolInput,
    ToolArgument,
    CodeTool,
    InputNodeSelector,
    StepOutputSelector,
    StringFormatter,
    InputSelector,
    TwoValueOperator,
    IsDefined,
    If,
    JoinOperator,
    FirstOperator, BasenameOperator, AliasSelector, ForEachSelector, FilterNullOperator, WildcardSelector,
    IndexOperator, RangeOperator, LengthOperator,
)
from janis_core.translations.janis import JanisTranslator
from janis_core.types.common_data_types import *
from janis_core.tool.commandtool import CommandTool
from janis_core.workflow.workflow import WorkflowBase, InputNode
from janis_core.translations import TranslatorBase

# TO avoid format errors when unwrapping StringFormatter, we'll use the following dictionary
# that returns the key if it's missing:
class DefaultDictionary(dict):
    def __missing__(self, key):
        return key

class HailBatchTranslator:
    @classmethod
    def get_method_name_for_id(cls, identifier):
        return f"add_{identifier}_step"

    @classmethod
    def translate_workflow(
        cls,
        workflow: WorkflowBase,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> str:

        kwargs, kwargs_with_defaults = [], []
        step_definitions = []
        step_calls = []
        inputs_to_read = []

        for inp in workflow.input_nodes.values():
            dt = inp.datatype
            kwarg, has_default = cls.get_kwarg_from_value(inp.id(), dt, inp.default)
            (kwargs_with_defaults if has_default else kwargs).append(kwarg)
            if dt.is_base_type(File) or (
                isinstance(dt, Array) and dt.fundamental_type().is_base_type(File)
            ):
                inputs_to_read.append(inp)

        for stp in workflow.step_nodes.values():
            if isinstance(stp.tool, WorkflowBase):
                raise NotImplementedError("Hail Batch doesn't support subworkflows yet")
            if isinstance(stp.tool, CodeTool):
                raise NotImplementedError("No support for code tool yet")

            connections = stp.tool.connections
            foreach = stp.foreach
            if stp.scatter:
                if len(stp.scatter.fields) == 1:
                    foreach = stp.scatter.fields[0]
                    connections[foreach] = ForEachSelector()
                else:
                    Logger.warn("Batch doesn't support scattering by fields at the moment")

            step_definitions.append(cls.translate_tool_internal(stp.tool, stp.id()))

            inner_connections = {}
            for k, con in stp.tool.connections.items():
                inner_connections[k] = cls.unwrap_expression(con)

            inner_connection_str = ", ".join(
                f"{k}={v}" for k, v in inner_connections.items()
            )
            inner_call = f"{cls.get_method_name_for_id(stp.id())}(b, {inner_connection_str})"
            if foreach is not None:

                foreach_str = cls.unwrap_expression(foreach)
                call = f"""
{stp.id()} = []
for idx in {foreach_str}:
    {stp.id()}.append({inner_call})"""
            else:
                call = f"{stp.id()} = {inner_call}"
            step_calls.append(call)

        # for outp in workflow.output_nodes.values():
        #     step_calls.append("b.write_output({source}, {name})")

        pd = 4 * " "
        inputs_to_read_str = "\n".join(
            f"{pd}{inp.id()} = "
            + cls.prepare_input_read_for_inp(inp.datatype, inp.id())
            for inp in inputs_to_read
        )
        step_calls_str = "\n".join(indent(s, pd) for s in step_calls)
        step_definitions_str = "\n\n".join(step_definitions)
        retval = f"""\
import hailtop.batch as hb

def main({', '.join([*kwargs, *kwargs_with_defaults])}):
    b = hb.Batch('{workflow.id()}')

{inputs_to_read_str}

{step_calls_str}
    return b
    
{step_definitions_str}

if __name__ == "__main__":
    main({", ".join(k + "=None" for k in kwargs)})
"""

        try:
            import black

            try:
                return black.format_str(retval, mode=black.FileMode(line_length=82))
            except black.InvalidInput:
                Logger.warn(
                    "Check the generated Janis code carefully, as there might be a syntax error. You should report this error along with the workflow you're trying to generate from"
                )
        except ImportError:
            Logger.debug(
                "Janis can automatically format generated Janis code if you install black: https://github.com/psf/black"
            )

        return retval

    @classmethod
    def prepare_input_read_for_inp(cls, dt: DataType, reference_var, batch_var="b"):
        if isinstance(dt, Array):
            inner_ref = f"inner_{reference_var}"
            return f"[{cls.prepare_input_read_for_inp(dt._t, inner_ref)} for {inner_ref} in {reference_var}]"
        if not dt.is_base_type(File):
            return reference_var

        if isinstance(dt, File) and dt.secondary_files():
            # we need to build a reference group
            ext = (dt.extension or "root").replace(".", "")
            dsec = {}
            for sec in dt.secondary_files():
                sec_key = sec.replace("^", "").replace(".", "")
                if "^" in sec:
                    sec_without_pattern = sec.replace("^", "")
                    dsec[
                        sec_key
                    ] = f'{reference_var}[:-{len(ext) + 1}] + "{sec_without_pattern}"'
                else:
                    dsec[sec_key] = f'{reference_var} + "{sec}"'

            dsec_str = ", ".join(f"{k}={v}" for k, v in dsec.items())
            return f"{batch_var}.read_input_group({ext}={reference_var}, {dsec_str})"
        else:
            return f"{batch_var}.read_input({reference_var})"

    @classmethod
    def translate_tool_internal(
        cls,
        tool: CommandTool,
        step_id: str,
        allow_empty_container=False,
        container_override: dict = None,
    ):

        # we're only going to bind on tool inputs that are REQUIRED or are in the connections, or have a default

        kwargs = []
        kwargs_with_defaults = []
        pd = 4 * " "

        args_to_bind: List[ToolArgument] = [*(tool.arguments() or [])]
        for inp in tool.inputs():
            is_required = not inp.input_type.optional
            is_specified = tool.connections and inp.id() in tool.connections
            has_default = inp.default is not None
            is_filename = isinstance(inp.input_type, Filename)

            if not any([is_required, is_specified, has_default, is_filename]):
                continue

            if inp.position is not None or inp.prefix:
                args_to_bind.append(inp)

            if not is_filename or is_specified:
                kwarg, has_default = cls.get_kwarg_from_value(
                    inp.id(), inp.input_type, inp.default
                )
                (kwargs_with_defaults if has_default else kwargs).append(kwarg)

        # outputs
        command_extras = ""
        output_collectors = []
        for outp in tool.outputs():
            if isinstance(outp.selector, Stdout):
                command_extras += f"> {{j.{outp.id()}}}"
            if isinstance(outp.selector, Stderr):
                command_extras += f"2> {{j.{outp.id()}}}"
            if isinstance(outp.selector, Selector):
                value = cls.unwrap_expression(outp.selector)
                output_collectors.append(f'ln "{{{value}}}" {{j.{outp.id()}}}')

        if tool.base_command() == ["sh", "script.sh"] and tool.files_to_create().get("script.sh") is not None:
            # In the WDL converter, instead of breaking down the script, we
            # just write it as StringFormatter to the files_to_create["script.sh"]
            # we can unwrap it here manually I think.
            script_str_formatter = tool.files_to_create().get("script.sh")
            command_constructor_str = f'''\
    j.command(f"""\
{cls.unwrap_expression(script_str_formatter).strip()}
""")
'''
        else:
            command = ""
            bc = tool.base_command()
            if bc is not None:
                if isinstance(bc, list):
                    command += " ".join(bc)
                else:
                    command += bc

            command_args = []
            precommand_args = []
            for arg in sorted(args_to_bind, key=lambda el: el.position or 0):
                commandline_param, other_required_statements = None, None
                if isinstance(arg, ToolInput):
                    (
                        commandline_param,
                        other_required_statements,
                    ) = cls.get_command_argument_for_tool_input(tool_input=arg, tool=tool)
                elif isinstance(arg, ToolArgument):
                    (
                        commandline_param,
                        other_required_statements,
                    ) = cls.get_command_argument_for_tool_argument(arg)
                else:
                    commandline_param = str(arg)

                precommand_args.extend(other_required_statements)
                if commandline_param:
                    command_args.append(commandline_param)
            command_args_str = "\n".join(
                indent(s, pd) for s in [*precommand_args, *command_args]
            )
            command_extras_str = ("\\\\\n  " + command_extras) if command_extras else ""

            command_constructor_str = f"""\
    command_args = []
{command_args_str}
    nl = " \\\\\\n  "
    command = f'''
{command} {{"".join(nl + a for a in command_args)}} {command_extras_str}
    '''
"""
        # command += "".join(f" \\\\\\n    {s}" for s in command_args)

        kwargs_with_defaults.append(f'container="{tool.container()}"')



        output_collectors_str = "\n".join(
            f'{pd}j.command(f\'{o}\')' for o in output_collectors
        )

        return f"""\
def {cls.get_method_name_for_id(step_id or tool.id())}(b, {", ".join([*kwargs, *kwargs_with_defaults])}):
    j = b.new_job('{tool.id()}')
    
    j.image(container)
{command_constructor_str}
{output_collectors_str}
    j.memory("4Gb")
    
    return j
"""

    @classmethod
    def get_command_argument_for_tool_input(
        cls, tool_input: ToolInput, tool: CommandTool
    ):
        # quote entire string
        q, sq = '"', '\\"'
        escape_quotes = lambda el: el.replace(q, sq)
        quote_value = lambda el: f"{q}{escape_quotes(el)}{q}"

        other_required_statements = []
        # will resolve to: "if {check_condition}: append ..."
        check_condition = f"{tool_input.id()} is not None"

        intype = tool_input.input_type
        is_flag = isinstance(intype, Boolean)

        separate_value_from_prefix = tool_input.separate_value_from_prefix is not False
        prefix = tool_input.prefix if tool_input.prefix else ""
        tprefix = prefix

        # prepare values array, then join it all in the end

        if prefix and separate_value_from_prefix and not is_flag:
            tprefix += " "

        code_value = tool_input.id()
        if is_flag:
            if not tool_input.prefix:
                Logger.warn(
                    f"Tool input '{tool_input.id()}' was a flag, but didn't have prefix: skipping"
                )
                return None, []

            check_condition = f"{tool_input.id()} is True"
            code_value = quote_value(tool_input.prefix)

        elif intype.is_array():
            separator = (
                tool_input.separator if tool_input.separator is not None else " "
            )
            # should_quote = (
            #     isinstance(intype.subtype(), (String, File, Directory))
            #     and tool_input.shell_quote is not False
            # )
            if prefix:
                if tool_input.prefix_applies_to_all_elements:
                    other_required_statements.append(
                        f'[e for el in {tool_input.id()} for e in ["{tool_input.prefix}", el] ]'
                    )
                code_value = (
                    quote_value(prefix)
                    + " + "
                    + quote_value(separator)
                    + f".join({tool_input.id()})"
                )
            else:
                code_value = f'"{separator}".join({tool_input.id()})'
        # elif requires_quoting:
        #     pass
        else:
            if isinstance(intype, Filename):
                is_specified = tool_input.id() in tool.connections
                inner_value = intype
                if is_specified:
                    inner_value = FirstOperator(InputSelector(tool_input.id(), intype))

                check_condition = None
                expr = cls.unwrap_expression(inner_value, quote_strings=False)
                other_required_statements.append(
                    f'{tool_input.id()} = f"{escape_quotes(expr)}"'
                )
            if prefix:
                sep = " "
                if tool_input.separate_value_from_prefix is False:
                    sep = ""

                code_value = f'f"{prefix}{sep}{{{tool_input.id()}}}"'
            else:
                # code_value = tool_input.id()
                pass

        if check_condition is not None:
            slash = "\\"
            retval = f"""\
if {check_condition}:
    command_args.append({code_value})\
"""
        else:
            retval = f"command_args.append({code_value})"

        return retval, other_required_statements
        # old logic


    @classmethod
    def get_command_argument_for_tool_argument(cls, arg: ToolArgument):
        if "--java" in str(arg):
            return None, []
        requires_quotes = arg.shell_quote is not False

        # quote entire string
        q, sq = '"', '\\"'
        escape_quotes = lambda el: el.replace(q, sq)
        quote_value_if_required = lambda el: f"{q}{escape_quotes(el)}{q}" if requires_quotes else el


        value = cls.unwrap_expression(arg.value, quote_strings=False).replace("\\t",  "\\\\t")

        # we're quoting early, so we don't need to do any quoting later
        if isinstance(value, list):
            value = " ".join(quote_value_if_required(v) for v in value)
        else:
            value = quote_value_if_required(value)

        if arg.prefix is None:
            retval = value
        else:
            sep = ""
            if arg.separate_value_from_prefix or arg.separate_value_from_prefix is None:
                sep = " "
            retval = sep.join([arg.prefix, value])

        return f'command_args.append(f\'{escape_quotes(retval)}\')', []

    @classmethod
    def get_kwarg_from_value(
        cls, identifier: str, datatype: DataType, default: any
    ) -> Tuple[str, bool]:
        has_default = default is not None
        kwarg = identifier
        if isinstance(datatype, (String, Float, Int, Double)):
            # kwarg += ": {inp.input_type.toPythonPrim()}"
            pass
        if has_default:
            inner_default = None
            if isinstance(default, Selector):
                # can't
                Logger.warn(
                    f"Can't calculate default for identifier '{identifier}': '{default}' as batch won't support this"
                )
                has_default = False
            else:
                # useful to get_string_repr
                from janis_core.translations.janis import JanisTranslator

                inner_default = JanisTranslator.get_string_repr(default)

            if has_default:
                kwarg += f"={inner_default}"

        return kwarg, has_default

    @classmethod
    def unwrap_expression(cls, value, quote_strings=True) -> Union[str, List[str]]:
        if value is None:
            return "None"
        elif isinstance(value, list):
            return [cls.unwrap_expression(e) for e in value]
        elif isinstance(value, (str, int, float)):
            if quote_strings:
                return JanisTranslator.get_string_repr(value)
            else:
                return str(value)
        elif isinstance(value, InputNodeSelector):
            return value.id()
        elif isinstance(value, StepOutputSelector):
            return f"{value.node.id()}.{value.tag}"
        elif isinstance(value, InputSelector):
            return value.input_to_select
        elif isinstance(value, Filename):
            prefix = value.prefix
            if not isinstance(prefix, str):
                Logger.warn(
                    "Hail batch does not support filenames generated from intermediate filenames"
                )
                prefix = "generated"
            return value.generated_filename({"prefix": prefix})
        elif isinstance(value, StringFormatter):
            retval = value._format.format_map(
                DefaultDictionary({
                    k: f"{{{cls.unwrap_expression(v)}}}"
                    for k, v in value.kwargs.items()
                })
            )
            return retval
        # i don't know about this, might need to be moved to somewhere later
        elif isinstance(value, TwoValueOperator):
            arg1, arg2 = [cls.unwrap_expression(a) for a in value.args]
            return f"({arg1} {value.symbol()} {arg2})"
        elif isinstance(value, IsDefined):
            return f"({cls.unwrap_expression(value.args[0])} is not None)"
        elif isinstance(value, If):
            condition, iftrue, iffalse = [cls.unwrap_expression(a) for a in value.args]
            return f"({iftrue} if {condition} else {iffalse})"
        elif isinstance(value, JoinOperator):
            iterable, sep = [cls.unwrap_expression(a) for a in value.args]
            return f"{sep}.join({iterable})"
        elif isinstance(value, FilterNullOperator):
            iterable = JanisTranslator.get_string_repr(cls.unwrap_expression(value.args[0]))
            return f"[a for a in {iterable} if a is not None]"
        elif isinstance(value, FirstOperator):
            iterable = JanisTranslator.get_string_repr(cls.unwrap_expression(value.args[0]))
            return f"[a for a in {iterable} if a is not None][0]"
        elif isinstance(value, AliasSelector):
            return cls.unwrap_expression(value.inner_selector)
        elif isinstance(value, BasenameOperator):
            val = JanisTranslator.get_string_repr(cls.unwrap_expression(value.args[0]))

            return f"os.path.basename({val})"
        elif isinstance(value, WildcardSelector):
            gl = cls.unwrap_expression(value.wildcard)
            return f"glob({gl})"
        elif isinstance(value, IndexOperator):
            iterable, idx = [cls.unwrap_expression(a) for a in value.args]
            return f"{iterable}[{idx}]"
        elif isinstance(value, RangeOperator):
            iterable = cls.unwrap_expression(value.args[0])
            return f"range({iterable})"
        elif isinstance(value, LengthOperator):
            iterable = cls.unwrap_expression(value.args[0])
            return f"len({iterable})"
        elif isinstance(value, ForEachSelector):
            return "idx"



        else:
            # return str(value)
            raise NotImplementedError(
                f"Can't unwrap value '{value}' of type {type(value)}"
            )
