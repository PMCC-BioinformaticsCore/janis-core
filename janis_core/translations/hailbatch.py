"""
Not implemented yet:

- present_as
- secondaries_present_as

"""


import inspect
from textwrap import indent
from typing import Tuple, Union, Callable

from janis_core import (
    Selector,
    ToolInput,
    ToolArgument,
    CodeTool,
    InputNodeSelector,
    StepOutputSelector,
    StringFormatter,
    InputSelector,
    FirstOperator,
    AliasSelector,
    ForEachSelector,
    WildcardSelector,
    Operator,
    apply_secondary_file_format_to_filename,
    TInput,
)
from janis_core.tool.commandtool import CommandTool
from janis_core.translations.janis import JanisTranslator
from janis_core.types.common_data_types import *
from janis_core.workflow.workflow import WorkflowBase

SED_REMOVE_EXTENSION = "| sed 's/\\.[^.]*$//'"
REMOVE_EXTENSION = (
    lambda x, iterations: f"$(echo '{x}' {iterations * SED_REMOVE_EXTENSION})"
    if iterations > 0
    else x
)


class HailBatchTranslator:
    @classmethod
    def get_method_name_for_id(cls, identifier):
        return f"add_{identifier}_step"

    @classmethod
    def janis_type_to_py_annotation(cls, dt: DataType, skip_typing=False):
        annotation = None
        if isinstance(dt, Array):
            inner = cls.janis_type_to_py_annotation(dt.subtype())
            annotation = f"List[{inner}]"
        elif isinstance(dt, UnionType):
            inner = set(cls.janis_type_to_py_annotation(t) for t in dt.subtypes)
            if len(inner) == 1 or skip_typing:
                annotation = list(inner)[0]
            else:
                annotation = f"Union[{', '.join(inner)}]"
        elif dt.is_base_type((File, String, Directory)):
            annotation = "str"
        elif dt.is_base_type(Int):
            annotation = "int"
        elif dt.is_base_type((Float, Double)):
            annotation = "float"
        elif dt.is_base_type(Boolean):
            annotation = "bool"

        if annotation is None:
            Logger.info(f"Couldn't generate python type annotation for {dt.name}")
        elif dt.optional and not skip_typing:
            annotation = f"Optional[{annotation}]"
        return annotation

    @classmethod
    def translate_workflow(
        cls,
        workflow: WorkflowBase,
        allow_empty_container=False,
        container_override: dict = None,
        generate_click_cli=True,
    ) -> str:

        kwargs_no_annotations, kwargs, kwargs_with_defaults = [], [], []
        step_definitions = []
        step_calls = []
        inputs_to_read = []
        additional_expressions = []
        extra_imports = []

        for inp in workflow.input_nodes.values():
            dt = inp.datatype
            kwarg, has_default, ad_expr = cls.get_kwarg_from_value(
                inp.id(), dt, inp.default
            )
            additional_expressions.extend(ad_expr)

            if not has_default:
                kwargs.append(kwarg)
                kwarg_no_annotation, _, _ = cls.get_kwarg_from_value(
                    inp.id(), dt, None, include_annotation=False
                )
                kwargs_no_annotations.append(kwarg_no_annotation)
            else:
                kwargs_with_defaults.append(kwarg)

            if dt.is_base_type(File) or (
                isinstance(dt, Array) and dt.fundamental_type().is_base_type(File)
            ):
                inputs_to_read.append(inp)

        for stp in workflow.step_nodes.values():
            if isinstance(stp.tool, WorkflowBase):
                raise NotImplementedError("Hail Batch doesn't support subworkflows yet")
            if isinstance(stp.tool, CodeTool):
                raise NotImplementedError("No support for code tool yet")

            connections = stp.sources
            foreach = stp.foreach
            if stp.scatter:
                if len(stp.scatter.fields) == 1:
                    foreach = stp.scatter.fields[0]
                    connections[foreach] = ForEachSelector()
                else:
                    Logger.warn(
                        "Batch doesn't support scattering by fields at the moment"
                    )

            step_definitions.append(cls.translate_tool_internal(stp.tool, stp.id()))

            inner_connections = {}
            for k, con in connections.items():
                src = con.source().source
                inner_connections[k] = cls.unwrap_expression(src, code_environment=True)

            inner_connection_str = ", ".join(
                f"{k}={v}" for k, v in inner_connections.items()
            )
            inner_call = (
                f"{cls.get_method_name_for_id(stp.id())}(b, {inner_connection_str})"
            )
            if foreach is not None:

                foreach_str = cls.unwrap_expression(foreach, code_environment=True)
                call = f"""
{stp.id()} = []
for idx in {foreach_str}:
    {stp.id()}.append({inner_call})"""
            else:
                call = f"{stp.id()} = {inner_call}"

            if stp.when is not None:
                when_str = cls.unwrap_expression(stp.when, code_environment=True)
                call = f"""\
{stp.id()} = None
if {when_str}:
{indent(call, 4 * ' ')}
"""
            step_calls.append(call)

        # for outp in workflow.output_nodes.values():
        #     step_calls.append("b.write_output({source}, {name})")

        pd = 4 * " "
        additional_preparation_expressions_str = "\n".join(
            indent(t, pd) for t in additional_expressions
        )
        inputs_to_read_str = "\n".join(
            f"{pd}{inp.id()} = "
            + cls.prepare_input_read_for_inp(inp.datatype, inp.id())
            for inp in inputs_to_read
        )
        step_calls_str = "\n".join(indent(s, pd) for s in step_calls)
        step_definitions_str = "\n\n".join(step_definitions)

        name_equals_main_arg = None
        if generate_click_cli:
            function_name = "main_from_click"
            extra_imports.append("import click")
            name_equals_main_arg = f"""\
{cls.generate_click_function(workflow.tool_inputs(), to_call="main", click_function_name=function_name)}
if __name__ == "__main__":
    {function_name}()
"""
        else:
            name_equals_main_arg = f"""\
if __name__ == "__main__":
    main({", ".join(k + "=None" for k in kwargs_no_annotations)})
"""
        extra_imports_str = "\n" + "\n".join(extra_imports)
        retval = f"""\
import math, os
from typing import Union, Optional, List
{extra_imports_str}

import hailtop.batch as hb

def main({', '.join([*kwargs, *kwargs_with_defaults])}):
    b = hb.Batch('{workflow.id()}')
{additional_preparation_expressions_str}
{inputs_to_read_str}

{step_calls_str}
    return b
    
{step_definitions_str}

{inspect.getsource(apply_secondary_file_format_to_filename)}

{name_equals_main_arg}
"""

        try:
            import black

            try:
                return black.format_str(retval, mode=black.FileMode(line_length=82))
            except black.InvalidInput as e:
                Logger.warn(f"Couldn't format python code due to Black error: {e}")
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
            exts = [a for a in [dt.extension, *(dt.alternate_extensions or [])] if a]
            base = f"{reference_var}" + "".join(f'.replace("{e}", "")' for e in exts)
            dsec = {}
            for sec in dt.secondary_files():
                sec_key = sec.replace("^", "").replace(".", "")
                if "^" in sec:
                    sec_without_pattern = sec.replace("^", "")
                    dsec[sec_key] = f'{base} + "{sec_without_pattern}"'
                else:
                    dsec[sec_key] = f'{reference_var} + "{sec}"'

            dsec_str = ", ".join(f"{k}={v}" for k, v in dsec.items())
            return f"{batch_var}.read_input_group({ext}={reference_var}, {dsec_str})"
        else:
            return f"{batch_var}.read_input({reference_var})"

    @classmethod
    def prepare_read_group_dictionary_from_dt(cls, datatype: DataType):
        if not isinstance(datatype, File):
            return None
        secs = datatype.secondary_files()
        if not secs:
            return None
        # if all extension are just additions, like ".vcf" and ".vcf.idx", it's:
        #   {"vcf": "{root}", "idx": "{root}.idx"}
        # else for example if it's ref.fasta and ref.dict, it's:
        #   {"fasta": "{root}.fasta", "dict": "{root}.dict"}

        extension = datatype.extension or ""
        extension_without_dot = extension.replace(".", "")
        nameroot_value = "{root}"

        # this only works if there's one hat
        d = {}
        if any(s.startswith("^") for s in secs):
            # nameroot_value = "{root}.fasta"
            nameroot_value = "{root}" + extension_without_dot

            if any(s.startswith("^^") for s in secs):
                Logger.warn(
                    f"Secondary file patterns in '{datatype.name}' ({secs}) with two carats (^^) are not supported in Batch, will "
                )

            for s in secs:
                sname = s.replace("^", "").replace(".", "")
                if "^" in s:
                    d[sname] = "{root}" + s.replace("^", "")
                else:
                    d[sname] = nameroot_value + s
        else:
            d = {s.replace(".", ""): nameroot_value + s for s in secs}

        d["root"] = nameroot_value

        return d

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
        additional_expressions = []
        pd = 4 * " "

        tinputs = tool.inputs()
        tinputs_map = {i.id(): i for i in tinputs}
        args_to_bind: List[ToolArgument] = [*(tool.arguments() or [])]
        inputs_specified_in_arguments = set()
        args_to_check = [tool.cpus({}), tool.memory({}), tool.disk({})]
        for arg in tool.arguments() or []:
            args_to_check.append(arg)
        filestocreate = tool.files_to_create() or {}
        for fn, file in filestocreate.items():
            args_to_check.extend([fn, file])

        while len(args_to_check) > 0:
            arg = args_to_check.pop(0)
            if arg is None:
                continue
            elif isinstance(arg, InputSelector):
                inputs_specified_in_arguments.add(arg.input_to_select)
            elif isinstance(arg, Operator):
                args_to_check.extend(arg.get_leaves())

        for inp in tinputs:
            is_required = not inp.input_type.optional
            is_specified = (
                inp.id() in tool.connections
                or inp.id() in inputs_specified_in_arguments
            )
            has_default = inp.default is not None
            is_filename = isinstance(inp.input_type, Filename)

            if not any([is_required, is_specified, has_default, is_filename]):
                continue

            if inp.position is not None or inp.prefix:
                args_to_bind.append(inp)

            if not is_filename or is_specified:
                kwarg, has_default, ad_expr = cls.get_kwarg_from_value(
                    inp.id(), inp.input_type, inp.default, include_annotation=False
                )
                additional_expressions.extend(ad_expr)
                (kwargs_with_defaults if has_default else kwargs).append(kwarg)

            if inp.input_type.is_base_type(File) and inp.localise_file:
                additional_expressions.append(f'j.command(f"mv {{{inp.id()}}} .")')
                # do same with secondary files
                secs = inp.input_type.secondary_files()
                if secs:
                    for sec in secs:
                        initial_ext, n_carats = cls.split_secondary_file_carats(sec)
                        src = f"{{{inp.id()}}}"
                        if n_carats:
                            src = REMOVE_EXTENSION(src, n_carats)
                        additional_expressions.append(
                            f'j.command(f"{src}{initial_ext} .")'
                        )

        # outputs
        command_extras = ""
        output_collectors = []
        for outp in tool.outputs():
            (
                out_command_extras,
                outputs_to_collect,
                out_additional_expressions,
            ) = cls.translate_tool_output(outp)
            if out_command_extras:
                command_extras += out_command_extras
            if outputs_to_collect:
                output_collectors.extend(outputs_to_collect)
            if out_additional_expressions:
                additional_expressions.extend(out_additional_expressions)

        def get_resolved_input_selector(inp):
            if not isinstance(inp, InputSelector):
                raise Exception(
                    f"Internal error when unwrapped input selector {inp} (type({type(inp)})"
                )
            default = inp.input_to_select
            if inp.input_to_select not in tinputs_map:
                Logger.warn(
                    f"Couldn't find input ({inp.input_to_select}) in tool {tool.id()}"
                )
                return default

            dt = tinputs_map[inp.input_to_select].input_type
            if not isinstance(dt, File) and not dt.secondary_files():
                return default

            return f"{default}.base"

        if (
            tool.base_command() == ["sh", "script.sh"]
            and filestocreate.get("script.sh") is not None
        ):
            # In the WDL converter, instead of breaking down the script, we
            # just write it as StringFormatter to the files_to_create["script.sh"]
            # we can unwrap it here manually I think.
            script_str_formatter = tool.files_to_create().get("script.sh")
            code_block = (
                cls.unwrap_expression(
                    script_str_formatter,
                    code_environment=False,
                    input_selector_overrider=get_resolved_input_selector,
                )
                .strip()
                .replace("\\n", "\\\\\n")
            )
            command_constructor_str = f'''\
    j.command(f"""{code_block}""")
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
                    ) = cls.get_command_argument_for_tool_input(
                        tool_input=arg, tool=tool
                    )
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

        # container
        kwargs_with_defaults.append(f'container="{tool.container()}"')
        additional_expressions.append("j.image(container)")

        # memory
        tmemory = tool.memory({})
        if tmemory is not None:
            tmemory = cls.unwrap_expression(tmemory, code_environment=False)
            additional_expressions.append(f"j.memory(f'{tmemory}G')")
        tdisk = tool.disk({})
        if tdisk is not None:
            tdisk = cls.unwrap_expression(tdisk, code_environment=False)
            additional_expressions.append(f"j.storage(f'{tdisk}G')")

        additional_expressions_str = "\n".join(
            indent(t, pd) for t in additional_expressions
        )
        output_collectors_str = "\n".join(
            f"{pd}j.command({o})" for o in output_collectors
        )

        return f"""\
def {cls.get_method_name_for_id(step_id or tool.id())}(b, {", ".join([*kwargs, *kwargs_with_defaults])}):
    j = b.new_job('{step_id}')
{additional_expressions_str}
    
{command_constructor_str}
{output_collectors_str}
    
    return j
"""

    @classmethod
    def translate_tool_output(cls, outp) -> Tuple[Optional[str], List[str], List[str]]:
        command_extras = None
        output_collectors, additional_expressions = [], []
        secs = (
            outp.output_type.secondary_files()
            if isinstance(outp.output_type, File)
            else None
        ) or []
        if secs:
            # prepare resource group
            sec = JanisTranslator.get_string_repr(
                cls.prepare_read_group_dictionary_from_dt(outp.output_type)
            )
            additional_expressions.append(
                f"j.declare_resource_group({outp.id()}={sec})"
            )

        values = []
        dests = [f"j.{outp.id()}"]
        sec_presents_as = outp.secondaries_present_as or {}
        if isinstance(outp.selector, Stdout):
            command_extras += f" > {{j.{outp.id()}}}"
        elif isinstance(outp.selector, Stderr):
            command_extras += f" 2> {{j.{outp.id()}}}"
        elif isinstance(outp.selector, WildcardSelector):
            gl = cls.unwrap_expression(outp.selector.wildcard, code_environment=True)
            values.append(gl)
            for s in secs:
                initial_sec = sec_presents_as.get(s, s).replace("^", "")
                final_ext, final_iters = cls.split_secondary_file_carats(s)
                values.append(gl + initial_sec)

                if final_iters:
                    dest_pattern = REMOVE_EXTENSION(f"j.{outp.id()}", final_iters)
                    dests.append(f'"{dest_pattern + final_ext}"')
                else:
                    dests.append(f"f'{{j.{outp.id()}}}{final_ext}")
        elif isinstance(outp.selector, Operator):
            # maybe check leaves for wildcard, because that won't be supported
            value = cls.unwrap_expression(outp.selector, code_environment=True)
            values = [value]
            for s in secs:
                initial_ext, initial_iters = cls.split_secondary_file_carats(
                    sec_presents_as.get(s, s)
                )
                final_ext, final_iters = cls.split_secondary_file_carats(s)
                if initial_iters:
                    escaped_value = f"{{{value}}}"
                    values.append(
                        f'f"{REMOVE_EXTENSION(escaped_value, initial_iters)}{initial_ext}"'
                    )
                else:
                    values.append(f"f'{{{value}}}{initial_ext}")

                if final_iters:
                    GET_BASE_OP = REMOVE_EXTENSION(f"{{j.{outp.id()}}}", final_iters)
                    dests.append(f'f"{GET_BASE_OP}{final_ext}"')
                else:
                    dests.append(f"f'{{j.{outp.id()}}}{final_ext}'")
        elif isinstance(outp.selector, Selector):
            value = cls.unwrap_expression(outp.selector, code_environment=True)
            values = [value]
            for s in secs:
                initial_ext, initial_iters = cls.split_secondary_file_carats(
                    sec_presents_as.get(s, s)
                )
                final_ext, final_iters = cls.split_secondary_file_carats(s)
                if initial_iters:
                    escaped_value = f"{{{value}}}"
                    values.append(
                        f'f"{REMOVE_EXTENSION(escaped_value, initial_iters)}{initial_ext}"'
                    )
                else:
                    values.append(f'f"{{{value}}}{initial_ext}"')

                if final_iters:
                    GET_BASE_OP = REMOVE_EXTENSION(f"{{j.{outp.id()}}}", final_iters)
                    dests.append(f'f"{GET_BASE_OP}{final_ext}"')
                else:
                    dests.append(f'f"{{j.{outp.id()}}}{final_ext}"')
        else:
            Logger.warn("Couldn't translate output selector ")

        if values and dests:
            for value, dest in zip(values, dests):
                output_collectors.append(
                    f"'ln \"{{value}}\" {{dest}}'.format(value={value}, dest={dest})"
                )

        return command_extras, output_collectors, additional_expressions

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
                expr = cls.unwrap_expression(inner_value, code_environment=False)
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

        requires_quotes = arg.shell_quote is not False

        # quote entire string
        q, sq = '"', '\\"'
        escape_quotes = lambda el: el.replace(q, sq)
        quote_value_if_required = (
            lambda el: f"{q}{escape_quotes(el)}{q}" if requires_quotes else el
        )

        value = cls.unwrap_expression(arg.value, code_environment=False).replace(
            "\\t", "\\\\t"
        )

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

        return f"command_args.append(f'{escape_quotes(retval)}')", []

    @classmethod
    def get_kwarg_from_value(
        cls, identifier: str, datatype: DataType, default: any, include_annotation=True
    ) -> Tuple[str, bool, List[str]]:
        has_default = default is not None or datatype.optional
        kwarg = identifier
        extra_statements = []

        if include_annotation:
            annotation = cls.janis_type_to_py_annotation(datatype)
            if annotation is not None:
                kwarg += f": {annotation}"

        if has_default:
            inner_default = None
            if isinstance(default, Selector):
                # can't
                if isinstance(default, StepOutputSelector) or (
                    isinstance(default, Operator)
                    and any(
                        isinstance(l, StepOutputSelector) for l in default.get_leaves()
                    )
                ):
                    Logger.warn(
                        f"Can't calculate default for identifier '{identifier}': '{default}' as it relies on "
                        f"the output of a step, and batch won't support this"
                    )
                else:
                    # insert extra statement to start of function to evaluate this value
                    extra_statements.append(
                        f"{identifier} = {identifier} if {identifier} is not None else {cls.unwrap_expression(default, code_environment=True)}"
                    )
            else:
                # useful to get_string_repr
                has_default = True
                inner_default = cls.unwrap_expression(default, code_environment=True)

            if has_default:
                kwarg += f"={inner_default}"

        return kwarg, has_default, extra_statements

    @classmethod
    def unwrap_expression(
        cls,
        value,
        code_environment=False,
        stringify_list=True,
        input_selector_overrider: Callable[[InputSelector], str] = None,
    ) -> Union[str, List[str]]:
        uwkwargs = {
            "input_selector_overrider": input_selector_overrider,
            "stringify_list": stringify_list,
        }

        if value is None:
            return "None"
        elif isinstance(value, list):
            values = [
                cls.unwrap_expression(e, code_environment=True, **uwkwargs)
                for e in value
            ]
            if stringify_list:
                return "[" + ", ".join(values) + "]"
            return values
        elif isinstance(value, (str, int, float)):
            if code_environment:
                return JanisTranslator.get_string_repr(value)
            else:
                return str(value)
        elif isinstance(value, InputNodeSelector):
            return value.id()
        elif isinstance(value, StepOutputSelector):
            return f"{value.node.id()}.{value.tag}"
        elif isinstance(value, InputSelector):
            val = None
            if input_selector_overrider is not None:
                val = input_selector_overrider(value)
            return val or value.input_to_select
        elif isinstance(value, Filename):
            prefix = value.prefix
            if not isinstance(prefix, str):
                Logger.warn(
                    "Hail batch does not support filenames generated from intermediate filenames"
                )
                prefix = "generated"
            return value.generated_filename({"prefix": prefix})
        elif isinstance(value, StringFormatter):
            retval = value.to_python(
                unwrap_operator=lambda val: cls.unwrap_expression(
                    val, code_environment=False, **uwkwargs
                )
            )
            if code_environment:
                return f'f"{retval}"'
            return retval
        elif isinstance(value, AliasSelector):
            return cls.unwrap_expression(
                value.inner_selector, code_environment=code_environment, **uwkwargs
            )
        elif isinstance(value, WildcardSelector):
            raise Exception("Wildcard selectors are not valid within operators")

        # i don't know about this, might need to be moved to somewhere later
        elif isinstance(value, Operator):
            val = cls.unwrap_operator(value, **uwkwargs)
            if not code_environment:
                val = "{" + str(val) + "}"
            return val
        elif isinstance(value, ForEachSelector):
            if code_environment:
                return "idx"
            return "{idx}"
        else:
            # return str(value)
            raise NotImplementedError(
                f"Can't unwrap value '{value}' of type {type(value)}"
            )

    @classmethod
    def unwrap_operator(cls, value: Operator, **kwargs):
        # assume code_environment is True
        inner_unwrap = lambda a: cls.unwrap_expression(
            a, code_environment=True, **kwargs
        )
        return value.to_python(inner_unwrap, *value.args)

    @staticmethod
    def split_secondary_file_carats(secondary_annotation: str):
        fixed_sec = secondary_annotation.lstrip("^")
        leading = len(secondary_annotation) - len(fixed_sec)
        return secondary_annotation[leading:], leading

    @classmethod
    def generate_click_function(
        cls,
        inputs: List[TInput],
        to_call="main",
        click_function_name="main_from_click",
        help=None,
    ):

        escape_string = lambda s: s.replace("\n", "\\n").replace('"', '\\"')
        help_if_relevant = f'help="{escape_string(help)}' if help else ""
        options = []
        for inp in inputs:
            inner_args = [f'"--{inp.id()}"', f'"{inp.id()}"']
            inner_annotation = inp.intype
            if isinstance(inp.intype, Boolean):
                inner_args.append("is_flag=True")
                inner_annotation = None
            elif isinstance(inp.intype, Array):
                # add list flags
                inner_args.append("multiple=True")
                inner_annotation = inp.intype.subtype()
            if inner_annotation is not None:
                intype = inp.intype
                if intype.is_array():
                    intype = intype.fundamental_type()
                annotation = cls.janis_type_to_py_annotation(intype, skip_typing=True)
                inner_args.append(f"type={annotation}")
            if not inp.intype.optional:
                inner_args.append("required=True")

            if inp.default is not None and not (
                isinstance(inp.default, Selector)
                or (
                    isinstance(inp.default, list)
                    and any(isinstance(a, Selector) for a in inp.default)
                )
            ):
                inner_args.append(
                    f"default={JanisTranslator.get_string_repr(inp.default)}"
                )
            if inp.doc and inp.doc.doc:
                safer = escape_string(inp.doc.doc)
                inner_args.append(f'help="{safer}"')

            options.append(f'@click.option({", ".join(inner_args)})')

        nl = "\n"
        return f"""\
@click.command({help_if_relevant})
{nl.join(options)}
def {click_function_name}(*args, **kwargs):
    return {to_call}(*args, **kwargs)
"""
