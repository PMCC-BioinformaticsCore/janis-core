
from typing import Optional
from dataclasses import dataclass
from abc import ABC, abstractmethod
from janis_core.types import File, Directory, Array, DataType, Stdout
from janis_core import ToolOutput
from janis_core.operators.selectors import WildcardSelector, InputSelector
from janis_core.utils.secondary import apply_secondary_file_format_to_filename

from .. import utils


@dataclass
class ProcessOutput(ABC):
    name: str
    is_optional: bool

    @abstractmethod
    def get_string(self) -> str:
        ...
    
    @property
    def emit(self) -> str:
        return f', emit: {self.name}'
    
    @property
    def optional(self) -> str:
        if self.is_optional:
            return ', optional: true'
        return ''


@dataclass
class StdoutProcessOutput(ProcessOutput):

    def get_string(self) -> str:
        return f'stdout{self.optional}{self.emit}'


@dataclass
class ValProcessOutput(ProcessOutput):
    expression: str

    def get_string(self) -> str:
        return f'val {self.expression}{self.optional}{self.emit}'


@dataclass
class PathProcessOutput(ProcessOutput):
    expression: str

    def get_string(self) -> str:
        return f'path {self.expression}{self.optional}{self.emit}'


@dataclass
class TupleProcessOutput(ProcessOutput):
    qualifiers: list[str]
    expressions: list[str]

    @property
    def fields(self) -> str:
        out: str = ''
        for qual, expr in zip(self.qualifiers, self.expressions):
            out += f'{qual}({expr}), '
        out = out.rstrip(', ') # strip off the last comma & space
        return out
    
    def get_string(self) -> str:
        return f'tuple {self.fields}{self.optional}{self.emit}'

 


def create_outputs(out: ToolOutput) -> list[ProcessOutput]:
    datatype: DataType = out.output_type
    if isinstance(datatype, Array):
        return create_outputs_array(out)
    else:
        return create_outputs_single(out)

def create_outputs_array(out: ToolOutput) -> list[ProcessOutput]:
    basetype: Optional[DataType] = utils.get_base_type(out.output_type)
    assert(basetype)
    
    # secondaries array
    if isinstance(basetype, File) and basetype.has_secondary_files():
        # a path output per file type
        outputs: list[ProcessOutput] = []
        exts = utils.get_extensions(out.output_type, allow_symbols=True)
        for ext in exts:
            outputs.append(create_path_output_secondaries(out, ext))
        return outputs

    # file array
    if isinstance(basetype, (File, Directory)):
        return [create_path_output(out)]

    # nonfile array
    return [create_val_output(out)]

def create_outputs_single(out: ToolOutput) -> list[ProcessOutput]:
    basetype: Optional[DataType] = utils.get_base_type(out.output_type)
    assert(basetype)
    
    # stdout
    if isinstance(basetype, Stdout):
        return [create_stdout_output(out)]

    # file secondaries
    if isinstance(basetype, File) and basetype.has_secondary_files():
        outputs = [create_tuple_output_secondaries(out)]
        return outputs
    
    # file
    if isinstance(basetype, (File, Directory)):
        return [create_path_output(out)]
    
    # nonfile
    return [create_val_output(out)]



def create_stdout_output(out: ToolOutput) -> StdoutProcessOutput:
    # stdout
    optional = True if out.output_type.optional else False  # type: ignore
    new_output = StdoutProcessOutput(name=out.id(), is_optional=optional)
    return new_output

def create_path_output(out: ToolOutput) -> PathProcessOutput:
    # file
    optional = True if out.output_type.optional else False  # type: ignore
    expression = get_expression(out)
    new_output = PathProcessOutput(
        name=out.id(), 
        is_optional=optional, 
        expression=expression
    )
    return new_output

def create_val_output(out: ToolOutput) -> ValProcessOutput:
    # nonfile
    optional = True if out.output_type.optional else False  # type: ignore
    expression = get_expression(out)
    new_output = ValProcessOutput(
        name=out.id(), 
        is_optional=optional, 
        expression=expression
    )
    return new_output

def create_tuple_output_secondaries(out: ToolOutput) -> TupleProcessOutput:
    """
    secondaries
    eg BamBai:
        selector=WildcardSelector("*.bam"),
        secondaries_present_as={".bai": ".bai"},
    """
    assert(isinstance(out.output_type, File))
    optional = True if out.output_type.optional else False  # type: ignore
    qualifiers: list[str] = []
    expressions: list[str] = []
    
    primary_expr = get_expression(out)
    primary_expr_unquoted = primary_expr.strip('"')
    exts = utils.get_extensions(out.output_type, allow_symbols=True)
    for ext in exts:
        # primary file
        if ext not in out.secondaries_present_as:
            qual = 'path'
            expr = primary_expr
        # secondary file
        else:
            secondary_ext = out.secondaries_present_as[ext]
            secondary_expr: str = apply_secondary_file_format_to_filename(primary_expr_unquoted, secondary_ext)
            qual = 'path'
            expr = f'"{secondary_expr}"'
        qualifiers.append(qual)
        expressions.append(expr)

    new_output = TupleProcessOutput(
        name=out.id(), 
        is_optional=optional,
        qualifiers=qualifiers, 
        expressions=expressions
    )
    return new_output

def create_path_output_secondaries(out: ToolOutput, ext: str) -> PathProcessOutput:
    # array of secondaries
    print(new_output.get_string())
    raise NotImplementedError



def get_expression(out: ToolOutput) -> str:
    # WildcardSelector
    if out.selector is not None:
        if isinstance(out.selector, WildcardSelector):
            return get_expression_wildcard_selector(out.selector)
        elif isinstance(out.selector, InputSelector):
            return get_expression_input_selector(out.selector)    
    raise NotImplementedError

def get_expression_input_selector(selector: InputSelector) -> str:
    return selector.input_to_select

def get_expression_wildcard_selector(selector: WildcardSelector) -> str:
    if selector.select_first:
        raise NotImplementedError
    else:
        return f'"{selector.wildcard}"'
 


# def get_output_qualifier_for_outtype(
#     cls,
#     out_type: DataType,
# ) -> nfgen.OutputProcessQualifier:
#     """
#     Generate Nextflow output qualifier based on Janis output data type

#     :param out_type:
#     :type out_type:
#     :return:
#     :rtype:
#     """
#     # what ????
#     if isinstance(out_type, Array):
#         return nfgen.OutputProcessQualifier.tuple

#     # elif isinstance(out_type, Stdout):
#     #     return nfgen.OutputProcessQualifier.stdout

#     elif isinstance(out_type, (File, Directory, Stdout)):
#         return nfgen.OutputProcessQualifier.path

#     return nfgen.OutputProcessQualifier.val



# def gen_outputs_for_process(
#     cls, process_name: str, tool: Tool
# ) -> List[nfgen.ProcessOutput]:
#     """
#     Generate a list of tool outputs in the form of nfgen.ProcessOutput objects

#     :param process_name:
#     :type process_name:
#     :param tool:
#     :type tool:
#     :return:
#     :rtype:
#     """
#     outputs: List[ToolOutput] = tool.outputs()

#     nf_outputs = []
#     for o in outputs:
#         output_type = o.output_type
#         selector = o.selector

#         qual = cls.get_output_qualifier_for_outtype(output_type)
#         expression = cls.unwrap_expression(
#             selector, inputs_dict=tool.inputs_map(), tool=tool, for_output=True
#         )

#         if isinstance(output_type, Array):
#             if isinstance(output_type.subtype(), (File, Directory)):
#                 sub_qual = nfgen.OutputProcessQualifier.path
#             else:
#                 sub_qual = nfgen.OutputProcessQualifier.val

#             tuple_elements = expression.strip("][").split(",")
#             formatted_list = []
#             for expression in tuple_elements:
#                 sub_exp = nfgen.TupleElementForOutput(
#                     qualifier=sub_qual, expression=expression
#                 )
#                 formatted_list.append(sub_exp.get_string())

#             expression = ", ".join(formatted_list)
#         elif isinstance(output_type, Stdout):
#             expression = f'"{settings.TOOL_STDOUT_FILENAME}_{process_name}"'
#         elif isinstance(output_type, File) and output_type.has_secondary_files():
#             sub_qual = nfgen.OutputProcessQualifier.path
#             tuple_elements = [expression]

#             primary_ext = output_type.extension
#             secondary_ext = []

#             if o.secondaries_present_as is not None:
#                 secondaries_present_as = o.secondaries_present_as
#             else:
#                 secondaries_present_as = {}

#             for ext in output_type.secondary_files():
#                 if ext in secondaries_present_as:
#                     secondary_ext.append(secondaries_present_as[ext])
#                 else:
#                     secondary_ext.append(ext)

#             for ext in secondary_ext:
#                 replacement = primary_ext + ext
#                 if ext.startswith("^"):
#                     replacement = ext[1:]

#                 sec_exp = None
#                 if primary_ext in expression:
#                     sec_exp = expression.replace(primary_ext, replacement)
#                 elif ".name" in expression:
#                     sec_exp = expression.replace(
#                         ".name", f".baseName + '{replacement}'"
#                     )

#                 if sec_exp is not None:
#                     tuple_elements.append(sec_exp)

#             formatted_list = []
#             for sec_exp in tuple_elements:
#                 tuple_el = nfgen.TupleElementForOutput(
#                     qualifier=sub_qual, expression=sec_exp
#                 )
#                 formatted_list.append(tuple_el.get_string())

#             expression = ", ".join(formatted_list)
#             qual = nfgen.OutputProcessQualifier.tuple

#         out = nfgen.ProcessOutput(
#             qualifier=qual,
#             name=o.id(),
#             expression=expression,
#             # is_optional=output_type.optional, # disable this because nextflow doesn't allow workflow to point to optional output
#         )

#         nf_outputs.append(out)

#     return nf_outputs

# def gen_output_expression(cls, o: Union[TOutput, ToolOutput]):
#     """
#     Based on the Janis output type, we generate string expression to represent outputs in Nextflow Process.

#     :param tool:
#     :type tool:
#     :param o:
#     :type o:
#     :return:
#     :rtype:
#     """
#     if isinstance(o, TOutput):
#         output_type = o.outtype
#     elif isinstance(o, ToolOutput):
#         output_type = o.output_type
#     else:
#         raise Exception("Unknown output object")

#     if isinstance(output_type, File):
#         expression = f"'*{output_type.extension}'"
#         qual = nfgen.OutputProcessQualifier.path
#     elif isinstance(output_type, Array) and isinstance(output_type.subtype(), File):
#         expression = f"'*{output_type.subtype().extension}'"
#         qual = nfgen.OutputProcessQualifier.path
#     else:
#         qual = nfgen.OutputProcessQualifier.val
#         expression = f"file(\"$workDir/{settings.PYTHON_CODE_OUTPUT_FILENAME_PREFIX}{o.tag}\").text.replace('[', '').replace(']', '').split(', ')"

#     return qual, expression
