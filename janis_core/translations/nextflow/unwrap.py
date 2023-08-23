

from __future__ import annotations
# from typing import TYPE_CHECKING
# if TYPE_CHECKING:
#     from .parsing.process.VariableManager import VariableManager

from .variables import VariableManager
from .variables import VariableType
from .variables import Variable

from copy import deepcopy
from typing import Any, Optional, Type
NoneType = type(None)
import re

from janis_core import (
    CommandTool, 
    ToolInput, 
    Operator
)
from janis_core.graph.steptaginput import Edge, StepTagInput
from janis_core.operators.operator import IndexOperator
from janis_core.utils.scatter import ScatterMethod
from janis_core.workflow.workflow import StepNode, InputNode
from janis_core.types import (
    Filename,
    File,
    DataType
)

from janis_core.operators.logical import (
    IsDefined,
    If,
    AssertNotNull,
    FloorOperator,
    CeilOperator,
    RoundOperator
)
from janis_core.operators.standard import (
    ReadContents,
    ReadJsonOperator,
    JoinOperator,
    BasenameOperator,
    NamerootOperator,
    NameextOperator,
    TransposeOperator,
    LengthOperator,
    RangeOperator,
    FlattenOperator,
    ApplyPrefixOperator,
    FileSizeOperator,
    FirstOperator,
    FilterNullOperator,
    ReplaceOperator
)
from janis_core.operators.operator import (
    IndexOperator,
    AsStringOperator,
    AsBoolOperator,
    AsIntOperator,
    AsFloatOperator,
    TwoValueOperator
)
from janis_core.operators.selectors import (
    InputNodeSelector, 
    StepOutputSelector,
    AliasSelector, 
    InputSelector, 
    MemorySelector,
    CpuSelector,
    DiskSelector,
    TimeSelector,
    WildcardSelector, 
    Selector,
)
from janis_core.operators.stringformatter import StringFormatter
from janis_core import translation_utils as utils
from janis_core.translation_utils import DTypeType

from . import naming
from . import nfgen_utils

# from .plumbing import cartesian_cross_subname
from .expressions import stringformatter_matcher

"""
NOTE: 
only unwrapping simple references for step inputs. 

reason: 
InputNodeSelector(input_fasta).name
--/-> 
ch_input_fasta.name

nextflow wants to provide process inputs as channels.
cant do channel.name etc, doesn't make sense. 
this stuff is supposed to be done inside the process. 
"""

def unwrap_expression(
    val: Any,
    # scope: Scope,
    context: str='general',
    variable_manager: Optional[VariableManager]=None,
    tool: Optional[CommandTool]=None,
    
    scatter_target: bool=False,
    scatter_method: Optional[ScatterMethod]=None,

    in_shell_script: bool=False,
    quote_strings: Optional[bool]=None,
    ) -> Any:

    unwrapper = Unwrapper(
        # scope=scope,
        context=context,
        variable_manager=variable_manager,
        tool=tool,

        scatter_target=scatter_target,
        scatter_method=scatter_method,

        in_shell_script=in_shell_script,
        quote_strings=quote_strings,
    )
    return unwrapper.unwrap(val)



class Unwrapper:
    """
    The main logic to unwrap a janis expression and represent it in Nextflow translation
    """
    def __init__(
        self,
        # scope: Scope,
        context: str,
        variable_manager: Optional[VariableManager],
        tool: Optional[CommandTool],
         
        scatter_target: bool,
        scatter_method: Optional[ScatterMethod], 
        
        in_shell_script: bool, 
        quote_strings: Optional[bool],
    ) -> None:
        # self.scope = scope
        self.context = context
        self.vmanager = variable_manager
        self.tool = tool
        self.scatter_target = scatter_target
        self.scatter_method = scatter_method
        self.quote_strings = quote_strings
        self.in_shell_script = in_shell_script

        self.operator_stack: list[str] = []
        self.operator_stack_ignore: list[Type[Any]] = [
            WildcardSelector,
            AliasSelector,
            InputNodeSelector,
            StepOutputSelector,
            StringFormatter,
            MemorySelector,
            CpuSelector,
            DiskSelector,
            TimeSelector,
        ]

        self.func_switchboard = {
            # primitives
            NoneType: self.unwrap_null,
            str: self.unwrap_str,
            bool: self.unwrap_bool,
            int: self.unwrap_int,
            float: self.unwrap_float,
            list: self.unwrap_list,

            # logical operators
            IsDefined: self.unwrap_is_defined_operator,
            If: self.unwrap_if_operator,
            AssertNotNull: self.unwrap_assert_not_null_operator,
            FloorOperator: self.unwrap_floor_operator,
            CeilOperator: self.unwrap_ceil_operator,
            RoundOperator: self.unwrap_round_operator,

            # operator operators
            IndexOperator: self.unwrap_index_operator,
            AsStringOperator: self.unwrap_as_string_operator,
            AsBoolOperator: self.unwrap_as_bool_operator,
            AsIntOperator: self.unwrap_as_int_operator,
            AsFloatOperator: self.unwrap_as_float_operator,
            
            # standard operators
            ReadContents: self.unwrap_read_contents_operator,
            ReadJsonOperator: self.unwrap_read_json_operator,
            JoinOperator: self.unwrap_join_operator,
            BasenameOperator: self.unwrap_basename_operator,
            NamerootOperator: self.unwrap_nameroot_operator,
            NameextOperator: self.unwrap_nameext_operator,
            TransposeOperator: self.unwrap_transpose_operator,
            LengthOperator: self.unwrap_length_operator,
            RangeOperator: self.unwrap_range_operator,
            FlattenOperator: self.unwrap_flatten_operator,
            ApplyPrefixOperator: self.unwrap_apply_prefix_operator,
            FileSizeOperator: self.unwrap_file_size_operator,
            FirstOperator: self.unwrap_first_operator,
            FilterNullOperator: self.unwrap_filter_null_operator,
            ReplaceOperator: self.unwrap_replace_operator,
            
            # selectors
            AliasSelector: self.unwrap_alias_selector,
            InputNodeSelector: self.unwrap_input_node_selector,
            WildcardSelector: self.unwrap_wildcard_selector,
            InputSelector: self.unwrap_input_selector,
            StepOutputSelector: self.unwrap_step_output_selector,
            MemorySelector: self.unwrap_memory_selector,
            CpuSelector: self.unwrap_cpu_selector,
            DiskSelector: self.unwrap_disk_selector,
            TimeSelector: self.unwrap_time_selector,

            # misc 
            StepTagInput: self.unwrap_step_tag_input,
            Edge: self.unwrap_edge,
            StringFormatter: self.unwrap_string_formatter,
            Filename: self.unwrap_filename,
            InputNode: self.unwrap_input_node,
        }
    
    def unwrap(self, val: Any) -> Any:
        """main public method"""

        # add to operator stack if entity is a Selector
        self.update_operator_stack(val, life_cycle='start')

        # do the unwrap  
        expr = self.unwrap_entity(val)

        # apply quotes if needed
        if self.should_quote(val, expr):
            expr = f'"{expr}"'

        # apply nextflow curly braces if needed
        if self.should_apply_curly_braces(val):
            expr = f'${{{expr}}}'

        # pop from operator stack if entity is a Selector
        self.update_operator_stack(val, life_cycle='end')
        return expr

    def unwrap_entity(self, entity: Any) -> Any:
        """switchboard to delegate unwrap based on entity type."""

        # most entities: custom unwrap function
        etype = type(entity)
        if etype in self.func_switchboard:
            func = self.func_switchboard[etype]
            expr = func(entity)

        # TwoValueOperators: shared unwrap function
        elif isinstance(entity, TwoValueOperator):
            expr = self.unwrap_two_value_operator(entity)

        # anything else with a .to_nextflow() method
        elif callable(getattr(entity, "to_nextflow", None)):
            expr = self.unwrap_operator(entity)
        
        # errors ---
        elif isinstance(entity, StepNode):
            raise Exception(
                f"The Step node '{entity.id()}' was found when unwrapping an expression, "
                f"you might not have selected an output."
            )
        else:
            raise Exception(
                "Could not detect type %s to convert to input value" % type(entity)
            )

        return expr


    ### OPERATOR STACK 
    
    """
    --- operator_stack ---

    keeps track of how deep we are into Janis operators. 
    This is needed to properly enclose curly braces. For example: 
        
    Janis:      InputSelector("read", remove_file_extension=True) + "_fastqc.zip"
    ->
    Nextflow:   "${read.simpleName}_fastqc.zip"
    
    The InputSelector function should appear in curly braces, while the "_fastqc.zip" should not.
    We apply curly braces for the outermost context of janis operators. 
    eg Operator(Operator(Operator)) would apply curly braces only for the outermost Operator, not for each.
    """
    
    def update_operator_stack(self, val: Any, life_cycle: str='start') -> None:
        if isinstance(val, Selector) and type(val) not in self.operator_stack_ignore:
            # start of main unwrap() function
            if life_cycle == 'start':
                self.operator_stack.append(val.__class__.__name__)
            
            # end of main unwrap() function
            elif life_cycle == 'end':
                if self.operator_stack[-1] == val.__class__.__name__:
                    self.operator_stack.pop(-1)
            else:
                raise NotImplementedError


    ### HELPERS ###
    def should_apply_curly_braces(self, val: Any) -> bool:
        if self.in_shell_script:
            if len(self.operator_stack) == 1:
                if self.operator_stack[0] == val.__class__.__name__:
                    return True
        return False

    def should_quote(self, val: Any, expr: Any) -> bool:
        # only quote strings
        if not isinstance(val, str) or not isinstance(expr, str):
            return False
        # dont quote if already quoted
        if expr.startswith('"') or expr.endswith('"'):
            return False
        # # dont quote if wrapped isn't string & context is workflow
        # if not isinstance(val, str) and self.context == 'workflow':
        #     return False
        # master override - set when calling unwrap_expression.
        # some sort of external context means the expr should be quoted. 
        if self.quote_strings == True:
            return True
        if self.quote_strings == False:
            return False
        # string within curly braces
        if len(self.operator_stack) > 0:
            return True
        # stringformatter within shell script
        elif self.in_shell_script and isinstance(val, StringFormatter):
            return True
        return False

    def get_input_by_id(self, input_id: str) -> ToolInput:
        assert(self.tool is not None)
        inputs = [x for x in self.tool.inputs() if x.id() == input_id]
        return inputs[0]
    
    def get_channel_expression(self, channel_name: str, upstream_dtype: DataType) -> str:
        # scatter
        if self.scatter_target:
            # ch_bams -> ch_cartesian_cross.bams
            if self.scatter_method == ScatterMethod.cross:
                raise NotImplementedError
                # return cartesian_cross_subname(channel_name)
            # ch_bams -> ch_bams.flatten()
            elif self.scatter_method == ScatterMethod.dot:
                if upstream_dtype.is_array():
                    return f'{channel_name}.flatten()'
        # everything else
        return channel_name

    def build_resources_dict(self) -> dict[str, Any]:
        assert(self.tool)
        return {
            'runtime_cpu': self.tool.cpus({}),
            'runtime_memory': self.tool.memory({}),
            'runtime_seconds': self.tool.time({}),
            'runtime_disk': self.tool.disk({}),
        }


    ### LOGIC ###

    # primitives
    def unwrap_str(self, val: str) -> str:
        return val
    
    def unwrap_null(self, val: None) -> None:
        return None
    
    def unwrap_bool(self, val: bool) -> str:
        if val == True:
            return 'true'
        return 'false'

    def unwrap_int(self, val: int) -> str:
        return str(val)
    
    def unwrap_float(self, val: float) -> str:
        return str(val)
    
    def unwrap_list(self, val: list[Any]) -> str:
        elements: list[Any] = []
        for elem in val:
            el = self.unwrap(val=elem)
            elements.append(str(el))
        list_representation = f"[{', '.join(elements)}]"
        return list_representation


    # logical operators
    def unwrap_is_defined_operator(self, op: IsDefined) -> str:
        # this is a little weird. not a 1:1 mapping. 
        # assume everything is defined and set to null at least. 
        arg = self.unwrap(op.args[0])
        # return f"{arg} != null"
        return f"{arg}"
        
    def unwrap_if_operator(self, op: If) -> str:
        cond = self.unwrap(op.args[0])
        v1 = self.unwrap(op.args[1])
        v2 = self.unwrap(op.args[2])
        return f"{cond} ? {v1} : {v2}"

    def unwrap_assert_not_null_operator(self, op: AssertNotNull) -> str:
        arg = self.unwrap(op.args[0])
        return f'assert {arg[0]} != null'

    def unwrap_floor_operator(self, op: FloorOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"Math.floor({arg})"

    def unwrap_ceil_operator(self, op: CeilOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"Math.ceil({arg})"

    def unwrap_round_operator(self, op: RoundOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"Math.round({arg})"

    def unwrap_two_value_operator(self, op: TwoValueOperator) -> str:
        # self.add_quotes_to_strs = True
        arg1 = self.unwrap(op.args[0])
        arg2 = self.unwrap(op.args[1])
        # self.add_quotes_to_strs = False
        return f'{arg1} {op.nextflow_symbol()} {arg2}'


    # operator operators
    def unwrap_index_operator(self, op: IndexOperator) -> str:
        obj = op.args[0]
        index: int = op.args[1]  # type: ignore

        # edge case: IndexOperator(InputSelector(SecondaryType), index))
        if isinstance(obj, InputSelector):
            assert(self.tool)
            inp = self.tool.inputs_map()[obj.input_to_select]
            dtt = utils.get_dtt(inp.intype)
            if dtt in [ 
                DTypeType.SECONDARY_ARRAY,
                DTypeType.FILE_PAIR_ARRAY,
                DTypeType.FILE_ARRAY,
            ]:
                expr = self.unwrap_input_selector(obj, index=index)
                return f"{expr}[{index}]"
            
            elif dtt in [ 
                DTypeType.SECONDARY,
                DTypeType.FILE_PAIR,
            ]:
                expr = self.unwrap_input_selector(obj, index=index)
                return expr
            
            else:
                raise NotImplementedError

        # normal case
        expr = self.unwrap(obj)
        return f"{expr}[{index}]"
    
    def unwrap_as_string_operator(self, op: AsStringOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"{arg}.toString()"
    
    def unwrap_as_bool_operator(self, op: AsBoolOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"{arg}.toBoolean()"
    
    def unwrap_as_int_operator(self, op: AsIntOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"{arg}.toInteger()"
    
    def unwrap_as_float_operator(self, op: AsFloatOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"Float.valueOf({arg})"
    
    
    # standard operators
    def unwrap_read_contents_operator(self, op: ReadContents) -> str:
        arg = self.unwrap(op.args[0])
        return f"{arg}.text"
        
    def unwrap_read_json_operator(self, op: ReadJsonOperator) -> str:
        arg = self.unwrap(op.args[0])
        if isinstance(arg, str) and arg.startswith('"') and arg.endswith('"'):
            arg = arg.strip('"')
        return f'jsonSlurper.parseText(file("${{task.workDir}}/{arg}").text)'
    
    def unwrap_join_operator(self, op: JoinOperator) -> str:
        iterable = self.unwrap(op.args[0])
        separator = self.unwrap(op.args[1])
        return f"{iterable}.join({separator})"
    
    def unwrap_basename_operator(self, op: BasenameOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"{arg}.name"
    
    def unwrap_nameroot_operator(self, op: NamerootOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"{arg}.simpleName"
    
    def unwrap_nameext_operator(self, op: NameextOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"{arg}.extension"
    
    def unwrap_transpose_operator(self, op: TransposeOperator) -> str:
        raise NotImplementedError
    
    def unwrap_length_operator(self, op: LengthOperator) -> str:
        arg = self.unwrap(op.args[0])
        return f"{arg}.size"
    
    def unwrap_range_operator(self, op: RangeOperator) -> str:
        ceil = self.unwrap(op.args[0])
        return f'0..{ceil}'
    
    def unwrap_flatten_operator(self, op: FlattenOperator) -> str:
        # TODO VALIDATE
        arg = self.unwrap(op.args[0])
        return f'{arg}.flatten()'
    
    def unwrap_apply_prefix_operator(self, op: ApplyPrefixOperator) -> str:
        prefix = self.unwrap(op.args[0])
        iterable = self.unwrap(op.args[1])
        return f"{iterable}.map{{item -> {prefix} + item}}"
    
    def unwrap_file_size_operator(self, op: FileSizeOperator) -> str:
        fbytes = self.unwrap(op.args[0])
        return f"({fbytes}.size / 1048576)"
    
    def unwrap_first_operator(self, op: FirstOperator) -> str:
        resolved_list = self.unwrap(op.args[0])
        return f'{resolved_list}.find{{ it != null }}'
    
    def unwrap_filter_null_operator(self, op: FilterNullOperator) -> str:
        iterable = self.unwrap(op.args[0])
        return f"{iterable}.filter{{item -> item != null}}"
    
    def unwrap_replace_operator(self, op: ReplaceOperator) -> str:
        base = self.unwrap(op.args[0])
        pattern = self.unwrap(op.args[1])
        replacement = self.unwrap(op.args[2])
        return f"{base}.replaceAll({pattern}, {replacement})"

    # other operators
    def unwrap_operator(self, op: Operator) -> str:
        unwrap_expression_wrap = lambda x: unwrap_expression(
            val=x,
            # scope=self.scope,
            context=self.context,
            variable_manager=self.vmanager,
            tool=self.tool,
            in_shell_script=self.in_shell_script,

            scatter_target=self.scatter_target,
            scatter_method=self.scatter_method
        )
        return op.to_nextflow(unwrap_expression_wrap, *op.args)

    # selectors
    def unwrap_alias_selector(self, val: AliasSelector) -> Any:
        return self.unwrap(val.inner_selector)
    
    # misc
    def unwrap_string_formatter(self, selector: StringFormatter) -> str:
        """
        Translate Janis StringFormatter data type to Nextflow
        """
        if self.context == 'process_script' or self.context == 'process_output':
            assert(self.tool)
            assert(self.vmanager)
            expr = self.unwrap_string_formatter_process(selector)
        elif self.context == 'workflow':
            expr = self.unwrap_string_formatter_workflow(selector)
        else:
            raise RuntimeError
        return expr
        
    def unwrap_string_formatter_process(self, selector: StringFormatter) -> str:
        # assert(self.tool)  # n
        if len(selector.kwargs) == 0:
            return self.unwrap(selector._format)

        kwarg_replacements: dict[str, Any] = {}

        for k, v in selector.kwargs.items():
            kwarg_replacements[k] = self.unwrap(v)

        arg_val = selector._format
        for k in selector.kwargs:
            arg_val = arg_val.replace(f"{{{k}}}", f"{kwarg_replacements[k]}")

        if self.in_shell_script:
            arg_val = arg_val.replace('\\', '\\\\')

        return self.unwrap(arg_val)
        
    def unwrap_string_formatter_workflow(self, selector: StringFormatter) -> str:
        if len(selector.kwargs) == 0:
            return self.unwrap(selector._format)

        kwarg_replacements: dict[str, Any] = {}
        for k, v in selector.kwargs.items():
            kwarg_replacements[k] = self.unwrap(v)

        # workaround for channels 
        # need .first() at the moment to grab the value from a channel which only has a single value, 
        # but should guarantee in future this is correct for the channel
        for key, val in kwarg_replacements.items():
            if val.startswith('ch_'):
                kwarg_replacements[key] = f'{val}.first()'

        # reformat the selector format to be correct groovy syntax for use in workflow scope
        text_format = self.reformat_stringformatter_format_for_workflow_scope(selector._format)

        # substitute in unwrapped var values
        for k in selector.kwargs:
            text_format = text_format.replace(f"{{{k}}}", f"{kwarg_replacements[k]}")
        
        return self.unwrap(text_format)

    def reformat_stringformatter_format_for_workflow_scope(self, old_format: str) -> str:
        # reformat the selector format to be correct groovy syntax for use in workflow scope
        # eg: '{tumor}--{normal}' -> '{tumor} + "--" + {normal}'
        new_format: str = ''
        matches = re.findall(stringformatter_matcher, old_format)

        # replace each segment of the old_format, adding '+' and double quotes if needed
        for filler_text, var_text in matches:
            if filler_text != '':
                new_format += f' + "{filler_text}"'
            elif var_text != '':
                new_format += f' + {var_text}'
            else:
                raise RuntimeError
        
        # remove any beginning whitespace and '+'
        new_format = new_format.lstrip(' +')
        return new_format

    def unwrap_filename(self, fn: Filename, varname: Optional[str]=None) -> str:
        """
        order:
        prefix ref suffix extension
        ${outputFilename}.fastq.gz
        ${outputFilename}.fastq.gz
        ${"generated.fastq.gz"}
        etc
        """
        prefix = self.unwrap(fn.prefix) or ''
        varname = varname or ''
        suffix = self.unwrap(fn.suffix) or ''
        extension = self.unwrap(fn.extension) or ''

        # special prefix formatting - where ref is present
        if varname != '' and prefix.strip('"') == 'generated':
            prefix = ''

        # special suffix formatting
        quote_suffix = True if suffix.startswith('"') and suffix.endswith('"') else False
        suffix = suffix.strip('"')
        if suffix != '':
            if str(suffix).startswith("."):
                suffix = str(suffix)
            else:
                suffix = f'-{suffix}'
        if quote_suffix:
            suffix = f'"{suffix}"'

        # inside curly braces (each component wrapped in string)
        if len(self.operator_stack) > 0:
            items = [prefix, varname, suffix, extension]
            grouped_words = self.group_quoted_strings_in_list(items)
            expr = ' + '.join(grouped_words)
            return expr

        # not inside curly braces
        if varname != '':
            return f'{prefix}${{{varname}}}{suffix}{extension}'
        else:
            return f'{prefix}{suffix}{extension}'

    def group_quoted_strings_in_list(self, the_list: list[str]) -> list[str]:
        groups: list[str] = []
        current_group: list[str] = []
        
        for word in the_list:
            if word != '':
                if word.startswith('"') and word.endswith('"'):
                    current_group.append(word.strip('"'))
                else:
                    if len(current_group) > 0:
                        groups.append(f"\"{''.join(current_group)}\"")
                    groups.append(word)
                    current_group = []
        
        # still words in current_group by time we reach end of list
        if len(current_group) > 0:
            groups.append(f"\"{''.join(current_group)}\"")
        
        return groups


    ### PROCESSES ###

    def unwrap_wildcard_selector(self, sel: WildcardSelector) -> str:
        return f'{sel.wildcard}'
    
    def unwrap_input_selector(self, sel: InputSelector, index: Optional[int]=None) -> Optional[str]:
        """
        Translate Janis InputSelector data type into Nextflow expression
        I hate this function, and I am sorry. 
        """
        if not sel.input_to_select:
            raise Exception("No input was selected for input selector: " + str(sel))

        ### defining required variables ###
        inp = self.get_input_by_id(sel.input_to_select)
        dtype: DataType = inp.input_type # type: ignore
        basetype = utils.get_base_type(dtype)
        basetype = utils.ensure_single_type(basetype)

        ### getting the current variable for this tinput ###
        # (copy is made so can be modified locally without altering original)
        var = self.unwrap_input_selector_get_var(inp, index)

        ### resolving the tool input expression ###
        expr = self.unwrap_input_selector_get_expr(sel, inp, var)
        return expr
    
    def unwrap_input_selector_get_var(self, inp: ToolInput, index: Optional[int]=None) -> Variable:
        # returns a Variable which is the local nextflow variable for this input. 
        # depends on the type of TInput. 
        # this is because certain types (secondaries, file pairs) get reassigned as tuple process inputs
        # additionally, they may have extra prescript processing (get_primary_files, collate etc) to put them in proper format
        # finally, sometimes we might want to refer to the original process input variable, or an array joined local var 
        # depending on how the tinput is being accessed. 

        assert(self.vmanager)
        dtt = utils.get_dtt(inp.input_type)
        
        if self.context == 'process_script':
            if dtt == DTypeType.SECONDARY_ARRAY:
                if index is None:
                    var = self.vmanager.get(inp.id()).items[2]  # bams_joined
                    var_copy = deepcopy(var)
                else:
                    var = self.vmanager.get(inp.id()).items[1]  # bams
                    var_copy = deepcopy(var)

            elif dtt == DTypeType.SECONDARY:
                var = self.vmanager.get(inp.id()).items[1]
                var_copy = deepcopy(var)
            
            elif dtt == DTypeType.FILE_PAIR_ARRAY:
                var = self.vmanager.get(inp.id()).items[2]  # always read_pairs_joined
                var_copy = deepcopy(var)
            
            elif dtt == DTypeType.FILE_PAIR:
                if index is None:
                    var = self.vmanager.get(inp.id()).items[1]   # reads_joined
                    var_copy = deepcopy(var)
                else:
                    var = self.vmanager.get(inp.id()).original   # ['reads1', 'reads2']
                    var_copy = deepcopy(var)
                    var_copy.value = var_copy.value[index]   # reads1 or reads2
            
            elif dtt == DTypeType.FILE_ARRAY:
                if index is None:
                    var = self.vmanager.get(inp.id()).items[1]   # file_array_joined
                    var_copy = deepcopy(var)
                else:
                    var = self.vmanager.get(inp.id()).original   # file_array
                    var_copy = deepcopy(var)
            
            else:
                var = self.vmanager.get(inp.id()).current
                var_copy = deepcopy(var)
        
        elif self.context == 'process_output':
            # always referring to the original process input
            var = self.vmanager.get(inp.id()).original
            var_copy = deepcopy(var)
            # FILENAME ??

            if dtt == DTypeType.SECONDARY_ARRAY:
                # determine actual index in flat_array using index * num_files
                assert(index is not None)
                exts = utils.get_extensions(inp.input_type)
                actual_index = index * len(exts)
                var_copy.value = f'{var_copy.value}[{actual_index}]'

            elif dtt == DTypeType.SECONDARY:
                assert(index is None)
                var_copy.value = f'{var_copy.value}[0]'

            elif dtt == DTypeType.FILE_PAIR_ARRAY:
                # determine actual index in flat_array using index * 2
                assert(index is not None)
                actual_index = index * 2
                var_copy.value = f'{var_copy.value}[{actual_index}]'

            elif dtt == DTypeType.FILE_PAIR:
                if index:
                    var_copy.value = f'{var_copy.value}[{index}]'
                else:
                    var_copy.value = f'{var_copy.value}[0]'

            elif dtt == DTypeType.FILE_ARRAY:
                if index:
                    var_copy.value = f'{var_copy.value}[{index}]'

            elif dtt == DTypeType.FILE:
                pass
                # if index:
                #     var_copy.value = f'{var_copy.value}[{index}]'
                # assert(index is None)

            else:
                pass
                
        else:
            raise RuntimeError
        
        assert(var_copy)
        return var_copy

    def unwrap_input_selector_get_expr(self, sel: InputSelector, inp: ToolInput, var: Variable) -> Any:
        dtype: DataType = inp.input_type # type: ignore
        basetype = utils.get_base_type(dtype)
        basetype = utils.ensure_single_type(basetype)

        # special case: filename
        if isinstance(basetype, Filename):
            # super edge case - filename type referencing another input to generate name
            if var.value:
                assert(isinstance(var.value, str))
                expr = self.unwrap_filename(basetype, varname=var.value)
            else:
                expr = self.unwrap(basetype)

        # tinputs which had var available in the current scope
        elif var.vtype in [VariableType.TASK_INPUT, VariableType.PARAM, VariableType.LOCAL]:
            expr = var.value

        # tinputs which have static value
        elif var.vtype == VariableType.STATIC:
            expr = self.unwrap(var.value)
            print()

        # tinputs which are ignored in process but have default value
        elif var.vtype == VariableType.IGNORED and inp.default is not None:
            expr = self.unwrap(inp.default)
            print()
        
        # tinputs which are ignored in process and have no default value
        else:
            expr = None
            print()

        ### applying modifiers ###
        # special case: remove file extension
        if isinstance(basetype, File) and sel.remove_file_extension:
            # expr = f'{expr}.simpleName'
            expr = f'{expr}.baseName'
        
        return expr
        
    def unwrap_memory_selector(self, sel: MemorySelector) -> Any:
        assert(self.tool)
        return self.tool.memory({})

    def unwrap_cpu_selector(self, sel: CpuSelector) -> Any:
        assert(self.tool)
        return self.tool.cpus({})

    def unwrap_disk_selector(self, sel: DiskSelector) -> Any:
        assert(self.tool)
        return self.tool.disk({})

    def unwrap_time_selector(self, sel: TimeSelector) -> Any:
        assert(self.tool)
        return self.tool.time({})


    ### WORKFLOW PLUMBING ###

    def unwrap_input_node_selector(self, sel: InputNodeSelector) -> Any:
        return self.unwrap(sel.input_node)
            
    def unwrap_input_node(self, node: InputNode) -> Any:
        if self.context != 'workflow':
            raise RuntimeError
        assert(self.vmanager)
        cvar = self.vmanager.get(node.id()).current
        if cvar.vtype == VariableType.STATIC:
            # TODO debug make breakpoint. test = test_static_inputs()
            should_quote = self.should_quote(cvar.value, cvar.value)
            return nfgen_utils.to_groovy(cvar.value, quote_override=should_quote)
        else:
            qs_temp = deepcopy(self.quote_strings)
            self.quote_strings = False
            expr = self.unwrap(cvar.value)
            self.quote_strings = qs_temp
            return expr

    def unwrap_step_tag_input(self, val: StepTagInput) -> Any:
        # TODO save state of self.quote string
        return self.unwrap(val.source_map[0])
    
    def unwrap_edge(self, val: Edge) -> Any:
        return self.unwrap(val.source)

    def unwrap_step_output_selector(self, sel: StepOutputSelector) -> Any:
        # if scatter & output type is Array, use .flatten()
        upstream_step: StepNode = sel.node
        upstream_out: str = sel.tag
        conn_out = [x for x in upstream_step.tool.tool_outputs() if x.tag == upstream_out][0]

        # arrays of secondaries
        # @secondariesarray
        if utils.is_secondary_array_type(conn_out.outtype):
            raise NotImplementedError('process outputs with format [[file1, file2]] (arrays of secondary files) not supported in nextflow translation')

        # everything else
        else:
            upstream_step_id = naming.constructs.gen_varname_process(upstream_step.id())
            channel_name: str = f'{upstream_step_id}.out.{upstream_out}'
            return self.get_channel_expression(
                channel_name=channel_name,
                upstream_dtype=conn_out.outtype,
            )
