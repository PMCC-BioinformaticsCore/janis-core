

from __future__ import annotations
# from typing import TYPE_CHECKING
# if TYPE_CHECKING:
#     from .parsing.process.VariableManager import VariableManager

from enum import Enum

from .variables import VariableManager
from .variables import VariableType
from .variables import Variable

from copy import deepcopy
from typing import Any, Optional
NoneType = type(None)
import re

from janis_core import (
    CommandToolBuilder, 
    ToolInput, 
    Operator,
    ReadJsonOperator,
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

from janis_core.messages import log_message
from janis_core.messages import ErrorCategory
from janis_core.operators.operator import (
    IndexOperator,
    AsStringOperator,
    AsBoolOperator,
    AsIntOperator,
    AsFloatOperator,
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


"""
BRACES 

1.  StringFormatter(format{token1} some text {token2}, token1=[1, 2, 3], token2=4)
    - expr: {token1} some text {token2} - apply_braces=False, inside_braces=False, outcome=False
    - expr: [1, 2, 3]                   - apply_braces=True,  inside_braces=False, outcome=True
    - expr: 1                           - apply_braces=False, inside_braces=True,  outcome=False

"""


def unwrap_expression(
    val: Any,
    context: str,
    variable_manager: Optional[VariableManager]=None,
    tool: Optional[CommandToolBuilder]=None,
    
    apply_braces: bool=False,
    strquote_override: Optional[bool]=None,
    
    scatter_target: bool=False,
    scatter_method: Optional[ScatterMethod]=None,
    ) -> Any:

    context_e = UnwrapContext(context)
    unwrapper = Unwrapper(
        context_e=context_e,
        variable_manager=variable_manager,
        tool=tool,
        apply_braces=apply_braces,
        strquote_override=strquote_override,
        scatter_target=scatter_target,
        scatter_method=scatter_method,
    )
    return unwrapper.unwrap(val)


class UnwrapContext(Enum):
    DEFAULT             = 'default'
    PROCESS_PRESCRIPT   = 'process_prescript'
    PROCESS_SCRIPT      = 'process_script'
    PROCESS_OUTPUT      = 'process_output'
    WORKFLOW            = 'workflow'


class Unwrapper:
    """
    The main logic to unwrap a janis expression and represent it in Nextflow translation
    """
    def __init__(
        self,
        context_e: UnwrapContext,
        variable_manager: Optional[VariableManager],
        tool: Optional[CommandToolBuilder],
        apply_braces: bool,
        strquote_override: Optional[bool],
        scatter_target: bool,
        scatter_method: Optional[ScatterMethod], 
    ) -> None:
        self.context_e = context_e
        self.vmanager = variable_manager
        self.tool = tool
        self.apply_braces = apply_braces
        self.strquote_override = strquote_override
        self.scatter_target = scatter_target
        self.scatter_method = scatter_method

        self.unwrap_overrides = {
            # primitives
            NoneType: self.unwrap_null,
            str: self.unwrap_str,
            bool: self.unwrap_bool,
            int: self.unwrap_int,
            float: self.unwrap_float,
            list: self.unwrap_list,

            # operator operators
            IndexOperator: self.unwrap_index_operator,
            AsStringOperator: self.unwrap_as_string_operator,
            AsBoolOperator: self.unwrap_as_bool_operator,
            AsIntOperator: self.unwrap_as_int_operator,
            AsFloatOperator: self.unwrap_as_float_operator,
            ReadJsonOperator: self.unwrap_read_json_operator,
            
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
        # do the unwrap  
        expr = self.unwrap_entity(val)

        # apply quotes if needed
        if self.should_quote(val, expr):
            expr = f'"{expr}"'

        # apply nextflow curly braces if needed
        if self.apply_braces:
            expr = f'${{{expr}}}'
        
        return expr
    
    def unwrap_entity(self, entity: Any) -> Any:
        """switchboard to delegate unwrap based on entity type."""

        # entity has custom unwrap function override
        etype = type(entity)
        if etype in self.unwrap_overrides:
            func = self.unwrap_overrides[etype]
            expr = func(entity)

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

    ### HELPERS ###

    def should_quote(self, val: Any, expr: Any) -> bool:
        if not isinstance(val, str):
            return False
        # only quote strings
        if self.strquote_override is not None:
            return self.strquote_override
        # dont quote if already quoted
        PATTERN = r'(^\'[\s\S]*\'$)|(^"[\s\S]*"$)'
        if re.match(PATTERN, expr):
            return False
        return True

    def get_input_by_id(self, input_id: str) -> ToolInput:
        assert(self.tool is not None)
        inputs = [x for x in self.tool.inputs() if x.id() == input_id]
        # if not inputs:
        #     msg = f"Could not find input with id '{input_id}' in tool '{self.tool.id()}'"
        #     log_message(self.tool.uuid, msg, ErrorCategory.PLUMBING)
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

    ### UNWRAP OVERRIDES ###

    # primitives
    def unwrap_str(self, val: str) -> str:
        return val
    
    def unwrap_null(self, val: None) -> str:
        return 'null'
    
    def unwrap_bool(self, val: bool) -> str:
        if val == True:
            return 'true'
        return 'false'

    def unwrap_int(self, val: int) -> str:
        return str(val)
    
    def unwrap_float(self, val: float) -> str:
        return str(val)
    
    def unwrap_list(self, val: list[Any]) -> str:
        ab_temp = deepcopy(self.apply_braces)
        self.apply_braces = False
        elements: list[Any] = []
        for elem in val:
            el = self.unwrap(val=elem)
            elements.append(str(el))
        list_representation = f"[{', '.join(elements)}]"
        self.apply_braces = ab_temp
        return list_representation

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
    
    def unwrap_read_json_operator(self, op: ReadJsonOperator) -> str:
        arg = unwrap_expression(
            op.args[0], 
            context=self.context_e.value,
            variable_manager=self.vmanager,
            tool=self.tool,
            apply_braces=True,
            strquote_override=True
        )
        return f'jsonSlurper.parseText(file("${{task.workDir}}/{arg}").text)'

    # other operators
    def unwrap_operator(self, op: Operator) -> str:

        unwrap_expression_wrap = lambda x: unwrap_expression(
            val=x,
            context=self.context_e.value,
            variable_manager=self.vmanager,
            tool=self.tool,

            scatter_target=self.scatter_target,
            scatter_method=self.scatter_method,
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
        return self.unwrap_string_formatter_process(selector)
        
    def unwrap_string_formatter_process(self, selector: StringFormatter) -> str:
        # assert(self.tool)  # n
        if len(selector.kwargs) == 0:
            return self.unwrap(selector._format)

        kwarg_replacements: dict[str, Any] = {}

        for k, v in selector.kwargs.items():
            # no need to add curly braces, don't quote the string. 
            if isinstance(v, str):
                kwarg_replacements[k] = v

            # must add curly braces & unwrap
            else:
                kwarg_replacements[k] = unwrap_expression(
                    v, 
                    context=self.context_e.value, 
                    variable_manager=self.vmanager,
                    tool=self.tool,
                    apply_braces=True,
                    strquote_override=True,
                )

        arg_val = selector._format
        for k in selector.kwargs:
            arg_val = arg_val.replace(f"{{{k}}}", f"{kwarg_replacements[k]}")

        if self.context_e == UnwrapContext.PROCESS_SCRIPT:
            arg_val = arg_val.replace('\\', '\\\\')

        return arg_val

    def unwrap_filename(self, fn: Filename, varname: Optional[str]=None) -> str:
        """
        order:
        prefix ref suffix extension
        ${outputFilename}.fastq.gz
        ${outputFilename}.fastq.gz
        ${"generated.fastq.gz"}
        etc
        """
        if any([isinstance(x, Selector) for x in [fn.prefix, fn.suffix, fn.extension]]):
            kwargs = {}
            fmt = ''
            if fn.prefix is not None:
                if isinstance(fn.prefix, str):
                    fmt += f'{fn.prefix}'
                else:
                    fmt += '{prefix}'
                    kwargs['prefix']=fn.prefix
            if fn.suffix is not None:
                if isinstance(fn.suffix, str):
                    fmt += f'{fn.suffix}'
                else:
                    fmt += '{suffix}'
                    kwargs['suffix']=fn.suffix
            if fn.extension is not None:
                if isinstance(fn.extension, str):
                    fmt += f'{fn.extension}'
                else:
                    fmt += '{extension}'
                    kwargs['extension']=fn.extension

            val = StringFormatter(fmt, **kwargs)
            expr = unwrap_expression(
                val=val,
                context=self.context_e.value,
                variable_manager=self.vmanager,
                tool=self.tool,
            )
        else:
            expr = fn.generated_filename()
        return expr
    

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
        
        if self.context_e == UnwrapContext.PROCESS_SCRIPT:
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
            
            elif dtt == DTypeType.FILE_PAIR:
                if index is None:
                    var = self.vmanager.get(inp.id()).items[1]   # reads_joined
                    var_copy = deepcopy(var)
                else:
                    var = self.vmanager.get(inp.id()).original   # ['reads1', 'reads2']
                    var_copy = deepcopy(var)
                    var_copy.value = var_copy.value[index]   # reads1 or reads2
            
            elif dtt in [DTypeType.FILE_PAIR_ARRAY, DTypeType.FILE_ARRAY]:
                var = self.vmanager.get(inp.id()).original   
                var_copy = deepcopy(var)
            
            else:
                var = self.vmanager.get(inp.id()).current
                var_copy = deepcopy(var)
        
        elif self.context_e == UnwrapContext.PROCESS_PRESCRIPT or self.context_e == UnwrapContext.PROCESS_OUTPUT:
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
            # print()

        # tinputs which are ignored in process but have default value
        elif var.vtype == VariableType.IGNORED and inp.default is not None:
            expr = self.unwrap(inp.default)
            # print()
        
        # tinputs which are ignored in process and have no default value
        else:
            expr = None
            # print()

        ### applying modifiers ###
        # special case: remove file extension
        if isinstance(basetype, File) and sel.remove_file_extension:
            # expr = f'{expr}.baseName'
            expr = f'{expr}.simpleName'
        
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
        if self.context_e != UnwrapContext.WORKFLOW:
            raise RuntimeError
        assert(self.vmanager)
        cvar = self.vmanager.get(node.id()).current
        if cvar.vtype == VariableType.STATIC:
            return nfgen_utils.to_groovy(cvar.value, quote_override=self.strquote_override)
        else:
            qs_temp = deepcopy(self.strquote_override)
            self.strquote_override = False
            expr = self.unwrap(cvar.value)
            self.strquote_override = qs_temp
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
