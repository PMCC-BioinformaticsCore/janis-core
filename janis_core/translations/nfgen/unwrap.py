
from typing import Any, Optional

from janis_core import Logger
from janis_core import (
    CommandTool, 
    TInput, 
    Operator, 
    AliasSelector, 
    InputSelector, 
    FirstOperator,
    WildcardSelector, 
    StringFormatter
)
from janis_core.graph.steptaginput import Edge, StepTagInput
from janis_core.operators.operator import IndexOperator
from janis_core.utils.scatter import ScatterMethod
from janis_core.workflow.workflow import StepNode, InputNode
from janis_core.types import (
    Filename,
    File,
    Directory,
    DataType
)
from janis_core.operators.selectors import InputNodeSelector, StepOutputSelector

from . import settings
from . import channels
from . import params
from . import nfgen_utils

from .scatter import cartesian_cross_subname
from .casefmt import to_case
from .process.inputs import create_inputs


def unwrap_expression(
    value: Any,
    
    tool: Optional[CommandTool]=None,
    quote_string: bool=True,
    for_output: bool=False,
    skip_inputs_lookup: bool=False,
    in_shell_script: bool=False, 
    
    sources: Optional[dict[str, Any]]=None,
    process_inputs: Optional[set[str]]=None,
    param_inputs: Optional[set[str]]=None,
    internal_inputs: Optional[set[str]]=None,

    scatter_target: bool=False,
    scatter_method: Optional[ScatterMethod]=None,
    
    # input_in_selectors: Optional[dict[str, Any]]=None,
    # inputs_dict: Optional[dict[str, TInput]]=None,
    # var_indicator: Optional[str]=None,
    # step_indicator: Optional[str]=None,
    # **debugkwargs: Any,
) -> Any:
    unwrapper = Unwrapper(
        # input_in_selectors=input_in_selectors,
        tool=tool,
        # inputs_dict=inputs_dict,
        quote_string=quote_string,
        for_output=for_output,
        skip_inputs_lookup=skip_inputs_lookup,
        in_shell_script=in_shell_script,

        sources=sources,
        process_inputs=process_inputs,
        param_inputs=param_inputs,
        internal_inputs=internal_inputs,

        scatter_target=scatter_target,
        scatter_method=scatter_method
    )
    return unwrapper.unwrap(value)


class Unwrapper:
    """
    The main logic to unwrap a janis expression and represent it in Nextflow translation
    """
    def __init__(
        self,
        # input_in_selectors: Optional[dict[str, Any]]=None,
        tool: Optional[CommandTool]=None,
        # inputs_dict: Optional[dict[str, TInput]]=None,
        quote_string: bool=True,
        for_output: bool=False,
        skip_inputs_lookup: bool=False,
        in_shell_script: bool=False, 

        sources: Optional[dict[str, Any]]=None,
        process_inputs: Optional[set[str]]=None,
        param_inputs: Optional[set[str]]=None,
        internal_inputs: Optional[set[str]]=None,

        scatter_target: bool=False,
        scatter_method: Optional[ScatterMethod]=None,
    ) -> None:
        # self.input_in_selectors = input_in_selectors if input_in_selectors else {}
        self.tool = tool
        # self.inputs_dict = inputs_dict
        self.quote_string = quote_string
        self.for_output = for_output
        self.skip_inputs_lookup = skip_inputs_lookup
        self.in_shell_script = in_shell_script

        self.sources = sources
        self.process_inputs = process_inputs
        self.param_inputs = param_inputs
        self.internal_inputs = internal_inputs

        self.scatter_target = scatter_target
        self.scatter_method = scatter_method

    def unwrap(self, val: Any) -> Any:
        ### ERRORS 

        if isinstance(val, StepNode):
            raise Exception(
                f"The Step node '{val.id()}' was found when unwrapping an expression, "
                f"you might not have selected an output."
            )

        ### INTERMEDIATE NODES 

        # step source
        if isinstance(val, StepTagInput):
            self.quote_string = False
            return self.unwrap(val.source_map[0])

        # edge
        elif isinstance(val, Edge):
            return self.unwrap(val.source)

        # alias
        elif isinstance(val, AliasSelector):
            return self.unwrap(val.inner_selector)
        
        # first
        elif isinstance(val, FirstOperator):
            # TODO implement properly
            return self.unwrap(val.args[0][0])
        
        # index
        elif isinstance(val, IndexOperator):
            return self.unwrap_index_operator(val)

        ### LEAF NODES

        # None
        elif val is None:
            if self.quote_string:
                return "null"
            return None

        # str
        elif isinstance(val, str):
            if self.quote_string:
                return f'"{val}"'
            return val

        # bool
        elif isinstance(val, bool):
            if self.quote_string:
                return f'"{val}"'
            return val

        # numeric
        elif isinstance(val, int) or isinstance(val, float):
            return str(val)

        # list
        if isinstance(val, list):
            # toolid = str(debugkwargs.get("tool_id", "unwrap_list_expression"))
            elements: list[Any] = []
            for elem in val:
                el = self.unwrap(val=elem)
                elements.append(el)
            list_representation = f"[{', '.join(elements)}]"
            return list_representation

        # Filename
        elif isinstance(val, Filename):
            formatted = val.generated_filename()
            return formatted

        elif isinstance(val, StringFormatter):
            assert(self.tool)
            # self.in_shell_script = False
            # self.skip_inputs_lookup = False
            return self.translate_string_formatter(val)

        # tool input
        elif isinstance(val, InputSelector):
            """
            TODO: edge case for secondaries 

            ToolInput("reads", FastqGzPair),
            ToolInput(
                "read1",
                FastqGz(optional=True),
                default=IndexOperator(InputSelector("reads"), 0),
                position=5,
            ),
            ToolInput(
                "read2",
                FastqGz(optional=True),
                default=IndexOperator(InputSelector("reads"), 1),
                position=6,
            ),

            janis: ToolInput("reads", FastqGzPair)
            nextflow: tuple path(reads1), path(reads2)

            unwrap(IndexOperator(InputSelector("reads"), 0)) -> "reads1"
            unwrap(IndexOperator(InputSelector("reads"), 1)) -> "reads2"

            """
            if self.for_output:
                return self.unwrap_input_selector_output(sel=val)
            else:
                return self.unwrap_input_selector(sel=val)

        # workflow input
        elif isinstance(val, InputNodeSelector):
            return self.unwrap_input_node(val.input_node)

        # step output
        elif isinstance(val, StepOutputSelector):
            return self.unwrap_connection(val)
            # step_name = to_case(val.node.id(), settings.NEXTFLOW_PROCESS_CASE)
            # return f'{step_name}.out.{val.tag}'

        elif isinstance(val, WildcardSelector):
            return f'"{val.wildcard}"'

        # any other Operators
        elif isinstance(val, Operator):
            return self.unwrap_operator(val)
        
        # anything else with a .to_nextflow() method
        elif callable(getattr(val, "to_nextflow", None)):
            return val.to_nextflow()
            # if var_indicator is not None and step_indicator is not None:
            #     return value.to_nextflow(
            #         var_indicator=var_indicator, step_indicator=step_indicator
            #     )
            # else:

        raise Exception(
            "Could not detect type %s to convert to input value" % type(val)
        )

    ### HELPER METHODS

    def get_src_variable(self, inp: TInput) -> Optional[str]:
        assert(self.process_inputs is not None)
        assert(self.param_inputs is not None)
        assert(self.internal_inputs is not None)
        
        if inp.id() in self.process_inputs:
            src = self.get_src_process_input(inp)
        
        elif inp.id() in self.param_inputs:
            src = self.get_src_param_input(inp)
        
        elif inp.id() in self.internal_inputs and inp.default is not None:
            src = self.unwrap(inp.default)

        elif inp.id() in self.internal_inputs:
            src = None
        
        else:
            raise NotImplementedError

        return src

    def get_src_process_input(self, inp: TInput) -> str:
        # data fed via process input
        dtype = inp.intype # type: ignore
        basetype = nfgen_utils.get_base_type(dtype)
        # secondary files (name mapped to ext of primary file)
        if isinstance(basetype, File) and basetype.has_secondary_files():
            exts = nfgen_utils.get_extensions(basetype)
            name = exts[0]
        # everything else
        else:
            name = inp.id()
        return name

    def get_src_param_input(self, inp: TInput) -> str: 
        # data fed via global param
        assert(self.sources is not None)
        src = self.sources[inp.id()]
        sel = src.source_map[0].source
        param = params.get(sel.input_node.uuid)
        return f'params.{param.name}'

    def get_input_by_id(self, inname: str) -> TInput:
        assert(self.tool is not None)
        inputs_dict = self.tool.inputs_map()
        inp = inputs_dict[inname]
        return inp

    def unwrap_index_operator(self, op: IndexOperator) -> Any:
        self.quote_string = False
        
        # special case: janis secondary -> nextflow tuple
        if isinstance(op.args[0], InputSelector):
            sel = op.args[0]
            inp = self.get_input_by_id(sel.input_to_select)
            if nfgen_utils.is_secondary_type(inp.intype):
                process_in = create_inputs(inp)[0]
                print()

        # special case: janis array secondary -> multiple nextflow path
        # maybe this is just in unwrap_input_selector()

        # everything else
        index = op.args[1]
        expr = self.unwrap(op.args[0])
        return f"{expr}[{index}]"

    def unwrap_connection(self, sel: StepOutputSelector) -> Any:
        # if scatter & output type is Array, use .flatten()
        args: list[str] = []
        upstream_step: StepNode = sel.node
        upstream_out: str = sel.tag
        conn_out = [x for x in upstream_step.tool.tool_outputs() if x.tag == upstream_out][0]

        # arrays of secondaries
        if nfgen_utils.is_array_secondary_type(conn_out.outtype):
            raise NotImplementedError
        else:
            upstream_step_id = to_case(upstream_step.id(), settings.NEXTFLOW_PROCESS_CASE)
            channel_name: str = f'{upstream_step_id}.out.{upstream_out}'
            arg = self.get_channel_expression(
                    channel_name=channel_name,
                    upstream_dtype=conn_out.outtype,
                )
            args.append(arg)
        return args

    def unwrap_input_node(self, node: InputNode) -> Any:
        if channels.exists(janis_uuid=node.uuid):
            return self.unwrap_channel(node)
        
        elif params.exists(janis_uuid=node.uuid):
            param = params.get(janis_uuid=node.uuid)
            return f'params.{param.name}'

    def unwrap_channel(self, node: InputNode) -> Any:
        """
        ch_name                     = same type (most cases)
        ch_name.collect()           = singles -> array (workflow input array channel creation)
        ch_name.flatten()           = array -> singles (scatter.dot)
        cartesian_cross.ch_subname  = scatter.cross  
        """
        args: list[str] = []
        relevant_channels = channels.getall(janis_uuid=node.uuid)
        
        # arrays of secondaries
        if relevant_channels and len(relevant_channels) > 1:
            for ch in relevant_channels:
                arg = self.get_channel_expression(
                    channel_name=ch.name,
                    upstream_dtype=node.datatype,
                )
                args.append(arg)
        
        # everything else
        elif relevant_channels and len(relevant_channels) == 1:
            arg = self.get_channel_expression(
                channel_name=relevant_channels[0].name,
                upstream_dtype=node.datatype,
            )
            args.append(arg)
        return args

    def get_channel_expression(self, channel_name: str, upstream_dtype: DataType) -> str:
        # scatter
        if self.scatter_target:
            # ch_bams -> ch_cartesian_cross.bams
            if self.scatter_method == ScatterMethod.cross:
                return cartesian_cross_subname(channel_name)
            # ch_bams -> ch_bams.flatten()
            elif self.scatter_method == ScatterMethod.dot:
                if upstream_dtype.is_array():
                    return f'{channel_name}.flatten()'
        # everything else
        return channel_name

    def unwrap_operator(self, operator: Operator) -> Any:
        unwrap_expression_wrap = lambda x: unwrap_expression(
            value=x,
            tool=self.tool,
            # inputs_dict=self.inputs_dict,
            quote_string=self.quote_string,
            for_output=self.for_output,
            skip_inputs_lookup=self.skip_inputs_lookup,
            in_shell_script=self.in_shell_script,

            sources=self.sources,
            process_inputs=self.process_inputs,
            param_inputs=self.param_inputs,
            internal_inputs=self.internal_inputs,

            scatter_target=self.scatter_target,
            scatter_method=self.scatter_method

            # # input_in_selectors=self.input_in_selectors,
            # quote_string=self.quote_string,
            # tool=self.tool,
            # for_output=self.for_output,
            # # inputs_dict=self.inputs_dict,
            # skip_inputs_lookup=self.skip_inputs_lookup,
            # in_shell_script=self.in_shell_script,
            # # **debugkwargs,
        )
        return operator.to_nextflow(unwrap_expression_wrap, *operator.args)


    def translate_string_formatter(self, selector: StringFormatter) -> str:
        """
        Translate Janis StringFormatter data type to Nextflow
        """
        if len(selector.kwargs) == 0:
            return str(selector)

        kwarg_replacements: dict[str, Any] = {}
        for k, v in selector.kwargs.items():
            kwarg_replacements[k] = self.unwrap(v)

        arg_val = selector._format
        for k in selector.kwargs:
            arg_val = arg_val.replace(f"{{{k}}}", f"${{{str(kwarg_replacements[k])}}}")

        if self.in_shell_script:
            arg_val = arg_val.replace("\\", "\\\\")

        return arg_val


    def unwrap_input_selector_output(self, sel: InputSelector) -> Optional[str]:
        """
        Generate a string expression to represent a filename in Nextflow
        """
        inp = self.get_input_by_id(sel.input_to_select)
        src = self.get_src_variable(inp)
        if src is None:
            return None

        dtype: DataType = inp.intype # type: ignore

        if nfgen_utils.is_array_secondary_type(dtype):
            raise NotImplementedError
        if nfgen_utils.is_secondary_type(dtype):
            raise NotImplementedError
        
        # HERE
        if not self.tool:
            raise NotImplementedError
            return f"${sel.input_to_select}.name"
            # raise Exception(
            #     f"Couldn't generate filename as an internal error occurred (self.inputs_dict did not contain {inp.input_to_select})"
            # )
        
        # if sel.input_to_select not in inputs_dict:
        #     raise Exception(
        #         f"The InputSelector '{sel.input_to_select}' did not select a valid input"
        #     )

        if isinstance(dtype, File):
            base = src

            if dtype.has_secondary_files():
                base = f"{base}[0]"

            if sel.remove_file_extension and dtype.get_extensions():
                base = f"{base}.simpleName"

            # elif hasattr(inp, "localise_file") and inp.localise_file:
            #     base = f"{base}.name"

        elif isinstance(dtype, Filename):
            self.for_output = True
            base = f"{self.unwrap(dtype)}"
        
        # elif (
        #     dtype.is_array()
        #     and isinstance(dtype.fundamental_type(), (File, Directory))
        #     and hasattr(inp, "localise_file")
        #     and inp.localise_file
        # ):
        #     base = f"{src}.map{{ el.name }}"
        
        else:
            base = f"{src}"

        # if intype.optional:
        #     replacement = f'{base}, optional: true'
        #     # default = '"generated"'
        #     # if isinstance(intype, Filename):
        #     #     default = base
        #     # replacement = f"({sel.input_to_select} && {sel.input_to_select} != 'None' && {sel.input_to_select} != '' ? {sel.input_to_select} : {default})"
        # else:
        # replacement = f'{base}'

        # return f"\"${{{replacement}}}\""
        return base


    def unwrap_input_selector(self, sel: InputSelector) -> Optional[str]:
        """
        Translate Janis InputSelector data type into Nextflow expressions
        """
        inp = self.get_input_by_id(sel.input_to_select)
        src = self.get_src_variable(inp)
        if src is None:
            return None

        # HERE
        inputs_dict = self.tool.inputs_map()

        # if self.tool.versioned_id() not in self.input_in_selectors:
        #     self.input_in_selectors[self.tool.versioned_id()] = set()

        # self.input_in_selectors[self.tool.versioned_id()].add(sel.input_to_select)

        expr: str = sel.input_to_select
        if not expr:
            raise Exception(
                "No input was selected for input selector: " + str(sel)
                # Why is this here???? GH
            )

        skip_lookup = self.skip_inputs_lookup or expr.startswith("runtime_")

        if not skip_lookup:

            if inputs_dict is None:
                raise Exception(
                    f"An internal error occurred when translating input selector '{expr}': the inputs dictionary was None"
                )
            if sel.input_to_select not in inputs_dict:
                raise Exception(
                    f"Couldn't find the input '{expr}' for the InputSelector(\"{expr}\")"
                )

            tinp: TInput = inputs_dict[sel.input_to_select]

            intype = tinp.intype

            if intype.is_base_type((File, Directory)):
                if intype.has_secondary_files():
                    # TODO HERE
                    expr = f"{expr}[0]"

            if sel.remove_file_extension:
                if intype.is_base_type((File, Directory)):

                    potential_extensions = (
                        intype.get_extensions() if intype.is_base_type(File) else None
                    )
                    if sel.remove_file_extension and potential_extensions:
                        expr = f"{expr}.simpleName"

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
                    expr = f"{expr}.map(function(el) {{ return {inner_sel}; }})"
                else:
                    Logger.warn(
                        f"InputSelector {expr} is requesting to remove_file_extension but it has type {tinp.intype.id()}"
                    )

        if self.in_shell_script:
            expr = f"${{{expr}}}"

        return expr