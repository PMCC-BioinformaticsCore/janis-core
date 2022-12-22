

from typing import Any, Optional
from collections import defaultdict

from janis_core.graph.steptaginput import Edge, StepTagInput
from janis_core.operators.operator import IndexOperator, Operator
from janis_core.types import Filename
from janis_core.operators.logical import (
    IsDefined,
    If,
    AssertNotNull,
    FloorOperator,
    CeilOperator,
    RoundOperator,
)
from janis_core.operators.standard import (
    ReadContents,
    ReadJsonOperator,
    JoinOperator,
    BasenameOperator,
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
)
from janis_core.operators.stringformatter import StringFormatter
from janis_core import CommandTool

def trace_janis_entities(entity: Any, tool: Optional[CommandTool]=None) -> dict[str, int]:
    tracer = EntityTracer(tool)
    tracer.trace(entity)
    return tracer.counter


class EntityTracer:
    
    def __init__(self, tool: Optional[CommandTool]=None):
        self.tool = tool
        self.counter: dict[str, int] = defaultdict(int)
        self.single_arg_trace_types = {
            IsDefined,
            AssertNotNull,
            FloorOperator,
            CeilOperator,
            RoundOperator,
            AsStringOperator,
            AsBoolOperator,
            AsIntOperator,
            AsFloatOperator,
            ReadContents,
            ReadJsonOperator,
            BasenameOperator,
            TransposeOperator,
            LengthOperator,
            RangeOperator,
            FlattenOperator,
            FileSizeOperator,
            FirstOperator,
            FilterNullOperator,
        }
        self.multi_arg_trace_types = {
            If,
            IndexOperator,
            JoinOperator,
            ApplyPrefixOperator,
            ReplaceOperator,
        }
        self.custom_trace_funcs = {

            # logical operators
            # IsDefined: self.is_defined_operator,
            # If: self.if_operator,
            # AssertNotNull: self.assert_not_null_operator,
            # FloorOperator: self.floor_operator,
            # CeilOperator: self.ceil_operator,
            # RoundOperator: self.round_operator,

            # operator operators
            # IndexOperator: self.index_operator,
            # AsStringOperator: self.as_string_operator,
            # AsBoolOperator: self.as_bool_operator,
            # AsIntOperator: self.as_int_operator,
            # AsFloatOperator: self.as_float_operator,
            
            # standard operators
            # ReadContents: self.read_contents_operator,
            # ReadJsonOperator: self.read_json_operator,
            # JoinOperator: self.join_operator,
            # BasenameOperator: self.basename_operator,
            # TransposeOperator: self.transpose_operator,
            # LengthOperator: self.length_operator,
            # RangeOperator: self.range_operator,
            # FlattenOperator: self.flatten_operator,
            # ApplyPrefixOperator: self.apply_prefix_operator,
            # FileSizeOperator: self.file_size_operator,
            # FirstOperator: self.first_operator,
            # FilterNullOperator: self.filter_null_operator,
            # ReplaceOperator: self.replace_operator,

            # primitives
            list: self.trace_list,
            
            # selectors
            AliasSelector: self.alias_selector,
            InputNodeSelector: self.input_node_selector,
            WildcardSelector: self.wildcard_selector,
            InputSelector: self.input_selector,
            StepOutputSelector: self.step_output_selector,
            MemorySelector: self.memory_selector,
            CpuSelector: self.cpu_selector,
            DiskSelector: self.disk_selector,
            TimeSelector: self.time_selector,

            # misc 
            StepTagInput: self.step_tag_input,
            Edge: self.edge,
            StringFormatter: self.string_formatter,
            Filename: self.filename,
            # InputNode: self.input_node,
        }

    
    def trace(self, entity: Any) -> None:
        etype = type(entity)
        ename = entity.__class__.__name__
        self.counter[ename] += 1

        if etype in self.single_arg_trace_types:
            self.operator_single_arg_trace(entity)
        
        elif etype in self.multi_arg_trace_types:
            self.operator_multi_arg_trace(entity)

        elif etype in self.custom_trace_funcs:
            func = self.custom_trace_funcs[etype]
            func(entity)
        


    def operator_single_arg_trace(self, entity: Operator) -> None:
        self.trace(entity.args[0])
    
    def operator_multi_arg_trace(self, entity: Operator) -> None:
        for arg in entity.args:
            self.trace(arg)

    def trace_list(self, entity: list[Any]) -> None:
        for item in entity:
            self.trace(item)


    # def is_defined_operator(self, entity: IsDefined) -> None:
    #     self.trace(entity.args[0])

    # def if_operator(self, entity: If) -> None:
    #     for arg in entity.args:
    #         self.trace(arg)

    # def assert_not_null_operator(self, entity: AssertNotNull) -> None:
    #     self.trace(entity.args[0])

    # def floor_operator(self, entity: FloorOperator) -> None:
    #     self.trace(entity.args[0])

    # def ceil_operator(self, entity: CeilOperator) -> None:
    #     self.trace(entity.args[0])

    # def round_operator(self, entity: RoundOperator) -> None:
    #     self.trace(entity.args[0])

    # def index_operator(self, entity: IndexOperator) -> None:
    #     for arg in entity.args:
    #         self.trace(arg)

    # def as_string_operator(self, entity: AsStringOperator) -> None:
    #     self.trace(entity.args[0])

    # def as_bool_operator(self, entity: AsBoolOperator) -> None:
    #     self.trace(entity.args[0])

    # def as_int_operator(self, entity: AsIntOperator) -> None:
    #     self.trace(entity.args[0])

    # def as_float_operator(self, entity: AsFloatOperator) -> None:
    #     self.trace(entity.args[0])

    # def read_contents_operator(self, entity: ReadContents) -> None:
    #     self.trace(entity.args[0])

    # def read_json_operator(self, entity: ReadJsonOperator) -> None:
    #     self.trace(entity.args[0])

    # def join_operator(self, entity: JoinOperator) -> None:
    #     for arg in entity.args:
    #         self.trace(arg)

    # def basename_operator(self, entity: BasenameOperator) -> None:
    #     self.trace(entity.args[0])

    # def transpose_operator(self, entity: TransposeOperator) -> None:
    #     self.trace(entity.args[0])

    # def length_operator(self, entity: LengthOperator) -> None:
    #     self.trace(entity.args[0])

    # def range_operator(self, entity: RangeOperator) -> None:
    #     self.trace(entity.args[0])

    # def flatten_operator(self, entity: FlattenOperator) -> None:
    #     self.trace(entity.args[0])

    # def apply_prefix_operator(self, entity: ApplyPrefixOperator) -> None:
    #     for arg in entity.args:
    #         self.trace(arg)

    # def file_size_operator(self, entity: FileSizeOperator) -> None:
    #     self.trace(entity.args[0])

    # def first_operator(self, entity: FirstOperator) -> None:
    #     self.trace(entity.args[0])

    # def filter_null_operator(self, entity: FilterNullOperator) -> None:
    #     self.trace(entity.args[0])

    # def replace_operator(self, entity: ReplaceOperator) -> None:
    #     for arg in entity.args:
    #         self.trace(arg)

    def alias_selector(self, entity: AliasSelector) -> None:
        self.trace(entity.inner_selector)
        self.trace(entity.data_type)

    def input_node_selector(self, entity: InputNodeSelector) -> None:
        self.trace(entity.input_node)

    def wildcard_selector(self, entity: WildcardSelector) -> None:
        # TODO check
        self.trace(entity.wildcard)

    def input_selector(self, entity: InputSelector) -> None:
        assert(self.tool)
        tinput = self.tool.inputs_map()[entity.input_to_select]
        
        # the string used as InputSelector reference
        self.trace(entity.input_to_select)
        # the datatype if it is a Filename type
        if isinstance(tinput.intype, Filename):
            self.trace(tinput.intype)

    def step_output_selector(self, entity: StepOutputSelector) -> None:
        # TODO check
        self.trace(entity.node)
        self.trace(entity.tag)

    def memory_selector(self, entity: MemorySelector) -> None:
        self.trace(entity.resource_to_select)
        self.trace(entity.resource_type)
        self.trace(entity.default)

    def cpu_selector(self, entity: CpuSelector) -> None:
        self.trace(entity.resource_to_select)
        self.trace(entity.resource_type)
        self.trace(entity.default)

    def disk_selector(self, entity: DiskSelector) -> None:
        self.trace(entity.resource_to_select)
        self.trace(entity.resource_type)
        self.trace(entity.default)

    def time_selector(self, entity: TimeSelector) -> None:
        self.trace(entity.resource_to_select)
        self.trace(entity.resource_type)
        self.trace(entity.default)

    def step_tag_input(self, entity: StepTagInput) -> None:
        for src in entity.source_map:
            self.trace(src)

    def edge(self, entity: Edge) -> None:
        self.trace(entity.source)

    def string_formatter(self, entity: StringFormatter) -> None:
        # trace the value of each keyword
        for item in entity.kwargs.values():
            self.trace(item)

    def filename(self, entity: Filename) -> None:
        self.trace(entity.prefix)
        self.trace(entity.suffix)
        self.trace(entity.extension)


 




