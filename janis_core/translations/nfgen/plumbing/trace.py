

from abc import ABC, abstractmethod
from typing import Any, Optional
from collections import defaultdict


from janis_core.graph.steptaginput import Edge, StepTagInput
from janis_core.operators.operator import IndexOperator, Operator
from janis_core.types import Filename, DataType, Stdout
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
    ResourceSelector
)
from janis_core.operators.stringformatter import StringFormatter
from janis_core import CommandTool, ToolInput, TInput, ToolArgument, ToolOutput

from .. import nfgen_utils



def trace_entity_counts(entity: Any, tool: Optional[CommandTool]=None) -> dict[str, int]:
    tracer = EntityCountTracer(tool)
    tracer.trace(entity)
    return tracer.counter

def trace_source_datatype(entity: Any, tool: Optional[CommandTool]=None) -> Optional[DataType]:
    tracer = SourceDatatypeTracer(tool)
    tracer.trace(entity)
    src_types = tracer.datatypes
    if len(src_types) == 0:
        return None
    elif len(src_types) == 1:
        return src_types[0]
    elif len(src_types) > 1:
        # remove types from Stdout as these are just janis guesses
        src_types = [x for x in src_types if not isinstance(x, Stdout)]
        # return first remaining type
        return src_types[0]
    return 

def trace_source_scatter(entity: Any, tool: Optional[CommandTool]=None) -> bool:
    tracer = SourceScatterTracer(tool)
    tracer.trace(entity)
    return tracer.source_scatter

def trace_referenced_variables(entity: Any, tool: Optional[CommandTool]) -> set[str]:
    tracer = ReferencedVariableTracer(tool)
    
    if isinstance(entity, ToolInput):
        # value, prefix, datatype
        items_to_trace = [entity.value, entity.prefix, entity.input_type]
        items_to_trace = [x for x in items_to_trace if x is not None]
        for item in items_to_trace:
            tracer.trace(item)
    
    elif isinstance(entity, ToolArgument):
        tracer.trace(entity.value)
    
    elif isinstance(entity, ToolOutput):
        tracer.trace(entity.selector)
    
    else:
        raise NotImplementedError

    return tracer.variables




class Tracer(ABC):
    
    def __init__(self, tool: Optional[CommandTool]=None):
        self.tool = tool
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
            # primitives
            list: self.trace_list,

            # datatypes
            Filename: self.filename,

            # selectors
            AliasSelector: self.alias_selector,
            InputNodeSelector: self.input_node_selector,
            WildcardSelector: self.wildcard_selector,
            InputSelector: self.input_selector,
            StepOutputSelector: self.step_output_selector,
            MemorySelector: self.resource_selector,
            CpuSelector: self.resource_selector,
            DiskSelector: self.resource_selector,
            TimeSelector: self.resource_selector,

            # misc 
            StepTagInput: self.step_tag_input,
            Edge: self.edge,
            StringFormatter: self.string_formatter,
            # InputNode: self.input_node,
        }

    @abstractmethod
    def trace(self, entity: Any) -> None:
        ...

    def operator_single_arg_trace(self, entity: Operator) -> None:
        self.trace(entity.args[0])
    
    def operator_multi_arg_trace(self, entity: Operator) -> None:
        for arg in entity.args:
            self.trace(arg)

    def trace_list(self, entity: list[Any]) -> None:
        for item in entity:
            self.trace(item)
    
    # def tool_input(self, entity: ToolInput) -> None:
    #     raise NotImplementedError
    
    # def tool_output(self, entity: ToolOutput) -> None:
    #     raise NotImplementedError
    
    # def tool_argument(self, entity: ToolArgument) -> None:
    #     raise NotImplementedError

    def alias_selector(self, entity: AliasSelector) -> None:
        self.trace(entity.inner_selector)
        self.trace(entity.data_type)

    def input_node_selector(self, entity: InputNodeSelector) -> None:
        self.trace(entity.input_node)

    def wildcard_selector(self, entity: WildcardSelector) -> None:
        # TODO check
        self.trace(entity.wildcard)

    def input_selector(self, entity: InputSelector) -> None:
        # a tool input
        assert(self.tool)
        tinput = self.tool.inputs_map()[entity.input_to_select]
        self.trace(tinput)
        
    def tool_input(self, entity: ToolInput | TInput) -> None:
        # the toolinput name
        self.trace(entity.id())
        dtype = entity.input_type if isinstance(entity, ToolInput) else entity.intype
        # the datatype if it is a Filename type
        if isinstance(dtype, Filename):
            self.trace(dtype)

    def step_output_selector(self, entity: StepOutputSelector) -> None:
        # TODO check
        self.trace(entity.node)
        self.trace(entity.tag)

    def resource_selector(self, entity: ResourceSelector) -> None:
        if hasattr(entity, 'resource_to_select'):
            self.trace(entity.resource_to_select)
        if hasattr(entity, 'resource_type'):
            self.trace(entity.resource_type)
        if hasattr(entity, 'default'):
            self.trace(entity.default)

    # def cpu_selector(self, entity: CpuSelector) -> None:
    #     if hasattr(entity, 'resource_to_select'):
    #         self.trace(entity.resource_to_select)
    #     if hasattr(entity, 'resource_type'):
    #         self.trace(entity.resource_type)
    #     if hasattr(entity, 'default'):
    #         self.trace(entity.default)

    # def disk_selector(self, entity: DiskSelector) -> None:
    #     if hasattr(entity, 'resource_to_select'):
    #         self.trace(entity.resource_to_select)
    #     if hasattr(entity, 'resource_type'):
    #         self.trace(entity.resource_type)
    #     if hasattr(entity, 'default'):
    #         self.trace(entity.default)

    # def time_selector(self, entity: TimeSelector) -> None:
    #     if hasattr(entity, 'resource_to_select'):
    #         self.trace(entity.resource_to_select)
    #     if hasattr(entity, 'resource_type'):
    #         self.trace(entity.resource_type)
    #     if hasattr(entity, 'default'):
    #         self.trace(entity.default)

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


 

class EntityCountTracer(Tracer):

    def __init__(self, tool: Optional[CommandTool]=None):
        super().__init__(tool)
        self.counter: dict[str, int] = defaultdict(int)
    
    def trace(self, entity: Any) -> None:
        etype = type(entity)
        ename = entity.__class__.__name__
        self.counter[ename] += 1

        if etype in self.custom_trace_funcs:
            func = self.custom_trace_funcs[etype]
            func(entity)
        
        elif etype in self.single_arg_trace_types:
            self.operator_single_arg_trace(entity)
        
        elif etype in self.multi_arg_trace_types:
            self.operator_multi_arg_trace(entity)
        
        else:
            pass




class SourceDatatypeTracer(Tracer):

    def __init__(self, tool: Optional[CommandTool]=None):
        super().__init__(tool)
        self.datatypes: list[DataType] = []
    
    def trace(self, entity: Any) -> None:
        # reached a leaf node (a data source)
        if self.have_reached_source(entity):
            self.handle_source(entity)

        # edge case: IndexOperator, where the target is the source
        elif isinstance(entity, IndexOperator):
            target: Any = entity.args[0]  # type: ignore
            if self.have_reached_source(target):
                self.handle_source(target, array_to_single=True)

        # other nodes: continue tracing
        else:
            etype = type(entity)
            if etype in self.custom_trace_funcs:
                func = self.custom_trace_funcs[etype]
                func(entity)
            
            elif etype in self.single_arg_trace_types:
                self.operator_single_arg_trace(entity)
            
            elif etype in self.multi_arg_trace_types:
                self.operator_multi_arg_trace(entity)
            
            else:
                pass

    def have_reached_source(self, entity: Any) -> bool:
        if isinstance(entity, InputNodeSelector) or isinstance(entity, StepOutputSelector):
            return True
        return False

    def handle_source(self, entity: InputNodeSelector | StepOutputSelector, array_to_single: bool=False) -> None:
        # reached an InputNodeSelector: get source datatype from InputNode
        if isinstance(entity, InputNodeSelector):
            dtype = entity.input_node.datatype
        
        # reached a StepOutputSelector: get source datatype from ToolOutput
        elif isinstance(entity, StepOutputSelector):
            step = entity.node
            outtag = entity.tag
            tout = step.tool.outputs_map()[outtag]
            dtype = tout.outtype
        
        else:
            raise RuntimeError
        
        # if the datatype is an Array(File()), but we are selecting a specific index
        # (via IndexOperator), the datatype will actually be File()
        assert(dtype)
        if array_to_single and dtype.is_array():
            dtype = nfgen_utils.get_base_type(dtype)
        
        self.datatypes.append(dtype)




class SourceScatterTracer(Tracer):
    
    def __init__(self, tool: Optional[CommandTool]=None):
        super().__init__(tool)
        self.source_scatter: bool = False

    def trace(self, entity: Any) -> None:
        # only steps have scatter, so only StepOutputSelector needs to be checked. 
        if isinstance(entity, StepOutputSelector):
            self.handle_step_output_selector(entity)

        # other nodes: continue tracing
        else:
            etype = type(entity)
            if etype in self.custom_trace_funcs:
                func = self.custom_trace_funcs[etype]
                func(entity)
            
            elif etype in self.single_arg_trace_types:
                self.operator_single_arg_trace(entity)
            
            elif etype in self.multi_arg_trace_types:
                self.operator_multi_arg_trace(entity)
            
            else:
                pass

    def handle_step_output_selector(self, entity: StepOutputSelector) -> None:
        # we have found the step which feeds the data.
        # check if the step is scattered. 
        step = entity.node
        if step.scatter:
            self.source_scatter = True
        


class ReferencedVariableTracer(Tracer):

    def __init__(self, tool: Optional[CommandTool]=None):
        super().__init__(tool)
        self.variables: set[str] = set()
    
    def trace(self, entity: Any) -> None:
        # reached a leaf node (variable reference to tool input)
        if isinstance(entity, (ToolInput, TInput)):
            # we need to check whether the TInput is a deadend.
            # in the case of a Filename type, can always guarantee a value? I hope? 
            dtype = entity.input_type if isinstance(entity, ToolInput) else entity.intype
            if isinstance(dtype, Filename):
                pass
            else:
                self.variables.add(entity.id())
        
        # other nodes: continue tracing
        else:
            etype = type(entity)
            if etype in self.custom_trace_funcs:
                func = self.custom_trace_funcs[etype]
                func(entity)
            
            elif etype in self.single_arg_trace_types:
                self.operator_single_arg_trace(entity)
            
            elif etype in self.multi_arg_trace_types:
                self.operator_multi_arg_trace(entity)
            
            else:
                # for things we dont care about. 
                # dead paths which don't lead anywhere. 
                pass

