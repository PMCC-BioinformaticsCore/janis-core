
from typing import Optional
from copy import deepcopy

from janis_core import (
    InputSelector,
    WildcardSelector
)
from janis_core import (
    DataType,
    Stdout,
    Array,
    String,
    Int, 
    Float,
    Boolean,
    File,
    Directory
)
from janis_core.ingestion.galaxy.gx.gxtool.param import OutputParam, CollectionOutputParam, DataOutputParam

from janis_core.ingestion.galaxy import tags
from janis_core.ingestion.galaxy.model.workflow.input import WorkflowInput
from janis_core.ingestion.galaxy.model.workflow.step.outputs import StepOutput
from janis_core.ingestion.galaxy.gx.command.components import (
    InputComponent,
    OutputComponent,
    RedirectOutput,
    InputOutput,
    WildcardOutput,
)

datatype_map = {
    'String': String(),
    'Int': Int(),
    'Float': Float(),
    'Boolean': Boolean(),
    'File': File(),
    'Directory': Directory(),
}

DATATYPE_COMPONENT = InputComponent | OutputComponent | WorkflowInput | StepOutput

def to_janis_datatype(component: DATATYPE_COMPONENT) -> DataType:
    if isinstance(component, StepOutput):
        component = component.tool_output
    
    # the underlying datatype, regardless of stdout / array / optionality
    dtype = deepcopy(datatype_map[component.datatype.classname])
    
    # modifying dtype in array case
    if component.array:
        dtype = Array(dtype)
    
    # modifying dtype in stdout case
    if isinstance(component, RedirectOutput):
        dtype = Stdout(dtype)

    # modifying dtype in optional case
    if component.optional:
        dtype.optional = True
    return dtype


def to_janis_selector(component: OutputComponent) -> Optional[InputSelector | WildcardSelector]:
    # # redirect outputs don't get a selector as they are Stdout type
    if isinstance(component, RedirectOutput):
        return None
    
    # if isinstance(component, RedirectOutput):
    #     input_comp_tag = tags.get(component.input_component.uuid)
    #     return InputSelector(input_comp_tag)
    
    pattern = get_collection_pattern(component)
    if pattern is not None:
        return WildcardSelector(pattern)
    
    elif isinstance(component, InputOutput):
        input_comp_tag = tags.get(component.input_component.uuid)
        return InputSelector(input_comp_tag)
    
    else:
        return WildcardSelector('unknown_collection_pattern')

def get_collection_pattern(component: OutputComponent) -> Optional[str]:
    # linked to CollectionOutputParam
    if isinstance(component.gxparam, CollectionOutputParam):
        if component.gxparam.discover_pattern:
            return component.gxparam.discover_pattern
    
    # linked to DataOutputParam
    if isinstance(component.gxparam, DataOutputParam):
        if component.gxparam.from_work_dir:
            return component.gxparam.from_work_dir
        if component.gxparam.discover_pattern:
            return component.gxparam.discover_pattern
    
    # InputOutput, but the actual galaxy param being linked is something in the tool
    # <outputs> section, rather than in the <inputs> section
    if isinstance(component, InputOutput):
        
        # linked to CollectionOutputParam
        if isinstance(component.gxparam, OutputParam):
            if isinstance(component.gxparam, CollectionOutputParam):
                if component.gxparam.discover_pattern:
                    return component.gxparam.discover_pattern
            
            # linked to DataOutputParam
            if isinstance(component.gxparam, DataOutputParam):
                if component.gxparam.from_work_dir:
                    return component.gxparam.from_work_dir
                if component.gxparam.discover_pattern:
                    return component.gxparam.discover_pattern
            
            # galaxy <output> param doesn't have a collector, so galaxy autotemplates it. 
            # can just provide any value, but will provide something that make sense 
            # using the <output> param's name & the first ext if exists
            if component.gxparam.formats:
                return f'{component.gxparam.name}.{component.gxparam.formats[0]}'
            else:
                return component.gxparam.name

    # can't get a wildcard pattern    
    return None

        
    
    


