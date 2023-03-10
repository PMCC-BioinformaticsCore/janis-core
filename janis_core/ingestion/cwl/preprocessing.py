

from typing import Any
from .types import cast_cwl_type_to_python


def handle_inline_cltool_identifiers(cwl_workflow: Any) -> Any:
    """
    alters tool ids in a cwl_workflow. 
    workflows with inline CommandLineTool definitions (rather than those in a dedicated tool.cwl file)
    will have random ids starting with _: for those inline tool definitions. 
    in these cases, this function overrides the random id, and replaces it with an id based on the step id.
    """
    for step in cwl_workflow.steps:
        tool = step.run
        if isinstance(tool, str):
            pass
        elif tool.id.startswith('_:'):
            step_name = step.id.rsplit('#', 1)[1]
            workdir = step.id.rsplit('#', 1)[0].rsplit('/', 1)[0]
            new_id = f'{workdir}/{step_name}_tool.cwl'
            tool.id = new_id
    return cwl_workflow



"""
This could be more thorough. Should check that each field in 
https://www.commonwl.org/v1.2/CommandLineTool.html and 
https://www.commonwl.org/v1.2/Workflow.html
is being converted. 
"""
def convert_cwl_types_to_python(entity: Any, cwl_utils: Any) -> Any:
    converter = CWLTypeConverter(cwl_utils)
    return converter.convert(entity)

class CWLTypeConverter:
    def __init__(self, cwl_utils: Any) -> None:
        self.cwl_utils = cwl_utils

    def convert(self, entity: Any) -> Any:       
        if isinstance(entity, self.cwl_utils.Workflow):
            entity = self._convert_types_workflow(entity)
        elif isinstance(entity, self.cwl_utils.CommandLineTool):
            entity = self._convert_types_clt(entity)
        else:
            raise RuntimeError
        return entity

    def _convert_types_workflow(self, entity: Any) -> Any:
        entity.label = cast_cwl_type_to_python(entity.label)
        entity.doc = cast_cwl_type_to_python(entity.doc)
        
        if entity.hints:
            entity.hints = [cast_cwl_type_to_python(req) for req in entity.hints]
        
        if entity.requirements:
            entity.requirements = [cast_cwl_type_to_python(req) for req in entity.requirements]

        for inp in entity.inputs:
            self._convert_types_workflow_input(inp)
        for stp in entity.steps:
            self._convert_types_workflow_step(stp)
        for out in entity.outputs:
            self._convert_types_workflow_output(out)
        return entity

    def _convert_types_workflow_input(self, entity: Any) -> None:
        entity.type = cast_cwl_type_to_python(entity.type)
        entity.secondaryFiles = cast_cwl_type_to_python(entity.secondaryFiles)
        entity.format = cast_cwl_type_to_python(entity.format)
        entity.default = cast_cwl_type_to_python(entity.default)
        entity.label = cast_cwl_type_to_python(entity.label)
        entity.doc = cast_cwl_type_to_python(entity.doc)

        if hasattr(entity, 'inputBinding') and entity.inputBinding:
            self._convert_types_inputbinding(entity)

    def _convert_types_inputbinding(self, entity: Any) -> None:
        entity.inputBinding.position = cast_cwl_type_to_python(entity.inputBinding.position)
        entity.inputBinding.prefix = cast_cwl_type_to_python(entity.inputBinding.prefix)
        entity.inputBinding.separate = cast_cwl_type_to_python(entity.inputBinding.separate)
        entity.inputBinding.itemSeparator = cast_cwl_type_to_python(entity.inputBinding.itemSeparator)
        entity.inputBinding.shellQuote = cast_cwl_type_to_python(entity.inputBinding.shellQuote)
        entity.inputBinding.valueFrom = cast_cwl_type_to_python(entity.inputBinding.valueFrom)

    def _convert_types_workflow_step(self, entity: Any) -> None:
        # step 
        if entity.scatter:
            print()
        entity.doc = cast_cwl_type_to_python(entity.doc)
        entity.scatter = cast_cwl_type_to_python(entity.scatter)
        entity.scatterMethod = cast_cwl_type_to_python(entity.scatterMethod)
        # run?
        # entity.run = cast_cwl_type_to_python(entity.run)

        # step inputs
        for inp in entity.in_:
            inp.default = cast_cwl_type_to_python(inp.default)
            inp.valueFrom = cast_cwl_type_to_python(inp.valueFrom)

    def _convert_types_workflow_output(self, entity: Any) -> None:
        entity.doc = cast_cwl_type_to_python(entity.doc)
        entity.format = cast_cwl_type_to_python(entity.format)
        entity.label = cast_cwl_type_to_python(entity.label)
        entity.type = cast_cwl_type_to_python(entity.type)
        entity.secondaryFiles = cast_cwl_type_to_python(entity.secondaryFiles)
        entity.outputSource = cast_cwl_type_to_python(entity.outputSource)
        
        if hasattr(entity, 'outputBinding') and entity.outputBinding:
            self._convert_types_outputbinding(entity)

    def _convert_types_outputbinding(self, entity: Any) -> None:
        entity.outputBinding.glob = cast_cwl_type_to_python(entity.outputBinding.glob)
        entity.outputBinding.outputEval = cast_cwl_type_to_python(entity.outputBinding.outputEval)

    def _convert_types_clt(self, entity: Any) -> Any:
        entity.baseCommand = cast_cwl_type_to_python(entity.baseCommand)
        entity.stderr = cast_cwl_type_to_python(entity.stderr)
        entity.stdout = cast_cwl_type_to_python(entity.stdout)
        entity.stdin = cast_cwl_type_to_python(entity.stdin)
        entity.label = cast_cwl_type_to_python(entity.label)
        entity.doc = cast_cwl_type_to_python(entity.doc)
        
        if entity.hints:
            entity.hints = [cast_cwl_type_to_python(req) for req in entity.hints]
        
        if entity.requirements:
            entity.requirements = [cast_cwl_type_to_python(req) for req in entity.requirements]

        if entity.arguments:
            for i, arg in enumerate(entity.arguments):
                if isinstance(arg, self.cwl_utils.CommandLineBinding):
                    self._convert_types_clt_argument(arg)
                else:
                    entity.arguments[i] = cast_cwl_type_to_python(arg)
                    
        for inp in entity.inputs:
            self._convert_types_clt_input(inp)
        for out in entity.outputs:
            self._convert_types_clt_output(out)
        
        return entity

    def _convert_types_clt_argument(self, entity: Any) -> None:
        entity.position = cast_cwl_type_to_python(entity.position)
        entity.valueFrom = cast_cwl_type_to_python(entity.valueFrom)
        entity.prefix = cast_cwl_type_to_python(entity.prefix)
        entity.separate = cast_cwl_type_to_python(entity.separate)
        entity.shellQuote = cast_cwl_type_to_python(entity.shellQuote)

    def _convert_types_clt_input(self, entity: Any) -> None:
        entity.type = cast_cwl_type_to_python(entity.type)
        entity.secondaryFiles = cast_cwl_type_to_python(entity.secondaryFiles)
        entity.default = cast_cwl_type_to_python(entity.default)
        entity.doc = cast_cwl_type_to_python(entity.doc)
        
        if hasattr(entity, 'inputBinding') and entity.inputBinding:
            self._convert_types_inputbinding(entity)

    def _convert_types_clt_output(self, entity: Any) -> None:
        entity.type = cast_cwl_type_to_python(entity.type)
        entity.secondaryFiles = cast_cwl_type_to_python(entity.secondaryFiles)
        
        if hasattr(entity, 'outputBinding') and entity.outputBinding:
            self._convert_types_outputbinding(entity)