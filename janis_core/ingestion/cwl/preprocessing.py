

from typing import Any
from .types import cast_cwl_type_to_python


def handle_inline_cltool_identifiers(cwl_workflow: Any, cwl_utils: Any) -> Any:
    """
    alters tool ids in a cwl_workflow. 
    workflows with inline CommandLineTool definitions (rather than those in a dedicated tool.cwl file)
    will have random ids starting with _: for those inline tool definitions. 
    in these cases, this function overrides the random id, and replaces it with an id based on the step id.
    """
    for step in cwl_workflow.steps:
        # handle ids for each step task (tool / workflow)
        task = step.run
        if isinstance(task, str):  # referencing external file, id is fine, ignore
            pass
        elif task.id.startswith('_:'):
            step_name = step.id.rsplit('#', 1)[1]
            workdir = step.id.rsplit('#', 1)[0].rsplit('/', 1)[0]
            if isinstance(task, cwl_utils.CommandLineTool) or isinstance(task, cwl_utils.ExpressionTool):
                new_id = f'{workdir}/{step_name}_tool.cwl'
            elif isinstance(task, cwl_utils.Workflow):
                new_id = f'{workdir}/{step_name}_workflow.cwl'
            else:
                raise NotImplementedError
            task.id = new_id

        # if task is a subworkflow, handle identifiers in the subworkflow too
        if isinstance(task, cwl_utils.Workflow):
            handle_inline_cltool_identifiers(task, cwl_utils)
    
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
        elif isinstance(entity, self.cwl_utils.ExpressionTool):
            entity = self._convert_types_etool(entity)
        return entity

    def _convert_types_workflow(self, wf: Any) -> Any:
        wf.label = cast_cwl_type_to_python(wf.label)
        wf.doc = cast_cwl_type_to_python(wf.doc)
        
        if wf.hints:
            wf.hints = [cast_cwl_type_to_python(req) for req in wf.hints]
        
        if wf.requirements:
            wf.requirements = [cast_cwl_type_to_python(req) for req in wf.requirements]

        for inp in wf.inputs:
            self._convert_types_workflow_input(inp)
        for stp in wf.steps:
            self._convert_types_workflow_step(stp)
        for out in wf.outputs:
            self._convert_types_workflow_output(out)
        return wf

    def _convert_types_workflow_input(self, inp: Any) -> None:
        inp.type = cast_cwl_type_to_python(inp.type)
        inp.secondaryFiles = cast_cwl_type_to_python(inp.secondaryFiles)
        inp.format = cast_cwl_type_to_python(inp.format)
        inp.default = cast_cwl_type_to_python(inp.default)
        inp.label = cast_cwl_type_to_python(inp.label)
        inp.doc = cast_cwl_type_to_python(inp.doc)

        if hasattr(inp, 'inputBinding') and inp.inputBinding:
            self._convert_types_inputbinding(inp)

    def _convert_types_inputbinding(self, binding: Any) -> None:
        binding.inputBinding.position = cast_cwl_type_to_python(binding.inputBinding.position)
        binding.inputBinding.prefix = cast_cwl_type_to_python(binding.inputBinding.prefix)
        binding.inputBinding.separate = cast_cwl_type_to_python(binding.inputBinding.separate)
        binding.inputBinding.itemSeparator = cast_cwl_type_to_python(binding.inputBinding.itemSeparator)
        binding.inputBinding.shellQuote = cast_cwl_type_to_python(binding.inputBinding.shellQuote)
        binding.inputBinding.valueFrom = cast_cwl_type_to_python(binding.inputBinding.valueFrom)

    def _convert_types_workflow_step(self, step: Any) -> None:
        # step 
        if step.scatter:
            print()
        step.doc = cast_cwl_type_to_python(step.doc)
        step.scatter = cast_cwl_type_to_python(step.scatter)
        step.scatterMethod = cast_cwl_type_to_python(step.scatterMethod)
        # run?
        # entity.run = cast_cwl_type_to_python(entity.run)

        # step inputs
        for inp in step.in_:
            inp.default = cast_cwl_type_to_python(inp.default)
            inp.valueFrom = cast_cwl_type_to_python(inp.valueFrom)

    def _convert_types_workflow_output(self, out: Any) -> None:
        out.doc = cast_cwl_type_to_python(out.doc)
        out.format = cast_cwl_type_to_python(out.format)
        out.label = cast_cwl_type_to_python(out.label)
        out.type = cast_cwl_type_to_python(out.type)
        out.secondaryFiles = cast_cwl_type_to_python(out.secondaryFiles)
        out.outputSource = cast_cwl_type_to_python(out.outputSource)
        
        if hasattr(out, 'outputBinding') and out.outputBinding:
            self._convert_types_outputbinding(out)

    def _convert_types_outputbinding(self, binding: Any) -> None:
        binding.outputBinding.glob = cast_cwl_type_to_python(binding.outputBinding.glob)
        binding.outputBinding.outputEval = cast_cwl_type_to_python(binding.outputBinding.outputEval)

    def _convert_types_clt(self, clt: Any) -> Any:
        clt.baseCommand = cast_cwl_type_to_python(clt.baseCommand)
        clt.stderr = cast_cwl_type_to_python(clt.stderr)
        clt.stdout = cast_cwl_type_to_python(clt.stdout)
        clt.stdin = cast_cwl_type_to_python(clt.stdin)
        clt.label = cast_cwl_type_to_python(clt.label)
        clt.doc = cast_cwl_type_to_python(clt.doc)
        
        if clt.hints:
            clt.hints = [cast_cwl_type_to_python(req) for req in clt.hints]
        
        if clt.requirements:
            clt.requirements = self._convert_types_clt_requirements(clt.requirements)

        if clt.arguments:
            for i, arg in enumerate(clt.arguments):
                if isinstance(arg, self.cwl_utils.CommandLineBinding):
                    self._convert_types_clt_argument(arg)
                else:
                    clt.arguments[i] = cast_cwl_type_to_python(arg)
                    
        for inp in clt.inputs:
            self._convert_types_clt_input(inp)
        for out in clt.outputs:
            self._convert_types_clt_output(out)
        
        return clt

    def _convert_types_clt_requirements(self, reqs: Any) -> None:
        for req in reqs:
            if isinstance(req, self.cwl_utils.ResourceRequirement):
                req.coresMax = cast_cwl_type_to_python(req.coresMax)
                req.coresMin = cast_cwl_type_to_python(req.coresMin)
                req.outdirMax = cast_cwl_type_to_python(req.outdirMax)
                req.outdirMin = cast_cwl_type_to_python(req.outdirMin)
                req.ramMax = cast_cwl_type_to_python(req.ramMax)
                req.ramMin = cast_cwl_type_to_python(req.ramMin)
                req.tmpdirMax = cast_cwl_type_to_python(req.tmpdirMax)
                req.tmpdirMin = cast_cwl_type_to_python(req.tmpdirMin)
            
            elif isinstance(req, self.cwl_utils.DockerRequirement):
                req.class_ = cast_cwl_type_to_python(req.class_)
                req.dockerFile = cast_cwl_type_to_python(req.dockerFile)
                req.dockerImageId = cast_cwl_type_to_python(req.dockerImageId)
                req.dockerImport = cast_cwl_type_to_python(req.dockerImport)
                req.dockerLoad = cast_cwl_type_to_python(req.dockerLoad)
                req.dockerOutputDirectory = cast_cwl_type_to_python(req.dockerOutputDirectory)
                req.dockerPull = cast_cwl_type_to_python(req.dockerPull)
            
            elif hasattr(req, 'timelimit') and isinstance(req, self.cwl_utils.ToolTimeLimit):
                req.timelimit = cast_cwl_type_to_python(req.timelimit)

            elif isinstance(req, self.cwl_utils.InitialWorkDirRequirement):
                if isinstance(req.listing, str):
                    req.listing = cast_cwl_type_to_python(req.listing)
                elif isinstance(req.listing, list):
                    for i, item in enumerate(req.listing):
                        if isinstance(item, str):
                            req.listing[i] = cast_cwl_type_to_python(item)
                        elif isinstance(item, self.cwl_utils.Dirent):
                            dirent: self.cwl_utils.Dirent = item
                            dirent.entry = cast_cwl_type_to_python(dirent.entry)
                            dirent.entryname = cast_cwl_type_to_python(dirent.entryname)
                            dirent.writable = cast_cwl_type_to_python(dirent.writable)
                        else:
                            raise NotImplementedError
                else:
                    raise NotImplementedError

            elif isinstance(req, self.cwl_utils.EnvVarRequirement):
                for envdef in req.envDef:
                    envdef.envName = cast_cwl_type_to_python(envdef.envName)
                    envdef.envValue = cast_cwl_type_to_python(envdef.envValue)
        return reqs
        
    def _convert_types_clt_argument(self, arg: Any) -> None:
        arg.position = cast_cwl_type_to_python(arg.position)
        arg.valueFrom = cast_cwl_type_to_python(arg.valueFrom)
        arg.prefix = cast_cwl_type_to_python(arg.prefix)
        arg.separate = cast_cwl_type_to_python(arg.separate)
        arg.shellQuote = cast_cwl_type_to_python(arg.shellQuote)

    def _convert_types_clt_input(self, inp: Any) -> None:
        inp.type = cast_cwl_type_to_python(inp.type)
        inp.secondaryFiles = cast_cwl_type_to_python(inp.secondaryFiles)
        inp.default = cast_cwl_type_to_python(inp.default)
        inp.doc = cast_cwl_type_to_python(inp.doc)
        
        if hasattr(inp, 'inputBinding') and inp.inputBinding:
            self._convert_types_inputbinding(inp)

    def _convert_types_clt_output(self, out: Any) -> None:
        out.type = cast_cwl_type_to_python(out.type)
        out.secondaryFiles = cast_cwl_type_to_python(out.secondaryFiles)
        
        if hasattr(out, 'outputBinding') and out.outputBinding:
            self._convert_types_outputbinding(out)

    def _convert_types_etool(self, etool: Any) -> Any:
        etool.expression = cast_cwl_type_to_python(etool.expression)
        etool.label = cast_cwl_type_to_python(etool.label)
        etool.doc = cast_cwl_type_to_python(etool.doc)
        
        if etool.hints:
            etool.hints = [cast_cwl_type_to_python(req) for req in etool.hints]
        
        if etool.requirements:
            etool.requirements = self._convert_types_clt_requirements(etool.requirements)

        for inp in etool.inputs:
            self._convert_types_clt_input(inp)
        for out in etool.outputs:
            self._convert_types_clt_output(out)
        
        return etool
