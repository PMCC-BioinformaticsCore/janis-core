
import copy
from typing import Any, Optional

from janis_core import ToolInput, ToolArgument, ToolOutput, WildcardSelector, CommandToolBuilder, CommandTool
from janis_core import settings

from ..types import ingest_cwl_type
from ..types import cast_cwl_type_to_python
from ..identifiers import get_id_entity
from ..identifiers import get_id_filename
from ..expressions import parse_basic_expression
from .common import EntityParser
from .common import RequirementsParser

class CLTToolParser(EntityParser):

    def parse(self, entity: Any) -> CommandTool:
        identifier = get_id_filename(entity.id)

        inputs = [self.ingest_command_tool_input(inp) for inp in entity.inputs]
        outputs = [self.ingest_command_tool_output(out) for out in entity.outputs]
        arguments = [self.ingest_command_tool_argument(arg) for arg in (entity.arguments or [])]

        jtool = CommandToolBuilder(
            tool=identifier,
            base_command=entity.baseCommand,
            inputs=inputs,
            outputs=outputs,
            arguments=arguments,
            version="v0.1.0",
            doc=entity.doc,
            friendly_name=entity.label,
            container="ubuntu:latest"
        )

        extra_args = self.ingest_io_streams(entity, jtool)
        # addressing cwl stdin: stdout: stderr: file naming
        jtool._arguments += extra_args

        # requirements
        req_parser = RequirementsParser(self.cwl_utils)
        requirements = req_parser.parse(entity, jtool.uuid)
        
        jtool._files_to_create = requirements['files_to_create'] or None
        jtool._env_vars = requirements['env_vars'] or None
        jtool._container = requirements['container']
        jtool._memory = requirements['memory']
        jtool._cpus = requirements['cpus']
        jtool._time = requirements['time']

        self.log_messages(jtool.uuid)
        return jtool
    
    def ingest_command_tool_argument(self, arg: Any) -> ToolArgument:
        parser = CLTArgumentParser(self.cwl_utils)
        return parser.parse(arg)

    def ingest_command_tool_input(self, inp: Any) -> ToolInput:
        parser = CLTInputParser(self.cwl_utils)
        return parser.parse(inp)
        
    def ingest_command_tool_output(self, out: Any, is_expression_tool: bool=False) -> ToolOutput:  
        parser = CLTOutputParser(self.cwl_utils)
        return parser.parse(out)
    
    def ingest_io_streams(self, entity: Any, jtool: CommandTool) -> list[ToolArgument]:
        # n = last position
        # stderr = n + 1
        # stdout = n + 2
        # stdin = n + 3
        out: list[ToolArgument] = []
        n = self.get_last_input_position(jtool)
        
        if entity.stderr:
            res, success = parse_basic_expression(entity.stderr)
            if not success:
                self.error_msgs.append('untranslated javascript expression in stderr filename')
            arg = ToolArgument(prefix='2>', value=res, position=n + 1)
            out.append(arg)
        
        if entity.stdout:
            res, success = parse_basic_expression(entity.stdout)
            if not success:
                self.error_msgs.append('untranslated javascript expression in stdout filename')
            arg = ToolArgument(prefix='1>', value=res, position=n + 2)
            out.append(arg)
        
        if entity.stdin:
            res, success = parse_basic_expression(entity.stdin)
            if not success:
                self.error_msgs.append('untranslated javascript expression in stdin filename')
            arg = ToolArgument(prefix='<', value=res, position=n + 3)
            out.append(arg)

        return out

    def get_last_input_position(self, jtool: CommandTool) -> int:
        max_pos = 0
        inputs: list[ToolInput | ToolArgument] = copy.deepcopy(jtool.inputs())
        if jtool.arguments():
            inputs += copy.deepcopy(jtool.arguments())
        for inp in inputs:
            if inp.position and inp.position > max_pos:
                max_pos = inp.position
        return max_pos



class CLTArgumentParser(EntityParser):

    def parse(self, entity: Any) -> Any:
        # I don't know when a clt argument would be a string
        if isinstance(entity, str):
            res, success = parse_basic_expression(entity)
            if not success:
                self.error_msgs.append('untranslated javascript expression')
            arg = ToolArgument(res)
        
        # normal case
        else:
            valueFrom = cast_cwl_type_to_python(entity.valueFrom)
            res, success = parse_basic_expression(valueFrom)
            if not success:
                self.error_msgs.append('untranslated javascript expression')
            arg = ToolArgument(
                value=res,
                position=cast_cwl_type_to_python(entity.position),
                prefix=cast_cwl_type_to_python(entity.prefix),
                separate_value_from_prefix=cast_cwl_type_to_python(entity.separate),
                shell_quote=cast_cwl_type_to_python(entity.shellQuote),
            )
        
        self.log_messages(arg.uuid)
        return arg



class CLTInputParser(EntityParser):

    def parse(self, entity: Any) -> Any:
        identifier = get_id_entity(entity.id)

        value: Optional[str] = None
        inpBinding = entity.inputBinding
        if inpBinding and inpBinding.valueFrom:
            if settings.ingest.cwl.INGEST_JAVASCRIPT_EXPRESSIONS:
                expr = cast_cwl_type_to_python(inpBinding.valueFrom)
                value, success = parse_basic_expression(expr)
                if not success:
                    msg = f'untranslated javascript expression in inputBinding.valueFrom'
                    self.error_msgs.append(msg)
            else:
                pass
                # j.Logger.warn(
                #     f"Won't translate the expression for input {entity.id}: {inpBinding.valueFrom}"
                # )

        etype = cast_cwl_type_to_python(entity.type)
        secondary_files = cast_cwl_type_to_python(entity.secondaryFiles)
        dtype, error_msgs = ingest_cwl_type(etype, self.cwl_utils, secondary_files=secondary_files)
        self.error_msgs += error_msgs

        position = cast_cwl_type_to_python(inpBinding.position) if inpBinding else None
        prefix = cast_cwl_type_to_python(inpBinding.prefix) if inpBinding else None
        separate_value_from_prefix = cast_cwl_type_to_python(inpBinding.separate) if inpBinding else None
        separator = cast_cwl_type_to_python(inpBinding.itemSeparator) if inpBinding else None
        shell_quote = cast_cwl_type_to_python(inpBinding.shellQuote) if inpBinding else None
        default = cast_cwl_type_to_python(entity.default)
        
        tinput = ToolInput(
            tag=identifier,
            input_type=dtype,
            position=position,
            prefix=prefix,
            separate_value_from_prefix=separate_value_from_prefix,
            separator=separator,
            shell_quote=shell_quote,
            default=default,
            value=value
        )

        self.log_messages(tinput.uuid)
        return tinput



class CLTOutputParser(EntityParser):

    def parse(self, entity: Any) -> Any:
        # tag
        identifier = get_id_entity(entity.id)
        
        # datatype
        dtype, error_msgs = ingest_cwl_type(entity.type, self.cwl_utils, secondary_files=entity.secondaryFiles)
        self.error_msgs += error_msgs

        # selector
        selector = None
        if hasattr(entity, 'janis_collection_override'):
            selector = entity.janis_collection_override

        elif entity.outputBinding:
            if entity.outputBinding.glob:
                res, success = parse_basic_expression(entity.outputBinding.glob)
                selector = WildcardSelector(res)
                if not success:
                    self.error_msgs.append('untranslated javascript expression in output collection')
            
            elif entity.outputBinding.outputEval:
                res, success = parse_basic_expression(entity.outputBinding.outputEval)
                selector = res
                if not success:
                    self.error_msgs.append('untranslated javascript expression in output collection')

        toutput = ToolOutput(identifier, dtype, selector)
        self.log_messages(toutput.uuid)
        return toutput

