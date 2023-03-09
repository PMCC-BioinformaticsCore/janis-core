
from typing import Any, Optional

from janis_core import ToolInput, ToolArgument, ToolOutput, WildcardSelector, CommandToolBuilder, CommandTool
from janis_core import settings

from ..types import ingest_cwl_type
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

        # requirements
        req_parser = RequirementsParser(self.cwl_utils)
        requirements = req_parser.parse(entity, jtool.uuid)
        
        jtool.files_to_create = requirements['files_to_create'] or None
        jtool.env_vars = requirements['env_vars'] or None
        jtool.container = requirements['container']
        jtool.memory = requirements['memory']
        jtool.cpus = requirements['cpus']
        jtool.time = requirements['time']

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
            res, success = parse_basic_expression(entity.valueFrom)
            if not success:
                self.error_msgs.append('untranslated javascript expression')
            arg = ToolArgument(
                value=res,
                position=entity.position,
                prefix=entity.prefix,
                separate_value_from_prefix=entity.separate,
                shell_quote=entity.shellQuote,
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
                expr = inpBinding.valueFrom
                value, success = parse_basic_expression(expr)
                if not success:
                    msg = f'untranslated javascript expression in inputBinding.valueFrom'
                    self.error_msgs.append(msg)
            else:
                pass
                # j.Logger.warn(
                #     f"Won't translate the expression for input {entity.id}: {inpBinding.valueFrom}"
                # )

        dtype, error_msgs = ingest_cwl_type(entity.type, self.cwl_utils, secondary_files=entity.secondaryFiles)
        self.error_msgs += error_msgs
        
        tinput = ToolInput(
            tag=identifier,
            input_type=dtype,
            position=inpBinding.position if inpBinding else None,
            prefix=inpBinding.prefix if inpBinding else None,
            separate_value_from_prefix=inpBinding.separate if inpBinding else None,
            separator=inpBinding.itemSeparator if inpBinding else None,
            shell_quote=inpBinding.shellQuote if inpBinding else None,
            default=entity.default,
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

