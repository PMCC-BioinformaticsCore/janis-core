
import copy
from typing import Any, Optional
from dataclasses import dataclass, field
from abc import ABC, abstractmethod

from janis_core import ToolInput, ToolArgument, ToolOutput, WildcardSelector, CommandToolBuilder, CommandTool, InputSelector
from janis_core.types import File, Stdout, Stderr
from janis_core import settings
from janis_core.messages import log_warning

from ..types import ingest_cwl_type
from ..identifiers import get_id_entity
from ..identifiers import get_id_filename
from ..expressions import parse_basic_expression
from ..preprocessing import convert_cwl_types_to_python



@dataclass
class CLTEntityParser(ABC):
    cwl_utils: Any
    entity: Any
    is_expression_tool: bool = False
    success: bool = False
    uuid: Optional[str] = None
    error_msgs: list[str] = field(default_factory=list)
    failure_message: str = ''

    def parse(self) -> Any:
        # normal mode
        if settings.ingest.SAFE_MODE:
            try:
                j_entity = self.do_parse()
                self.success = True
            except Exception:
                j_entity = self.fallback()
        
        # dev mode
        else:
            j_entity = self.do_parse()
            self.success = True
        
        self.log_messages()
        return j_entity

    def log_messages(self) -> None:
        if self.success == False:
            self.error_msgs.append(self.failure_message)
        
        for msg in self.error_msgs:
            log_warning(self.uuid, msg)

    @abstractmethod
    def do_parse(self) -> Any:
        ...
    
    @abstractmethod
    def fallback(self) -> Any:
        ...



@dataclass
class CLTRequirementsParser(CLTEntityParser):

    def fallback(self) -> dict[str, Any]:
        return {}

    def do_parse(self) -> dict[str, Any]:  # type: ignore I FKN KNOW
        requirements = self.entity.requirements or []

        out: dict[str, Any] = {
            'container': None,
            'memory': None,
            'cpus': None,
            'time': None,
            'directories_to_create': {},
            'files_to_create': {},
            'env_vars': {}
        }
        
        out['container'] = self.get_container(requirements)
        out['memory'] = self.get_memory(requirements)
        out['cpus'] = self.get_cpus(requirements)
        out['time'] = self.get_time(requirements)
        out['directories_to_create'] = self.get_directories_to_create(requirements)
        out['files_to_create'] = self.get_files_to_create(requirements)
        out['env_vars'] = self.get_env_vars(requirements)
        
        return out

    def get_container(self, requirements: list[Any]) -> str:
        fallback = 'ubuntu:latest'
        if self.is_expression_tool:
            return 'node:latest'
        else:
            for req in requirements:
                if isinstance(req, self.cwl_utils.DockerRequirement):
                    return str(req.dockerPull)
        return fallback
    
    def get_memory(self, requirements: list[Any]) -> Optional[int]:
        for req in requirements:
            if isinstance(req, self.cwl_utils.ResourceRequirement):
                # maybe convert mebibytes to megabytes?
                memory, success = parse_basic_expression(req.ramMin or req.ramMax)
                if not success:
                    msg = 'untranslated javascript expression in task MEM requirement'
                    self.error_msgs.append(msg)
                return memory
        return None
    
    def get_cpus(self, requirements: list[Any]) -> Optional[int]:
        for req in requirements:
            if isinstance(req, self.cwl_utils.ResourceRequirement):
                cpus, success = parse_basic_expression(req.coresMin)
                if not success:
                    msg = 'untranslated javascript expression in task CPU requirement'
                    self.error_msgs.append(msg)
                return cpus
        return None

    def get_time(self, requirements: list[Any]) -> Optional[int]:
        for req in requirements:
            if hasattr(req, 'timelimit') and isinstance(req, self.cwl_utils.ToolTimeLimit):
                time, success = parse_basic_expression(req.timelimit)
                if not success:
                    msg = 'untranslated javascript expression in task TIME requirement'
                    self.error_msgs.append(msg)
                return time
        return None
    
    def get_files_to_create(self, requirements: list[Any]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        
        for req in requirements:
            if isinstance(req, self.cwl_utils.InitialWorkDirRequirement):
                "array<File | Directory | Dirent | string | Expression> | string | Expression"
                if isinstance(req.listing, str) and req.listing.startswith('$('):
                    # Expression
                    raise NotImplementedError
                
                elif isinstance(req.listing, str):
                    # string
                    raise NotImplementedError
                
                else:
                    # array<File | Directory | Dirent | string | Expression>
                    for item in req.listing:
                        # Expression
                        if isinstance(item, str):
                            if item.startswith('$('):
                                # most likely just specifying the file / directory needs to be staged. do nothing for now.
                                # future: detect whether its a file or directory. create directory if directory. 
                                pass  
                            
                            elif item.startswith('${'):
                                # most likely a js expr which returns a file or directory to stage. 
                                # future: detect whether its a file or directory. create directory if directory. 
                                pass  
                            
                            else:
                                # no idea yet. 
                                raise NotImplementedError
                        
                        # Dirent
                        elif isinstance(item, self.cwl_utils.Dirent):
                            dirent = item
                            if dirent.entry.startswith('$('):
                                if "class: 'Directory'" in dirent.entry:
                                    # directory to create (dir = dirent.entryname) 
                                    pass
                                elif "class: 'File'" in dirent.entry:
                                    # file to create
                                    raise NotImplementedError
                                else:
                                    # most likely just specifying the file / directory needs to be staged. do nothing for now.
                                    # future: detect whether its a file or directory. create directory if directory. 
                                    res, success = parse_basic_expression(dirent.entry)
                                    if isinstance(res, InputSelector):
                                        print()
                                    pass 

                            elif dirent.entry.startswith('${'):
                                # js code to generate required file
                                msg = f'{dirent.entryname}: js code to dynamically create runtime file. please address'
                                self.error_msgs.append(msg)
                                out[dirent.entryname] = dirent.entry

                            else:
                                assert(isinstance(dirent.entryname, str))
                                if self.is_expression_tool:
                                    # the .js code script to run
                                    out[dirent.entryname] = dirent.entry
                                else:
                                    # either a static script, or a script to template at runtime
                                    # script to template at runtime
                                    if '$(inputs.' in dirent.entry:
                                        msg = f'likely untranslated cwl / js in script: {dirent.entryname}'
                                        self.error_msgs.append(msg)
                                        out[dirent.entryname] = dirent.entry
                                    # static script
                                    else:
                                        out[dirent.entryname] = dirent.entry

                        elif isinstance(item, dict):
                            raise NotImplementedError

                        else:
                            raise RuntimeError
        
        return out
    
    def get_directories_to_create(self, requirements: list[Any]) -> dict[str, Any]:
        return {}

    def get_env_vars(self, requirements: list[Any]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        
        for req in requirements:
            if isinstance(req, self.cwl_utils.EnvVarRequirement):
                for envdef in req.envDef:
                    # entry name
                    name_expr, success = parse_basic_expression(envdef.envName)
                    if not success:
                        msg = 'untranslated javascript expression in environment variable name'
                        self.error_msgs.append(msg)
                    # entry 
                    entry_expr, success = parse_basic_expression(envdef.envValue)
                    if not success:
                        msg = 'untranslated javascript expression in environment variable value'
                        self.error_msgs.append(msg)
                    out[name_expr] = entry_expr
        
        return out


class CLTParser(CLTEntityParser):
    
    failure_message: str = 'error parsing task. returned minimal task definition as fallback'

    def fallback(self) -> CommandTool:
        # inputs, arguments, outputs, requirements have fallbacks. 
        # error must have occured with other clt fields (io streams, base command, doc etc). 
        identifier = get_id_filename(self.entity.id)  # this does not have a fallback. really hope the error isnt here :/
        inputs = [self.ingest_command_tool_input(inp) for inp in self.entity.inputs]
        outputs = [self.ingest_command_tool_output(out) for out in self.entity.outputs]
        arguments = [self.ingest_command_tool_argument(arg) for arg in (self.entity.arguments or [])]
        
        jtool = CommandToolBuilder(
            tool=identifier,
            base_command=None,
            inputs=inputs,
            outputs=outputs,
            arguments=arguments,
            version="v0.1.0",
            doc=None,
            friendly_name=None,
            container="ubuntu:latest"
        )
        # this must be set for error messages to be associated with this entity
        self.uuid = jtool.uuid
    
        # requirements
        req_parser = CLTRequirementsParser(
            cwl_utils=self.cwl_utils, 
            entity=self.entity, 
            uuid=self.uuid, 
            is_expression_tool=self.is_expression_tool
        )
        requirements = req_parser.parse()
        
        jtool._directories_to_create = requirements['directories_to_create'] or None
        jtool._files_to_create = requirements['files_to_create'] or None
        jtool._env_vars = requirements['env_vars'] or None
        jtool._container = requirements['container']
        jtool._memory = requirements['memory']
        jtool._cpus = requirements['cpus']
        jtool._time = requirements['time']
        return jtool

    def do_parse(self) -> CommandTool:
        # convert yaml datatypes to python datatypes
        entity = convert_cwl_types_to_python(self.entity, self.cwl_utils)

        identifier = get_id_filename(self.entity.id)
        inputs = [self.ingest_command_tool_input(inp) for inp in self.entity.inputs]
        outputs = [self.ingest_command_tool_output(out) for out in self.entity.outputs]
        arguments = [self.ingest_command_tool_argument(arg) for arg in (self.entity.arguments or [])]

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
        # this must be set for error messages to be associated with this entity
        self.uuid = jtool.uuid

        # arguments & selector patterns for io stream names
        extra_args = self.ingest_io_streams(entity, jtool)
        
        # addressing cwl stdin: stdout: stderr: file naming
        jtool._arguments += extra_args

        # requirements
        req_parser = CLTRequirementsParser(
            cwl_utils=self.cwl_utils, 
            entity=self.entity, 
            uuid=self.uuid, 
            is_expression_tool=self.is_expression_tool
        )
        requirements = req_parser.parse()
        
        jtool._directories_to_create = requirements['directories_to_create'] or None
        jtool._files_to_create = requirements['files_to_create'] or None
        jtool._env_vars = requirements['env_vars'] or None
        jtool._container = requirements['container']
        jtool._memory = requirements['memory']
        jtool._cpus = requirements['cpus']
        jtool._time = requirements['time']
        return jtool
    
    def ingest_command_tool_argument(self, arg: Any) -> ToolArgument:
        parser = CLTArgumentParser(cwl_utils=self.cwl_utils, entity=arg)
        return parser.parse()

    def ingest_command_tool_input(self, inp: Any) -> ToolInput:
        parser = CLTInputParser(cwl_utils=self.cwl_utils, entity=inp)
        return parser.parse()
        
    def ingest_command_tool_output(self, out: Any, is_expression_tool: bool=False) -> ToolOutput:  
        parser = CLTOutputParser(cwl_utils=self.cwl_utils, entity=out)
        return parser.parse()
    
    def ingest_io_streams(self, entity: Any, jtool: CommandTool) -> list[ToolArgument]:
        out: list[ToolArgument] = []

        # n = last position for clt inputs / arguments
        n = self.get_last_input_position(jtool)
        
        # stderr: n + 1
        if entity.stderr:
            filename, success = parse_basic_expression(entity.stderr)
            if not success:
                filename = 'stderr.txt'
                self.error_msgs.append('untranslated javascript expression in stderr filename. used stderr.txt as fallback')
            arg = ToolArgument(prefix='2>', value=filename, position=n + 1)
            out.append(arg)
            self.apply_collection_to_stderr_types(filename, jtool)
        
        elif self.clt_has_stderr_outputs(entity):
            filename = 'stderr.txt'
            arg = ToolArgument(prefix='2>', value=filename, position=n + 1)
            out.append(arg)
            self.apply_collection_to_stderr_types(filename, jtool)
        
        # stdout: n + 2
        if entity.stdout:
            filename, success = parse_basic_expression(entity.stdout)
            if not success:
                filename = 'stdout.txt'
                self.error_msgs.append('untranslated javascript expression in stdout filename. used stdout.txt as fallback')
            arg = ToolArgument(prefix='>', value=filename, position=n + 2)
            out.append(arg)
            self.apply_collection_to_stdout_types(filename, jtool)
        
        # stdin: n + 3
        if entity.stdin:
            filename, success = parse_basic_expression(entity.stdin)
            if not success:
                self.error_msgs.append('untranslated javascript expression in stdin filename')
            arg = ToolArgument(prefix='<', value=filename, position=n + 3)
            out.append(arg)

        return out
    
    def clt_has_stderr_outputs(self, entity: Any) -> bool:
        clt = entity
        for out in clt.outputs:
            dtype, error_msgs = ingest_cwl_type(out.type, self.cwl_utils, secondary_files=out.secondaryFiles)
            self.error_msgs += error_msgs
            if isinstance(dtype, Stderr):
                return True
        return False
    
    def apply_collection_to_stderr_types(self, filename: str, jtool: CommandTool) -> None: 
        for out in jtool.outputs():
            if isinstance(out.output_type, Stderr):
                out.selector = WildcardSelector(filename)
    
    def apply_collection_to_stdout_types(self, filename: str, jtool: CommandTool) -> None: 
        for out in jtool.outputs():
            if isinstance(out.output_type, Stdout):
                out.selector = WildcardSelector(filename)

    def get_last_input_position(self, jtool: CommandTool) -> int:
        max_pos = 0
        inputs: list[ToolInput | ToolArgument] = copy.deepcopy(jtool.inputs())
        if jtool.arguments():
            inputs += copy.deepcopy(jtool.arguments())
        for inp in inputs:
            if inp.position and inp.position > max_pos:
                max_pos = inp.position
        return max_pos



@dataclass
class CLTArgumentParser(CLTEntityParser):

    failure_message: str = "error parsing a section of the task script. substituted this section for '[error]' as fallback"

    def fallback(self) -> ToolArgument:
        return ToolArgument('[error]')

    def do_parse(self) -> ToolArgument: 
        # I don't know when a clt argument would be a string
        if isinstance(self.entity, str):
            res, success = parse_basic_expression(self.entity)
            if not success:
                self.error_msgs.append('untranslated javascript expression')
            arg = ToolArgument(res)
        
        # normal case
        else:
            res, success = parse_basic_expression(self.entity.valueFrom)
            if not success:
                self.error_msgs.append('untranslated javascript expression')
            arg = ToolArgument(
                value=res,
                position=self.entity.position,
                prefix=self.entity.prefix,
                separate_value_from_prefix=self.entity.separate,
                shell_quote=self.entity.shellQuote,
            )
        
        # this must be set for error messages to be associated with this entity
        self.uuid = arg.uuid
        return arg


@dataclass
class CLTInputParser(CLTEntityParser):

    failure_message: str = 'error parsing task input. returned generic optional File input as fallback'

    def fallback(self) -> ToolInput:
        identifier = get_id_entity(self.entity.id) # hope the error isnt here lol
        return ToolInput(
            tag=identifier,
            input_type=File(optional=True),
            position=None,
            prefix=None,
            separate_value_from_prefix=None,
            separator=None,
            shell_quote=None,
            default=None,
            value=None
        )

    def do_parse(self) -> ToolInput:
        identifier = get_id_entity(self.entity.id)

        value: Optional[str] = None
        inpBinding = self.entity.inputBinding
        if inpBinding and inpBinding.valueFrom:
            if settings.ingest.cwl.INGEST_JAVASCRIPT_EXPRESSIONS: # why? 
                value, success = parse_basic_expression(inpBinding.valueFrom)
                if not success:
                    msg = f'untranslated javascript expression in inputBinding.valueFrom'
                    self.error_msgs.append(msg)
            else:
                pass
                # j.Logger.warn(
                #     f"Won't translate the expression for input {self.entity.id}: {inpBinding.valueFrom}"
                # )

        dtype, error_msgs = ingest_cwl_type(self.entity.type, self.cwl_utils, secondary_files=self.entity.secondaryFiles)
        self.error_msgs += error_msgs

        # edge case - InputParameter.default is File / Directory provided as cwl dict
        if isinstance(self.entity.default, dict):
            if 'class' in self.entity.default and 'location' in self.entity.default:
                self.entity.default = self.entity.default['location']
            else:
                print()

        inp = ToolInput(
            tag=identifier,
            input_type=dtype,
            position=inpBinding.position if inpBinding else None,
            prefix=inpBinding.prefix if inpBinding else None,
            separate_value_from_prefix=inpBinding.separate if inpBinding else None,
            separator=inpBinding.itemSeparator if inpBinding else None,
            shell_quote=inpBinding.shellQuote if inpBinding else None,
            default=self.entity.default,
            value=value
        )
        # this must be set for error messages to be associated with this entity
        self.uuid = inp.uuid
        return inp


@dataclass 
class CLTOutputParser(CLTEntityParser):

    failure_message: str = 'error parsing task output. returned generic optional File output as fallback'
    
    def fallback(self) -> ToolOutput:
        identifier = get_id_entity(self.entity.id)  # hope the error isnt here lol
        return ToolOutput(
            tag=identifier, 
            output_type=File(optional=True), 
            selector=WildcardSelector('*')
        )

    def do_parse(self) -> ToolOutput:
        # tag
        identifier = get_id_entity(self.entity.id)
        
        # datatype
        dtype, error_msgs = ingest_cwl_type(self.entity.type, self.cwl_utils, secondary_files=self.entity.secondaryFiles)
        self.error_msgs += error_msgs

        # selector
        selector = None
        if hasattr(self.entity, 'janis_collection_override'):
            selector = self.entity.janis_collection_override

        elif self.entity.outputBinding:
            if self.entity.outputBinding.glob:
                res, success = parse_basic_expression(self.entity.outputBinding.glob)
                selector = WildcardSelector(res)
                if not success:
                    self.error_msgs.append('untranslated javascript expression in output collection')
            
            elif self.entity.outputBinding.outputEval:
                res, success = parse_basic_expression(self.entity.outputBinding.outputEval)
                selector = res
                if not success:
                    self.error_msgs.append('untranslated javascript expression in output collection')
              
        out = ToolOutput(
            tag=identifier, 
            output_type=dtype, 
            selector=selector
        )
        # this must be set for error messages to be associated with this entity
        self.uuid = out.uuid
        return out

