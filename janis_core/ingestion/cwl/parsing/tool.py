
import re 
import copy
from typing import Any, Optional, Tuple
from dataclasses import dataclass
from abc import ABC, abstractmethod
from functools import cached_property

from janis_core import ToolInput, ToolArgument, ToolOutput, WildcardSelector, CommandToolBuilder, CommandTool, Selector
from janis_core.types import File, Stdout, Stderr, Directory, DataType
from janis_core import settings
from janis_core.messages import log_error
from janis_core.messages import ErrorCategory

from ..types import ingest_cwl_type
from ..identifiers import get_id_entity
from ..identifiers import get_id_filename
from ..expressions import parse_expression


@dataclass
class CLTParser:
    cwl_utils: Any
    clt: Any
    entity: Any
    is_expression_tool: bool = False
    success: bool = False

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
        
        return j_entity

    def fallback(self) -> CommandTool:
        # inputs, arguments, outputs, requirements have fallbacks. 
        # error must have occured with other clt fields (io streams, base command, doc etc). 
        
        identifier = get_id_filename(self.entity.id) # this does not have a fallback. really hope no errors here :/
        jtool = CommandToolBuilder(
            tool=identifier,
            base_command=self.entity.baseCommand,
            inputs=[],
            outputs=[],
            arguments=[],
            version="v0.1.0",
            doc=None,
            friendly_name=None,
            container="ubuntu:latest"
        )
        # this must be set for error messages to be associated with this entity
        self.tool_uuid = jtool.uuid

        # log message
        msg = 'error parsing tool. returned minimal tool definition as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK)

        inputs = [self.ingest_command_tool_input(inp) for inp in self.entity.inputs]
        outputs = [self.ingest_command_tool_output(out) for out in self.entity.outputs]
        arguments = [self.ingest_command_tool_argument(arg) for arg in (self.entity.arguments or [])]
        
        inputs = [inp for inp in inputs if inp is not None]
        outputs = [out for out in outputs if out is not None]
        arguments = [arg for arg in arguments if arg is not None]

        jtool._inputs = inputs
        jtool._outputs = outputs
        jtool._arguments = arguments
    
        # requirements
        req_parser = CLTRequirementsParser(
            cwl_utils=self.cwl_utils,
            clt=self.clt,
            entity=self.entity, 
            tool_uuid=self.tool_uuid, 
            is_expression_tool=self.is_expression_tool
        )
        reqs = req_parser.parse()
        
        jtool._directories_to_create = reqs['directories_to_create'] or None # type: ignore
        jtool._files_to_create = reqs['files_to_create'] or None # type: ignore
        jtool._env_vars = reqs['env_vars'] or None # type: ignore
        jtool._container = reqs['container']
        jtool._memory = reqs['memory']
        jtool._cpus = reqs['cpus']
        jtool._time = reqs['time']
        return jtool

    def do_parse(self) -> CommandTool:
        identifier = get_id_filename(self.entity.id)
        jtool = CommandToolBuilder(
            tool=identifier,
            base_command=self.entity.baseCommand,
            inputs=[],
            outputs=[],
            arguments=[],
            version="v0.1.0",
            doc=self.entity.doc,
            friendly_name=self.entity.label,
            container="ubuntu:latest"
        )
        # this must be set for error messages to be associated with this entity
        self.tool_uuid = jtool.uuid

        inputs = [self.ingest_command_tool_input(inp) for inp in self.entity.inputs]
        outputs = [self.ingest_command_tool_output(out) for out in self.entity.outputs]
        arguments = [self.ingest_command_tool_argument(arg) for arg in (self.entity.arguments or [])]
        
        inputs = [inp for inp in inputs if inp is not None]
        outputs = [out for out in outputs if out is not None]
        arguments = [arg for arg in arguments if arg is not None]

        jtool._inputs = inputs
        jtool._outputs = outputs
        jtool._arguments = arguments

        if self.is_expression_tool:
            # ToolInput for script file, staging under correct name
            # ToolArgument to correctly supply inputs as argv to script
            # TODO HERE
            raise NotImplementedError

        # arguments & selector patterns for io stream names
        # addresses cwl stdin: stdout: stderr: file naming
        jtool._arguments += self.ingest_io_streams(self.entity, jtool)
        
        # requirements
        req_parser = CLTRequirementsParser(
            cwl_utils=self.cwl_utils, 
            clt=self.clt,
            entity=self.entity, 
            tool_uuid=self.tool_uuid, 
            is_expression_tool=self.is_expression_tool
        )
        requirements = req_parser.parse()
        
        jtool._directories_to_create = requirements['directories_to_create'] or None # type: ignore
        jtool._files_to_create = requirements['files_to_create'] or None # type: ignore
        jtool._env_vars = requirements['env_vars'] or None # type: ignore
        jtool._container = requirements['container']
        jtool._memory = requirements['memory']
        jtool._cpus = requirements['cpus']
        jtool._time = requirements['time']
        return jtool
    
    def ingest_command_tool_argument(self, arg: Any) -> Optional[ToolArgument]:
        parser = CLTArgumentParser(cwl_utils=self.cwl_utils, clt=self.clt, entity=arg, tool_uuid=self.tool_uuid)
        return parser.parse()

    def ingest_command_tool_input(self, inp: Any) -> Optional[ToolInput]:
        parser = CLTInputParser(cwl_utils=self.cwl_utils, clt=self.clt, entity=inp, tool_uuid=self.tool_uuid)
        return parser.parse()
        
    def ingest_command_tool_output(self, out: Any, is_expression_tool: bool=False) -> ToolOutput:  
        parser = CLTOutputParser(cwl_utils=self.cwl_utils, clt=self.clt, entity=out, tool_uuid=self.tool_uuid)
        return parser.parse()
    
    def ingest_io_streams(self, entity: Any, jtool: CommandTool) -> list[ToolArgument]:
        out: list[ToolArgument] = []

        # n = last position for clt inputs / arguments
        n = self.get_last_input_position(jtool)
        
        # stderr: n + 1
        if entity.stderr:
            filename, success = parse_expression(entity.stderr, self.tool_uuid)
            # if not success:
            #     msg = 'error parsing tool. returned minimal tool definition as fallback'
            #     log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK)
            #     filename = 'stderr.txt'
            #     self.error_msgs.append('untranslated javascript expression in stderr filename. used stderr.txt as fallback')
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
            filename, success = parse_expression(entity.stdout, self.tool_uuid)
            # if not success:
            #     filename = 'stdout.txt'
            #     self.error_msgs.append('untranslated javascript expression in stdout filename. used stdout.txt as fallback')
            arg = ToolArgument(prefix='>', value=filename, position=n + 2)
            out.append(arg)
            self.apply_collection_to_stdout_types(filename, jtool)
        
        # stdin: n + 3
        if entity.stdin:
            filename, success = parse_expression(entity.stdin, self.tool_uuid)
            # if not success:
            #     self.error_msgs.append('untranslated javascript expression in stdin filename')
            arg = ToolArgument(prefix='<', value=filename, position=n + 3)
            out.append(arg)

        return out
    
    def clt_has_stderr_outputs(self, entity: Any) -> bool:
        clt = entity
        for out in clt.outputs:
            dtype = ingest_cwl_type(out.type, self.cwl_utils, out, self.tool_uuid, secondaries=out.secondaryFiles)
            # self.error_msgs += error_msgs
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
class CLTEntityParser(ABC):
    cwl_utils: Any
    clt: Any
    entity: Any
    tool_uuid: str
    is_expression_tool: bool = False
    success: bool = False

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
        
        return j_entity

    @abstractmethod
    def do_parse(self) -> Any:
        ...
    
    @abstractmethod
    def fallback(self) -> Any:
        ...


@dataclass
class CLTRequirementsParser(CLTEntityParser):

    def fallback(self) -> dict[str, Any]:
        # log message
        msg = 'error parsing tool requirements. ignored requirements as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK)

        # fallback
        return {
            'container': None,
            'memory': None,
            'cpus': None,
            'time': None,
            'directories_to_create': [],
            'files_to_create': {},
            'env_vars': {}
        }

    def do_parse(self) -> dict[str, Any]:  
        self.handle_expression_libs()
        requirements = self.entity.requirements or []

        out: dict[str, Any] = {
            'container': None,
            'memory': None,
            'cpus': None,
            'time': None,
            'directories_to_create': [],
            'files_to_create': {},
            'env_vars': {}
        }
        
        out['container'] = self.get_container(requirements)
        out['memory'] = self.get_memory(requirements)
        out['cpus'] = self.get_cpus(requirements)
        out['time'] = self.get_time(requirements)
        out['files_to_create'], out['directories_to_create'] = self.get_files_directories_to_create(requirements)
        out['env_vars'] = self.get_env_vars(requirements)
        
        return out
        
    def handle_expression_libs(self) -> None:
        """
        goal is to mark expression libraries as messages at top of output file.
        achieve this by attempting to parse (will fail), causing a __TOKENX__ = "[JS EXPR]" to be created as a message. 
        """

        for req in self.entity.requirements:
            if req.__class__.__name__ == 'InlineJavascriptRequirement':
                if hasattr(req, 'expressionLib') and isinstance(req.expressionLib, list):
                    for exprlib in req.expressionLib:
                        assert isinstance(exprlib, str)
                        if not exprlib.startswith('$(') or not exprlib.startswith('$('):
                            exprlib = f'$({exprlib})'
                        res, success = parse_expression(exprlib, self.tool_uuid)

    def get_container(self, requirements: list[Any]) -> str:
        fallback = 'ubuntu:latest'
        if self.is_expression_tool:
            return 'node:latest'
        else:
            for req in requirements:
                if isinstance(req, self.cwl_utils.DockerRequirement):
                    return req.dockerPull
        return fallback
    
    def get_memory(self, requirements: list[Any]) -> Optional[int]:
        for req in requirements:
            if isinstance(req, self.cwl_utils.ResourceRequirement):
                # maybe convert mebibytes to megabytes?
                res, success = parse_expression(req.ramMin or req.ramMax, self.tool_uuid)
                return res
        return None
    
    def get_cpus(self, requirements: list[Any]) -> Optional[int]:
        for req in requirements:
            if isinstance(req, self.cwl_utils.ResourceRequirement):
                res, success = parse_expression(req.coresMin, self.tool_uuid)
                return res
        return None

    def get_time(self, requirements: list[Any]) -> Optional[int]:
        for req in requirements:
            if hasattr(req, 'timelimit') and isinstance(req, self.cwl_utils.ToolTimeLimit):
                res, success = parse_expression(req.timelimit, self.tool_uuid)
                return res
        return None
    
    def get_files_directories_to_create(self, requirements: list[Any]) -> Tuple[dict[str, Any], list[str | Selector]]:
        files_to_create: dict[str, Any] = {}
        directories_to_create: list[str | Selector] = []

        for req in requirements:
            if isinstance(req, self.cwl_utils.InitialWorkDirRequirement):
                # ensure array
                if isinstance (req.listing, str):
                    listing = [req.listing]
                else:
                    listing = req.listing
                
                for item in listing:
                    parser = InitialWorkDirRequirementParser(cwl_utils=self.cwl_utils, clt=self.entity, req=item, tool_uuid=self.tool_uuid, is_expression_tool=self.is_expression_tool)
                    parser.parse()
                    for fname, fcontents in parser.files_to_create:
                        files_to_create[fname] = fcontents
                    for dname in parser.directories_to_create:
                        directories_to_create.append(dname)
        
        return files_to_create, directories_to_create
    
    def get_env_vars(self, requirements: list[Any]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        
        for req in requirements:
            if isinstance(req, self.cwl_utils.EnvVarRequirement):
                for envdef in req.envDef:
                    name, success = parse_expression(envdef.envName, self.tool_uuid)
                    entry, success = parse_expression(envdef.envValue, self.tool_uuid)
                    out[name] = entry
        return out


class InitialWorkDirRequirementParser:
    """
    extracts files / directories to create from single InitialWorkDirRequirement entry
    requirement types: File | Directory | Dirent | string | Expression
    """
    def __init__(self, cwl_utils: Any, clt: Any, req: Any, tool_uuid: str, is_expression_tool: bool=False) -> None:
        self.cwl_utils = cwl_utils
        self.clt = clt
        self.req = req
        self.tool_uuid = tool_uuid
        self.is_expression_tool = is_expression_tool
        error_token_override = "error parsing InitialWorkDirRequirement"

        # resolving value, name
        self.r_value, self.r_value_ok = parse_expression(self.value, self.tool_uuid, error_token_override=error_token_override)
        if self.name is not None:
            self.r_name, self.r_name_ok = parse_expression(self.name, self.tool_uuid)
        else:
            self.r_name, self.r_name_ok = None, True
        
        # want to calculate these fields
        self.error_msgs: list[str] = []
        self.files_to_create: list[Tuple[str, str | Selector]] = []
        self.directories_to_create: list[str | Selector] = []

    def parse(self) -> None:
        # normal mode
        if settings.ingest.SAFE_MODE:
            try:
                j_entity = self.do_parse()
            except Exception:
                j_entity = self.fallback()
        # dev mode
        else:
            j_entity = self.do_parse()
        return j_entity
    
    def fallback(self) -> None:
        msg = 'error parsing InitialWorkDirRequirement. ignored as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK)
        self.files_to_create = [] 
        self.directories_to_create = [] 

    def do_parse(self) -> None:
        success = False
        parsers = [
            self.parse_as_textfile,
            self.parse_as_selector,
            self.parse_as_object,
            self.parse_as_directory_path,
        ]
        for parser in parsers:
            success = parser()
            if success:
                return None
        raise RuntimeError

    @cached_property
    def name(self) -> Optional[str]:
        name: Optional[str] = None
        if isinstance(self.req, self.cwl_utils.Dirent):
            name = self.req.entryname
        return name
    
    @cached_property
    def value(self) -> str:
        if isinstance(self.req, self.cwl_utils.Dirent):
            text: str = self.req.entry
        elif isinstance(self.req, str):
            text: str = self.req
        else:
            raise NotImplementedError
        return text
        
    def extract_input_type(self, sel: Selector) -> Optional[str]:
        pattern = r'inputs\.([a-zA-Z0-9_]+)'
        template = str(sel)
        match = re.match(pattern, template)
        if match:
            input_name = match.group(1)
            for inp in self.clt.inputs:
                name = get_id_entity(inp.id)
                if name == input_name:
                    dtype = ingest_cwl_type(inp.type, cwl_utils=self.cwl_utils, cwl_entity=inp, tool_uuid=self.tool_uuid, secondaries=inp.secondaryFiles)
                    if isinstance(dtype, Directory):
                        return 'dir'
                    else:
                        return 'file'
        return None

    def looks_like_expr(self, text: str) -> bool:
        if text.startswith('$(') or text.startswith('${'):
            return True
        pattern = r'^inputs\.([a-zA-Z0-9_]+)'
        match = re.match(pattern, text)
        if match:
            return True
        return False
    
    def parse_as_textfile(self) -> bool:
        fname_pattern = r'[a-zA-Z0-9_-]+\.(sh|py|r|R|rscript|Rscript|bash|txt|config|csv|tsv|json|yml|yaml|xml|html)'
        # must be dirent
        if not isinstance(self.req, self.cwl_utils.Dirent):
            return False
        
        # resolved entry and entryname must be strings
        if not isinstance(self.r_value, str) or not isinstance(self.r_name, str):
            return False
        
        # entryname and entry must be strings
        if self.looks_like_expr(self.req.entryname):
            return False
        if self.looks_like_expr(self.req.entry):
            return False
        
        # entryname must match a typical script file
        match = re.match(fname_pattern, self.req.entryname)
        if not match:
            return False
        
        # entry must have multiple lines 
        if '\n' not in self.req.entry:
            return False 
        
        # looks like script. do parse. 
        self.files_to_create.append((self.req.entryname, self.req.entry))
        return True
    
    def parse_as_selector(self) -> bool:
        # name and value must be resolvable
        if not self.r_name_ok or not self.r_value_ok:
            return False
        
        # resolved value must be a Selector
        if not isinstance(self.r_value, Selector):
            return False
        
        # resolved value must map to input
        input_type = self.extract_input_type(self.r_value)
        if input_type is None:
            return False 
        
        # file with entryname
        if input_type == 'file' and self.r_name is not None:
            self.files_to_create.append((str(self.r_name), self.r_value))
            return True 
        
        # file
        elif input_type == 'file':
            self.files_to_create.append((str(self.r_value), self.r_value))
            return True 
        
        # directory
        else:
            # TODO this is really weird 
            return False
            self.directories_to_create.append(str(self.r_value))
            return True 
    
    def parse_as_object(self) -> bool:
        # TODO?
        return False

    def parse_as_directory_path(self) -> bool:
        # name and value must be resolvable
        if not self.r_name_ok or not self.r_value_ok:
            return False
        
        # resolved value must be a string
        if not isinstance(self.r_value, str):
            return False
        
        # resolved name must be None
        if self.r_name is not None:
            return False
        
        # must look like directory
        # filename_pattern = r'^\/?([a-zA-Z0-9-_]+\/)*([a-zA-Z0-9-_.]+)$'
        dirname_pattern = r'^\/?([a-zA-Z0-9-_]+\/)*([a-zA-Z0-9-_]+\/?)$'
        if not re.match(dirname_pattern, self.r_value):
            return False 
        
        self.directories_to_create.append(self.r_value)
        return True 
    
    # def modify_js(self, script: str) -> str:
    #     text_to_remove = '"use strict";\nvar inputs=$(inputs);\nvar runtime=$(runtime);\n'
    #     text_to_substitute = '\nvar inputs = JSON.parse( process.argv[2] );\n\n'
    #     assert(text_to_remove in script)
    #     script = script.replace(text_to_remove, text_to_substitute)
    #     return script


@dataclass
class CLTArgumentParser(CLTEntityParser):


    def fallback(self) -> None:
        msg = 'error parsing CommandLineTool Argument. ignored as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK)
        return None

    def do_parse(self) -> Optional[ToolArgument]: 
        # I don't know when a clt argument would be a string
        if isinstance(self.entity, str):
            res, success = parse_expression(self.entity, self.tool_uuid)
            if res is None:
                return None
            arg = ToolArgument(res)
        
        # normal case
        else:
            res, success = parse_expression(self.entity.valueFrom, self.tool_uuid)
            if res is None:
                return None
            arg = ToolArgument(
                value=res,
                position=self.entity.position if self.entity.position else 0,
                prefix=self.entity.prefix,
                separate_value_from_prefix=self.entity.separate,
                shell_quote=self.entity.shellQuote,
            )
        
        # this must be set for error messages to be associated with this entity
        self.uuid = arg.uuid
        return arg


@dataclass
class CLTInputParser(CLTEntityParser):

    def fallback(self) -> ToolInput:
        identifier = get_id_entity(self.entity.id) # hope the error isnt here lol
        msg = 'error parsing CommandLineTool Input. returned generic optional File input as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK)

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
        tag = get_id_entity(self.entity.id)
        dtype = ingest_cwl_type(self.entity.type, self.cwl_utils, self.entity, self.tool_uuid, secondaries=self.entity.secondaryFiles)
        inp = ToolInput(
            tag=tag,
            input_type=dtype,
            position=self.parse_position(),
            prefix=self.parse_prefix(),
            separate_value_from_prefix=self.parse_separate(),
            separator=self.parse_separator(),
            shell_quote=self.parse_shell_quote(),
            default=self.parse_default(),
            value=self.parse_value()
        )
        return inp
    
    def parse_position(self) -> Optional[int]:
        if self.entity.inputBinding is not None:
            if self.entity.inputBinding.position is not None:
                return self.entity.inputBinding.position
        if self.entity.type is not None:
            if isinstance(self.entity.type, self.cwl_utils.InputArraySchema):
                if self.entity.type.inputBinding is not None:
                    if self.entity.type.inputBinding.position is not None:
                        return self.entity.type.inputBinding.position
        return None
    
    def parse_prefix(self) -> Optional[str]:
        # import cwl_utils.parser.cwl_v1_2 as cwlutils
        if self.entity.inputBinding is not None:
            if self.entity.inputBinding.prefix is not None:
                return self.entity.inputBinding.prefix
        if self.entity.type is not None:
            if isinstance(self.entity.type, self.cwl_utils.InputArraySchema):
                if self.entity.type.inputBinding is not None:
                    if self.entity.type.inputBinding.prefix is not None:
                        return self.entity.type.inputBinding.prefix
        return None
    
    def parse_separate(self) -> Optional[bool]:
        return self.entity.inputBinding.separate if self.entity.inputBinding else None
    
    def parse_separator(self) -> Optional[str]:
        return self.entity.inputBinding.itemSeparator if self.entity.inputBinding else None
    
    def parse_shell_quote(self) -> Optional[bool]:
        return self.entity.inputBinding.shellQuote if self.entity.inputBinding else None
    
    def parse_default(self) -> Any:
        # edge case - InputParameter.default is File / Directory provided as cwl dict
        if isinstance(self.entity.default, dict):
            if 'class' in self.entity.default and 'location' in self.entity.default:
                return self.entity.default['location']
        return self.entity.default
    
    def parse_value(self) -> Any:
        value = None
        if self.entity.inputBinding and self.entity.inputBinding.valueFrom:
            value, success = parse_expression(self.entity.inputBinding.valueFrom, self.tool_uuid)
        return value


@dataclass 
class CLTOutputParser(CLTEntityParser):

    def fallback(self) -> ToolOutput:
        # log message
        identifier = get_id_entity(self.entity.id) # hope the error isnt here lol
        msg = 'error parsing CommandLineTool output. returned generic File output as fallback'
        log_error(self.tool_uuid, msg, ErrorCategory.FALLBACK)

        # fallback
        return ToolOutput(
            tag=identifier, 
            output_type=File(optional=True), 
            selector=WildcardSelector('*')
        )

    def do_parse(self) -> ToolOutput:
        # tag
        identifier = get_id_entity(self.entity.id)
        # datatype
        dtype = ingest_cwl_type(self.entity.type, self.cwl_utils, self.entity, self.tool_uuid, secondaries=self.entity.secondaryFiles)

        if isinstance(dtype, Stdout):
            return self.do_parse_stdout(identifier, dtype)
        else:
            return self.do_parse_generic(identifier, dtype)

    def do_parse_stdout(self, identifier: str, dtype: DataType) -> ToolOutput:
        selector = None
        
        if self.clt.stdout is not None:
            expr, success = parse_expression(self.clt.stdout, self.tool_uuid)
            dtype = File()
            selector = expr
            
        out = ToolOutput(
            tag=identifier, 
            output_type=dtype, 
            selector=selector
        )
        return out
    
    def do_parse_generic(self, identifier: str, dtype: DataType) -> ToolOutput:
        # selector
        selector = None
        if hasattr(self.entity, 'janis_collection_override'):
            selector = self.entity.janis_collection_override

        elif self.entity.outputBinding:
            if self.entity.outputBinding.glob:
                res, success = parse_expression(self.entity.outputBinding.glob, self.tool_uuid)
                if isinstance(res, str) and not res.startswith('__TOKEN'):
                    selector = WildcardSelector(res)
                else:
                    selector = res
            
            elif self.entity.outputBinding.outputEval:
                res, success = parse_expression(self.entity.outputBinding.outputEval, self.tool_uuid)
                if isinstance(res, str) and not res.startswith('__TOKEN'):
                    selector = WildcardSelector(res)
                else:
                    selector = res
              
        out = ToolOutput(
            tag=identifier, 
            output_type=dtype, 
            selector=selector
        )
        return out

