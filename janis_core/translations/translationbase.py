
import os
from abc import ABC, abstractmethod
from typing import Tuple, Any, Optional
import functools

from janis_core import (
    CommandToolBuilder, 
    CodeTool,
    PythonTool, 
    WorkflowBuilder, 
    Tool, 
    InputSelector, 
    Selector,
    SupportedTranslation, 
    ToolInput,
    ToolType,
    File, 
    Directory, 
    Int,
)

from janis_core.translation_deps.exportpath import ExportPathKeywords
from janis_core.utils.logger import Logger
from janis_core import settings

from janis_core.messages import inject_messages
from .common.todisk import write_tool_to_console
from .common.todisk import write_tool_to_disk
from .common.todisk import write_workflow_to_console
from .common.todisk import write_workflow_to_disk

class TranslationError(Exception):
    def __init__(self, message: str, inner: Exception):
        super().__init__(message, inner)


kwargstoignore = {"container_override"}



def try_catch_translate(type):
    def try_catch_translate_inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except TranslationError:
                raise
            except Exception as e:

                components = ", ".join(
                    [
                        *[repr(a) for a in args],
                        *[
                            f"{k}={v}"
                            for k, v in kwargs.items()
                            if k not in kwargstoignore
                        ],
                    ]
                )
                message = f"Couldn't translate {type or ''} with ({components})"
                er = TranslationError(message, inner=e)
                Logger.log_ex(er)
                raise er

        return wrapper

    return try_catch_translate_inner


class TranslatorMeta(type(ABC)):
    def __repr__(cls):
        return cls.__name__

    def __str__(cls):
        return cls.__name__




class TranslatorBase(ABC):
    """
    Base class for translators. 
    Anything that inherits from this class will be picked up by the translation system.
    Anything which is abstract will need to be implemented by the inheriting class.

    To add a new translation, you will need to understand the internal Janis model 
    and the target language. tool and janis concepts:
        - inputs
        - outputs
        - steps + tools
        - secondary files
        - command line binding of tools
            - position, quoted
        - selectors
            - input selectors
            - wildcard glob selectors
            - cpu and memory selectors
        - propagated resource overrides

    In addition there are some modules which provide utilities.
    For example, the messaging module allows you to fetch any messages that were logged 
    during ingestion, so you can include these at the top of output files. 

    you can find tests for other translate units in /janis_core/tests/test_translation_*.py
    """

    __metaclass__ = TranslatorMeta
    basedir: str = ''

    # this must be implemented in the inheriting class if you want to change the 
    # output directory structure. defaults seen below. 
    OUTDIR_STRUCTURE: dict[str, str | None] = {
        'main': None,
        'subworkflows': 'subworkflows',
        'tools': 'tools',
        'inputs': None,
        'helpers': 'files',
        'resources': None,
    }

    def __init__(self, name: str):
        self.name: str = name                           # translation destination
        self.main: Optional[Tuple[Any, Any]] = None     # (internal, translated)
        
        # use the add_tool() and add_subworkflow() methods to add to these:
        self.subworkflows: list[Tuple] = []             # list of (internal, translated)
        self.tools: list[Tuple] = []                    # list of (internal, translated)
        
        self.inputs_file: Optional[str] = None               # contents of input config file
        self.resources_file: Optional[str] = None            # contents of resources config file
        self.helper_files: list[Tuple[str, str]] = []        # list of (filename: contents)

    
    ### IMPLEMENTED HELPER METHODS ###

    def add_tool(self, internal: CommandToolBuilder | CodeTool, translated: Any) -> None:
        """ensures unique tools - tool can be translated more than once"""
        for item in self.tools:
            if item[0].id() == internal.id():
                return None
        self.tools.append((internal, translated))
    
    def get_tool(self, query: CommandToolBuilder | CodeTool) -> Optional[CommandToolBuilder | CodeTool]:
        """used to check if tool already translated to save runtime"""
        for item in self.tools:
            if item[0].id() == query.id():
                return item[1]
        return None
    
    def add_subworkflow(self, internal: WorkflowBuilder, translated: Any) -> None:
        """ensures unique subworkflows - subworkflow can be translated more than once"""
        for item in self.subworkflows:
            if item[0].id() == internal.id():
                return None
        self.subworkflows.append((internal, translated))
    
    def get_subworkflow(self, query: WorkflowBuilder) -> Optional[WorkflowBuilder]:
        """used to check if subworkflow already translated to save runtime"""
        for item in self.subworkflows:
            if item[0].id() == query.id():
                return item[1]
        return None
    
    
    ### IMPLEMENTED PUBLIC METHODS ###

    def translate_workflow(self, wf: WorkflowBuilder) -> Any:
        """
        For workflow translation mode. 
        This method is only called once during execution. 
        Performs translation of a janis WorkflowBuilder.
        Calls abstract methods implemented in child class to do the translation. 
        Once translation is complete, writes to stdout and/or disk.
        Returns stringified translation purely for testing purposes
        """
        # do translation to output model
        self.translate_workflow_internal(wf)
        self.build_inputs_file(wf)     
        self.build_resources_file(wf) 
        self.build_helper_files(wf) 
        assert self.inputs_file is not None
        assert self.main is not None

        # stringify main workflow 
        internal, translated = self.main
        fn_main = self.workflow_filename(wf, is_main=True)
        str_main = self.stringify_translated_workflow(internal, translated)
        str_main = inject_messages(internal, str_main)
        tup_main = (fn_main, str_main)

        # stringify tools (commandtools, pythontools)
        tup_tools = []
        for internal, translated in self.tools:
            filename = self.tool_filename(internal)
            str_tool = self.stringify_translated_tool(internal, translated)
            str_tool = inject_messages(internal, str_tool)
            tup_tools.append((filename, str_tool))
        
        # stringify subworkflows
        tup_subworkflows = []
        for internal, translated in self.subworkflows:
            filename = self.workflow_filename(internal)
            str_subwf = self.stringify_translated_workflow(internal, translated)
            str_subwf = inject_messages(internal, str_subwf)
            tup_subworkflows.append((filename, str_subwf))

        # stringify input config file
        fn_inputs = self.inputs_filename(wf)
        str_inputs = self.inputs_file                # already a string
        tup_inputs = (fn_inputs, str_inputs)

        # stringify resource config file
        if self.resources_file is not None:
            fn_resources = self.resources_filename(wf)
            str_resources = self.resources_file      # already a string
            tup_resources = (fn_resources, str_resources)
        else:
            tup_resources = None

        # write to stdout / disk
        if settings.translate.TO_CONSOLE:
            write_workflow_to_console(tup_main, tup_tools, tup_inputs, tup_resources)
        if settings.translate.TO_DISK:
            write_workflow_to_disk(
                tup_main, 
                tup_subworkflows, 
                tup_tools, 
                tup_inputs, 
                self.helper_files, 
                tup_resources, 
                self.OUTDIR_STRUCTURE
            )

        return tup_main[1], tup_inputs[1], tup_subworkflows, tup_tools

    def translate_tool(self, internal: CommandToolBuilder) -> str:
        """
        For tool translation mode. 
        This method is only called once during execution. 
        Performs translation of a janis CommandToolBuilder.
        Calls abstract methods implemented in child class to do the translation. 
        Once translation is complete, writes to stdout and/or disk.
        Returns string purely for testing purposes
        """
        # do translation to output model
        self.translate_tool_internal(internal)
        assert len(self.tools) == 1
        translated = self.tools[0][1]
        
        # stringify output model
        str_tool = self.stringify_translated_tool(internal, translated)
        str_tool = inject_messages(internal, str_tool)
        
        # write to stdout / disk
        if settings.translate.TO_CONSOLE:
            write_tool_to_console(str_tool)
        if settings.translate.TO_DISK:
            filename = self.tool_filename(internal)
            self.build_helper_files(internal)
            write_tool_to_disk(str_tool, filename, self.helper_files)

        return str_tool 

    def translate_code_tool(self, internal: CodeTool):
        """
        For tool translation mode. 
        This method is only called once during execution. 
        Performs translation of a Janis CodeTool.
        Calls abstract methods implemented in child class to do the translation. 
        Once translation is complete, writes to stdout and/or disk.
        Returns string purely for testing purposes
        """
        self.translate_code_tool_internal(internal)
        assert len(self.tools) == 1
        translated = self.tools[0][1]
        
        str_tool = self.stringify_translated_tool(internal, translated)
        str_tool = inject_messages(internal, str_tool)

        if settings.translate.TO_CONSOLE:
            print(str_tool)

        if settings.translate.TO_DISK:
            d = ExportPathKeywords.resolve(
                settings.translate.EXPORT_PATH, workflow_spec=self.name, workflow_name=internal.id()
            )
            if not os.path.exists(d):
                os.makedirs(d)
            fn_tool = self.tool_filename(internal)
            with open(os.path.join(d, fn_tool), "w+") as wf:
                Logger.log(f"Writing {fn_tool} to disk")
                wf.write(str_tool)
                Logger.log(f"Wrote {fn_tool}  to disk")

        return str_tool

    # Resource overrides?? what for???
    @staticmethod
    def get_resource_override_inputs() -> list[ToolInput]:
        """not really sure what this does / why is needed"""
        return [
            ToolInput("runtime_cpu", Int(optional=True)),  # number of CPUs
            ToolInput("runtime_memory", Int(optional=True)),  # GB of memory
            ToolInput("runtime_seconds", Int(optional=True)),  # seconds of running time
            ToolInput("runtime_disk", Int(optional=True)),  # GB of storage required
        ]

    @staticmethod
    def get_container_override_for_tool(tool: CommandToolBuilder | CodeTool):
        """this is an implementation detail - have been delegated to child classes. """
        container_override = settings.translate.CONTAINER_OVERRIDE
        if not container_override:
            return None
        if tool.id().lower() in container_override:
            return container_override.get(tool.id().lower())
        elif tool.versioned_id().lower() in container_override:
            return container_override.get(tool.versioned_id().lower())
        elif "*" in container_override:
            return container_override["*"]

    def build_helper_files(self, tool: CommandToolBuilder | CodeTool | WorkflowBuilder) -> None:
        """
        Generate a dictionary of helper files to run Nextflow.
        Key of the dictionary is the filename, the value is the file content

        :param tool:
        :type tool:
        :return:
        :rtype:
        """
        code_files = self._gen_python_code_files(tool)
        template_files = self._gen_template_files(tool)
        helpers = template_files | code_files
        helpers = [(k, v) for k, v in helpers.items()]
        self.helper_files = helpers
    
    def _gen_python_code_files(self, tool: CommandToolBuilder | CodeTool | WorkflowBuilder) -> dict[str, str]:
        # Python files for Python code tools
        files: dict[str, str] = {}

        if isinstance(tool, PythonTool):
            # helpers["__init__.py"] = ""
            #helpers[f"{tool.versioned_id()}.py"] = self.gen_python_script(tool)
            filename = f'{tool.id()}.py'

            st = SupportedTranslation.from_str(self.name)
            contents = tool.prepared_script(st)
            assert contents is not None
            files[filename] = contents

        elif isinstance(tool, WorkflowBuilder):
            for step in tool.step_nodes.values():
                step_code_files = self._gen_python_code_files(step.tool)
                files = files | step_code_files # merge dicts
        
        return files

    def _gen_template_files(self, tool: CommandToolBuilder | CodeTool | WorkflowBuilder) -> dict[str, str]:
        # files from tool.files_to_create
        files: dict[str, str] = {}

        if isinstance(tool, CommandToolBuilder):
            if tool.files_to_create():
                for name, contents in tool.files_to_create().items():
                    if not isinstance(name, str):
                        # If name is a File or Directory, the entryname field overrides the value of basename of the File or Directory object 
                        raise NotImplementedError()
                    
                    if name == 'script.sh' and tool.is_shell_script:
                        # don't create externally passed script for shell tools
                        continue 
                    
                    if isinstance(contents, str):
                        files[name] = contents
                    
                    elif isinstance(contents, InputSelector):
                        tinput_name = contents.input_to_select
                        tinput = tool.inputs_map()[tinput_name]
                        if isinstance(tinput.intype, File | Directory):
                            raise NotImplementedError()
                        else:
                            raise NotImplementedError()
                    else:
                        raise NotImplementedError
        
        elif isinstance(tool, WorkflowBuilder):
            for step in tool.step_nodes.values():
                step_template_files = self._gen_template_files(step.tool)
                files = files | step_template_files # merge dicts

        return files


    ### JANIS ->  OUTPUT MODEL MAPPING ###
    # everything below this point is abstract. 
    # child classes must implement these.

    @abstractmethod
    def translate_workflow_internal(self, wf: WorkflowBuilder) -> Any:
        """
        Translate a janis model WorkflowBuilder to the target output model.
        This mapping is specific to the target language so must be implemented in child class. 
        Can be called 2+ times during execution in the case of subworkflows. 
        """
        pass

    @abstractmethod
    def translate_tool_internal(self, tool: CommandToolBuilder) -> Any:
        """
        Translate a janis model CommandToolBuilder to the target output model.
        This mapping is specific to the target language so must be implemented in child class. 
        Can be called 1+ times during execution. 
        """
        pass

    @abstractmethod
    def translate_code_tool_internal(self, tool: CodeTool) -> Any:
        """
        Translate a janis model CodeTool to the target output model.
        This mapping is specific to the target language so must be implemented in child class. 
        Can be called 1+ times during execution. 
        """
        pass
    
    @abstractmethod
    def build_inputs_file(self, entity: WorkflowBuilder | CommandToolBuilder | CodeTool) -> None:
        """
        Generates the inputs configuration file (as str) for the main translated workflow.
        (or tool in the case of single tool translation).
        For example:
            - CWL:      {id}-inp.yml
            - Nextflow: nextflow.config
            - WDL:      {id}-inp.json
        
        Assigns to self.inputs.
        """
        ...
    
    @abstractmethod
    def build_resources_file(self, entity: WorkflowBuilder | CommandToolBuilder | CodeTool) -> None:
        """
        Generates a separate resources configuration file (as str) for the main translated workflow.
        (or tool in the case of single tool translation).
        
        This is in addition to the inputs configuration file. 
        Some languages separate the inputs & resource files (eg CWL/WDL) whereas others do not (eg Nextflow).
        If you would like to merge these into 1 file, check that the settings.translate.MERGE_RESOURCES flag 
        is true somewhere in your implementation. 

        For example:
            - CWL:      {id}-resources.yml
            - WDL:      {id}-resources.json
        
        Assigns to self.resources.
        """
        ...

    @classmethod
    def _build_resources_dict(
        cls, 
        tool: Tool, 
        inputs: Optional[dict[str, Any]]=None,
        prefix: str=""
    ) -> dict[str, Any]:
        """
        A (mostly) shared implementation to get resources as structured dict for CWL and WDL. 
        Used in build_resources_file() and build_inputs_file().
        """
        inputs = inputs or {}
        hints = settings.translate.HINTS
        max_cores = settings.translate.MAX_CORES
        max_mem = settings.translate.MAX_MEM
        max_duration = settings.translate.MAX_DURATION

        if not tool.type() == ToolType.Workflow:
            cpus = inputs.get(f"{prefix}runtime_cpu", tool.cpus(hints))
            if cpus is None:
                cpus = 1
            mem = inputs.get(f"{prefix}runtime_memory", tool.memory(hints))
            disk = inputs.get(f"{prefix}runtime_disk", 20)
            seconds = inputs.get(f"{prefix}runtime_seconds", 86400)

            if max_cores is not None and cpus > max_cores:
                Logger.info(
                    f"Tool '{tool.id()}' exceeded ({cpus}) max number of cores ({max_cores}), "
                    "this was dropped to the new maximum"
                )
                cpus = max_cores
            if mem is not None and max_mem and mem > max_mem:
                Logger.info(
                    f"Tool '{tool.id()}' exceeded ({mem} GB) max amount of memory ({max_mem} GB), "
                    "this was dropped to the new maximum"
                )
                mem = max_mem

            if seconds is not None and max_duration and seconds > max_duration:
                Logger.info(
                    f"Tool '{tool.id()}' exceeded ({seconds} secs) max duration in seconds ({max_duration} secs), "
                    "this was dropped to the new maximum"
                )
                seconds = max_duration

            return {
                prefix + "runtime_memory": mem
                if not isinstance(mem, Selector)
                else None,
                prefix + "runtime_cpu": cpus
                if not isinstance(cpus, Selector)
                else None,
                prefix + "runtime_disk": disk
                if not isinstance(disk, Selector)
                else None,
                prefix + "runtime_seconds": seconds
                if not isinstance(seconds, Selector)
                else None,
            }

        new_inputs = {}
        for s in tool.step_nodes.values():
            step_inputs = cls._build_resources_dict(
                s.tool,
                prefix=prefix + s.id() + "_",
                inputs=inputs,
            )
            new_inputs.update(step_inputs)

        return new_inputs

    @classmethod
    @abstractmethod
    def unwrap_expression(cls, expression: Any) -> Any:
        """
        This is poor design. This is an implementation detail, should not be part of the interface.
        Why is it a classmethod too??? These class methods suck. 
        
        Anyway, this is supposed to be a recursive function to map Janis Operators to the target language.
        Eg [Janis -> CWL]: InputSelector('infile') -> ${inputs.infile}
        Eg [Janis -> NXF]: IndexOperator(InputSelector('infile'), 0) -> ${infile[0]}
        """
        ...
    
    
    ### STRINGIFY ###

    @abstractmethod
    def stringify_translated_workflow(self, internal: WorkflowBuilder, translated: Any) -> str:
        """
        Output Model -> String.
        Generates file contents for the main workflow (or subworkflow) so can be written to disk. 
        Messages to user should be injected here. 
        """
        ...

    @abstractmethod
    def stringify_translated_tool(self, internal: CommandToolBuilder | CodeTool, translated: Any) -> str:
        """
        Output Model -> String.
        Generates file contents for tool so can be written to disk. 
        Messages to user should be injected here. 
        """
        ...

    
    ### FILENAMES ###

    @staticmethod
    @abstractmethod
    def workflow_filename(workflow: WorkflowBuilder, is_main: Optional[bool]=False) -> str:
        """returns the filename for a workflow so file can be written to disk"""
        ...

    @staticmethod
    @abstractmethod
    def tool_filename(tool: CommandToolBuilder | CodeTool) -> str:
        """returns the filename for a tool so file can be written to disk"""
        ...

    @staticmethod
    @abstractmethod
    def inputs_filename(workflow: WorkflowBuilder) -> str:
        """returns the filename for the inputs configuration so file can be written to disk"""
        ...

    @staticmethod
    @abstractmethod
    def resources_filename(workflow: WorkflowBuilder) -> str:
        """returns the filename for the resources configuration so file can be written to disk"""
        ...

    
    ### VALIDATION ###
    # could be quite useful for user messages - eg "this tool doesn't seem valid"

    @staticmethod
    @abstractmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        pass





# DEPRECATED

# # do translation to output model
        # tr_workflow, tr_tools = self.translate_workflow_internal(wf)
        # tr_helpers = self.translate_helper_files(wf)
        # tr_inp = self.build_inputs_dict(wf)
        # tr_res = self.build_resources_input(wf)

        # # stringify main workflow file
        # fn_main = self.workflow_filename(wf)
        # str_main = self.stringify_translated_workflow(wf, tr_workflow)
        # tup_main = (fn_main, str_main)

        # # stringify subfiles (commandtools, pythontools & subworkflows)
        # # TODO separate tools and subworkflows
        # tup_tools = []
        # for t in tr_tools.keys():
        #     filename = self.tool_filename(t)
        #     str_tool = self.stringify_translated_workflow(tr_tools[t])
        #     tup_tools.append((filename, str_tool))

        # # stringify auxiliary files
        # tup_helpers = []
        # for filename in tr_helpers.keys():
        #     str_helper = tr_helpers[filename]
        #     tup_helpers.append((filename, str_helper))
        
        # # stringify input config file
        # fn_inputs = self.inputs_filename(wf)
        # str_inputs = self.stringify_translated_inputs(tr_inp)
        # tup_inputs = (fn_inputs, str_inputs)

        # # stringify resource config file
        # fn_resources = self.resources_filename(wf)
        # str_resources = self.stringify_translated_inputs(tr_res)
        # tup_resources = (fn_resources, str_resources)

        # # write to stdout / disk
        # if settings.translate.TO_CONSOLE:
        #     write_workflow_to_console(tup_main, tup_tools, tup_inputs, tup_resources)
        # if settings.translate.TO_DISK:
        #     tup_subworkflows = [] # TODO separate tools and subworkflows
        #     write_workflow_to_disk(
        #         tup_main, 
        #         tup_subworkflows, 
        #         tup_tools, 
        #         tup_inputs, 
        #         tup_helpers, 
        #         tup_resources, 
        #         self.OUTDIR_STRUCTURE
        #     )

        # return tup_main[1], tup_inputs[1], tup_tools





    # @classmethod
    # def validate_inputs(cls, inputs, allow_null_if_optional):
    #     return True
        # invalid = [
        #     i
        #     for i in inputs
        #     if not i.input.validate_value(
        #         allow_null_if_not_optional=allow_null_if_optional
        #     )
        # ]
        # if len(invalid) == 0:
        #     return True
        # raise TypeError(
        #     "Couldn't validate inputs: "
        #     + ", ".join(
        #         f"{i.id()} (expected: {i.input.data_type.id()}, "
        #         f"got: '{TranslatorBase.get_type(i.input.value)}')"
        #         for i in invalid
        #     )
        # )

    # @staticmethod
    # def get_type(t):
    #     if isinstance(t, list):
    #         q = set(TranslatorBase.get_type(tt) for tt in t)
    #         if len(q) == 0:
    #             return "empty array"
    #         val = q.pop() if len(q) == 1 else "Union[" + ", ".join(q) + "]"
    #         return f"Array<{val}>"

    #     return type(t).__name__


    # OUTPUTS

    # @classmethod
    # def filename(cls, tool):
    #     if tool.type() == ToolType.Workflow:
    #         return cls.workflow_filename(tool)
    #     return cls.tool_filename(tool)


    # @staticmethod
    # def dependencies_filename(workflow: WorkflowBuilder) -> str:
    #     return "tools.zip"



    # @abstractmethod
    # def build_inputs_file(self, entity: WorkflowBuilder | CommandToolBuilder | CodeTool) -> None:
    #     """
    #     Generates the input config file (as str) for this workflow or tool.*
    #     Assigns to self.inputs.

    #     *in the case of single tool translation
    #     """
    #     ...
        # ad = settings.translate.ADDITIONAL_INPUTS or {}
        # values_provided_from_tool = {}
        
        # # workflow translation
        # if entity.type() == ToolType.Workflow:
        #     values_provided_from_tool = {
        #         i.id(): i.value or i.default
        #         for i in entity.input_nodes.values()
        #         if i.value is not None
        #         or (i.default is not None and not isinstance(i.default, Selector))
        #     }

        # # tool translation
        # inp = {
        #     i.id(): ad.get(i.id(), values_provided_from_tool.get(i.id()))
        #     for i in entity.tool_inputs()
        #     if (i.default is not None and not isinstance(i.default, Selector))
        #     or not i.intype.optional
        #     or i.id() in ad
        #     or i.id() in values_provided_from_tool
        # }

        # if settings.translate.MERGE_RESOURCES:
        #     res_input = self.build_resources_file(entity)
        #     for k, v in res_input.items():
        #         inp[k] = ad.get(k, v)

        # return inp