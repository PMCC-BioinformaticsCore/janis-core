import os
from abc import ABC, abstractmethod
from typing import Tuple, List, Dict, Any, Optional
import functools
import shutil

from path import Path
from janis_core import CommandTool, CodeTool, WorkflowBase, Tool
from janis_core.code.codetool import CodeTool
from janis_core.tool.commandtool import ToolInput
from janis_core.tool.tool import ToolType
from janis_core.translation_deps.exportpath import ExportPathKeywords
from janis_core.types.common_data_types import Int
from janis_core.utils.logger import Logger
from janis_core.operators.selectors import Selector
from janis_core import settings

class TranslationError(Exception):
    def __init__(self, message: str, inner: Exception):
        super().__init__(message, inner)
        # self.inner = inner


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

    __metaclass__ = TranslatorMeta
    basedir: str = ''

    DIR_TOOLS: str = "tools"
    DIR_FILES: str = "files"
    SUBDIRS_TO_CREATE: list[str] = []  # this is if you want to write tools / workflows to subfolders

    def __init__(self, name: str):
        self.name: str = name

    def translate_workflow(self, wf: WorkflowBase):
        self.basedir = ExportPathKeywords.resolve(
            settings.translate.EXPORT_PATH, workflow_spec=self.name, workflow_name=wf.versioned_id()
        )
        str_tool, tr_tools, tr_helpers = None, [], {}

        # GENERATE MAIN FILE
        tr_workflow, tr_tools = self.translate_workflow_internal(wf)
        str_tool = self.stringify_translated_workflow(tr_workflow)

        # GENERATE SUBFILES - COMMANDTOOLS, PYTHONTOOLS & SUBWORKFLOWS
        # [filepath, filecontents] for subfiles (tools, subworkflows etc)
        str_tools = [
            (
                os.path.join(self.DIR_TOOLS, self.tool_filename(t)),
                self.stringify_translated_workflow(tr_tools[t]),
            )
            for t in tr_tools
        ]

        # GENERATE AUXILIARY FILES
        # {filepath: filecontents} for auxiliary files (PythonTool code.py files etc)
        tr_helpers = self.translate_helper_files(wf)
        str_helpers = [
            (os.path.join(self.DIR_FILES, filename), tr_helpers[filename])
            for filename in tr_helpers.keys()
        ]

        # GENERATE INPUT CONFIG
        # {name: value} for inputs config file (nextflow.config, inputs.yaml etc)
        tr_inp = self.build_inputs_dict(wf)
        
        # inputs config file contents
        str_inp = self.stringify_translated_inputs(tr_inp)

        # {name: value} for resource inputs config
        tr_res = self.build_resources_input(wf)
        # resource config file contents
        str_resources = self.stringify_translated_inputs(tr_res)

        # WRITING TO CONSOLE
        if settings.translate.TO_CONSOLE:
            print(str_tool)
            if settings.translate.TOOL_TO_CONSOLE:
                print("\n=== TOOLS ===")
                [print(f":: {t[0]} ::\n" + t[1]) for t in str_tools]
            print("\n=== INPUTS ===")
            print(str_inp)
            if not settings.translate.MERGE_RESOURCES and settings.translate.WITH_RESOURCE_OVERRIDES:
                print("\n=== RESOURCES ===")
                print(str_resources)

        # WRITING TO DISK        
        if settings.translate.TO_DISK:
            # setting filepaths
            basedir = self.basedir
            if os.path.isdir(basedir):
                shutil.rmtree(basedir)
            fn_workflow = self.workflow_filename(wf)
            fn_inputs = self.inputs_filename(wf)
            fn_resources = self.resources_filename(wf)

            # generating subfolders
            subfolders: list[str] = []
            subfolders.append(self.DIR_TOOLS)
            subfolders += self.SUBDIRS_TO_CREATE
            for subfolder in subfolders:
                path = os.path.join(basedir, subfolder)
                if not os.path.isdir(path):
                    os.makedirs(path)

            # writing inputs config file
            if settings.translate.WRITE_INPUTS_FILE:
                if not os.path.isdir(basedir):
                    os.makedirs(basedir)

                with open(os.path.join(basedir, fn_inputs), "w+") as f:
                    Logger.log(f"Writing {fn_inputs} to disk")
                    f.write(str_inp)
                    Logger.log(f"Written {fn_inputs} to disk")
            else:
                Logger.log("Skipping writing input (yaml) job file")

            # writing resources config file
            if not settings.translate.MERGE_RESOURCES and settings.translate.WITH_RESOURCE_OVERRIDES:
                print("\n=== RESOURCES ===")
                with open(os.path.join(basedir, fn_resources), "w+") as wf:
                    Logger.log(f"Writing {fn_resources} to disk")
                    wf.write(str_inp)
                    Logger.log(f"Wrote {fn_resources}  to disk")
                print(str_resources)

            # writing workflow / tool files
            Logger.info(f"Exporting tool files to '{basedir}'")

            # writing main workflow
            with open(os.path.join(basedir, fn_workflow), "w+") as wf:
                Logger.log(f"Writing {fn_workflow} to disk")
                wf.write(str_tool)
                Logger.log(f"Wrote {fn_workflow}  to disk")

            # writing tools, subworkflows
            for (fn_tool, disk_str_tool) in str_tools:
                path = os.path.join(basedir, fn_tool)
                with open(path, "w+") as toolfp:
                    Logger.log(f"Writing {fn_tool} to disk")
                    toolfp.write(disk_str_tool)
                    Logger.log(f"Written {fn_tool} to disk")
            
            # copying source files 
            if settings.general.SOURCE_FILES is not None:
                # create source folder in basedir
                source_dir = os.path.join(basedir, 'source')
                if not os.path.isdir(source_dir):
                    os.mkdir(source_dir)
                
                # copy files
                for src, dest in settings.general.SOURCE_FILES:
                    dest = os.path.join(basedir, 'source', dest)
                    if not os.path.isdir(os.path.dirname(dest)):
                        os.mkdir(os.path.dirname(dest))
                    shutil.copy2(src, dest)

            # writing helper files 
            for (fn_helper, disk_str_helper) in str_helpers:
                with open(os.path.join(basedir, fn_helper), "w+") as helperfp:
                    Logger.log(f"Writing {fn_helper} to disk")
                    helperfp.write(disk_str_helper)
                    Logger.log(f"Written {fn_helper} to disk")

            # zipping tools file
            import subprocess

            if settings.translate.SHOULD_ZIP:
                Logger.debug("Zipping tools")
                with Path(basedir):
                    FNULL = open(os.devnull, "w")
                    zip_result = subprocess.run(
                        ["zip", "-r", "tools.zip", "tools/"], stdout=FNULL
                    )
                    if zip_result.returncode == 0:
                        Logger.debug("Zipped tools")
                    else:
                        Logger.critical(str(zip_result.stderr.decode()))

            if settings.translate.SHOULD_VALIDATE:
                with Path(basedir):

                    Logger.info(f"Validating outputted {self.name}")

                    enved_vcs = [
                        (os.getenv(x[1:]) if x.startswith("$") else x)
                        for x in self.validate_command_for(
                            fn_workflow, fn_inputs, "tools/", "tools.zip"
                        )
                    ]

                    cwltool_result = subprocess.run(enved_vcs)
                    if cwltool_result.returncode == 0:
                        Logger.info(
                            "Exported tool was validated by: " + " ".join(enved_vcs)
                        )
                    else:
                        Logger.critical(str(cwltool_result.stderr))

        return str_tool, str_inp, str_tools

    def translate_tool(self, tool: CommandTool):
        tr_tool = self.translate_tool_internal(tool)
        tool_out = self.stringify_translated_tool(tr_tool)

        if settings.translate.TO_CONSOLE:
            print(tool_out)

        if settings.translate.TO_DISK:
            # set output folder
            basedir = ExportPathKeywords.resolve(
                settings.translate.EXPORT_PATH, workflow_spec=self.name, workflow_name=tool.id()
            )
            
            # create output folder
            if not os.path.exists(basedir):
                os.makedirs(basedir)

            # write tool file
            fn_tool = self.tool_filename(tool)
            with open(os.path.join(basedir, fn_tool), "w+") as wf:
                Logger.log(f"Writing {fn_tool} to disk")
                wf.write(tool_out)
                Logger.log(f"Wrote {fn_tool}  to disk")

            # copy source files to output folder
            if settings.general.SOURCE_FILES is not None:
                # create source folder in basedir
                source_dir = os.path.join(basedir, 'source')
                if not os.path.isdir(source_dir):
                    os.mkdir(source_dir)
                
                # copy files
                for src, dest in settings.general.SOURCE_FILES:
                    dest = os.path.join(basedir, 'source', dest)
                    if not os.path.isdir(os.path.dirname(dest)):
                        os.mkdir(os.path.dirname(dest))
                    shutil.copy2(src, dest)

        return tool_out

    def translate_code_tool(self, codetool: CodeTool):
        tr_tool = self.translate_code_tool_internal(codetool)
        tool_out = self.stringify_translated_tool(tr_tool)

        if settings.translate.TO_CONSOLE:
            print(tool_out)

        if settings.translate.TO_DISK:
            d = ExportPathKeywords.resolve(
                settings.translate.EXPORT_PATH, workflow_spec=self.name, workflow_name=codetool.id()
            )
            if not os.path.exists(d):
                os.makedirs(d)
            fn_tool = self.tool_filename(codetool)
            with open(os.path.join(d, fn_tool), "w+") as wf:
                Logger.log(f"Writing {fn_tool} to disk")
                wf.write(tool_out)
                Logger.log(f"Wrote {fn_tool}  to disk")

        return tool_out

    @classmethod
    def build_inputs_dict(cls, tool: CommandTool | CodeTool | WorkflowBase) -> dict[str, Any]:
        ad = settings.translate.ADDITIONAL_INPUTS or {}
        values_provided_from_tool = {}
        if tool.type() == ToolType.Workflow:
            values_provided_from_tool = {
                i.id(): i.value or i.default
                for i in tool.input_nodes.values()
                if i.value is not None
                or (i.default is not None and not isinstance(i.default, Selector))
            }

        inp = {
            i.id(): ad.get(i.id(), values_provided_from_tool.get(i.id()))
            for i in tool.tool_inputs()
            if (i.default is not None and not isinstance(i.default, Selector))
            or not i.intype.optional
            or i.id() in ad
            or i.id() in values_provided_from_tool
        }

        if settings.translate.MERGE_RESOURCES:
            res_input = cls.build_resources_input(tool)
            for k, v in res_input.items():
                inp[k] = ad.get(k, v)

        return inp

    @classmethod
    def build_resources_input(
        cls, 
        tool: Tool, 
        inputs: Optional[dict[str, Any]]=None,
        prefix: str=""
    ) -> dict[str, Any]:

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
            step_inputs = cls.build_resources_input(
                s.tool,
                prefix=prefix + s.id() + "_",
                inputs=inputs,
            )
            new_inputs.update(step_inputs)

        return new_inputs

    # Resource overrides
    @staticmethod
    def get_resource_override_inputs() -> List[ToolInput]:
        return [
            ToolInput("runtime_cpu", Int(optional=True)),  # number of CPUs
            ToolInput("runtime_memory", Int(optional=True)),  # GB of memory
            ToolInput("runtime_seconds", Int(optional=True)),  # seconds of running time
            ToolInput("runtime_disk", Int(optional=True)),  # GB of storage required
        ]

    @staticmethod
    def get_container_override_for_tool(tool: CommandTool | CodeTool):
        container_override = settings.translate.CONTAINER_OVERRIDE
        if not container_override:
            return None
        if tool.id().lower() in container_override:
            return container_override.get(tool.id().lower())
        elif tool.versioned_id().lower() in container_override:
            return container_override.get(tool.versioned_id().lower())
        elif "*" in container_override:
            return container_override["*"]

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

    # CHILD CLASS SPECIFIC

    @classmethod
    @abstractmethod
    def translate_workflow_internal(cls, wf: WorkflowBase) -> Tuple[Any, dict[str, Any]]:
        pass

    @classmethod
    @abstractmethod
    def translate_tool_internal(cls, tool: CommandTool) -> Any:
        pass

    @classmethod
    @abstractmethod
    def translate_code_tool_internal(cls, tool: CodeTool) -> Any:
        pass

    # overridden in NextflowTranslator
    @classmethod
    def translate_helper_files(cls, tool: CommandTool | CodeTool | WorkflowBase) -> Dict[str, str]:
        return {}

    @classmethod
    @abstractmethod
    def unwrap_expression(cls, expression: Any) -> Any:
        pass
    
    # STRINGIFY

    @staticmethod
    @abstractmethod
    def stringify_translated_workflow(wf: WorkflowBase) -> str:
        pass

    @staticmethod
    @abstractmethod
    def stringify_translated_tool(tool: CommandTool | CodeTool) -> str:
        pass

    @staticmethod
    @abstractmethod
    def stringify_translated_inputs(inputs: dict[str, Any]) -> str:
        pass

    # OUTPUTS

    # @classmethod
    # def filename(cls, tool):
    #     if tool.type() == ToolType.Workflow:
    #         return cls.workflow_filename(tool)
    #     return cls.tool_filename(tool)

    @staticmethod
    @abstractmethod
    def workflow_filename(workflow: WorkflowBase) -> str:
        pass

    @staticmethod
    @abstractmethod
    def inputs_filename(workflow: WorkflowBase) -> str:
        pass

    @staticmethod
    @abstractmethod
    def tool_filename(tool: CommandTool | CodeTool) -> str:
        pass

    # @staticmethod
    # def dependencies_filename(workflow: WorkflowBase) -> str:
    #     return "tools.zip"

    @staticmethod
    @abstractmethod
    def resources_filename(workflow: WorkflowBase) -> str:
        pass

    # VALIDATION

    @staticmethod
    @abstractmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        pass


