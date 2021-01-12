import os
from abc import ABC, abstractmethod
from typing import Tuple, List, Dict
import functools

from path import Path

from janis_core.code.codetool import CodeTool
from janis_core.tool.commandtool import ToolInput
from janis_core.tool.tool import ToolType
from janis_core.translationdeps.exportpath import ExportPathKeywords
from janis_core.types.common_data_types import Int
from janis_core.utils import lowercase_dictkeys
from janis_core.utils.logger import Logger
from janis_core.operators.selectors import Selector


class TranslationError(Exception):
    def __init__(self, message, inner: Exception):
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
    """
    So you're thinking about adding a new tWranslation :)

    This class will hopefully give you a pretty good indication
    on what's required to add a new translation, however what I
    can't fully talk about are all the concepts.

    To add a new translation, you should understand at least these
    tool and janis concepts:
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

    I'd ask that you keep your function sizes small and direct,
    and then write unit tests to cover each component of the translation
    and then an integration test of the whole translation on the related workflows.

    You can find these in /janis/tests/test_translation_*.py)
    """

    __metaclass__ = TranslatorMeta

    def __init__(self, name):
        self.name = name

    def translate(
        self,
        tool,
        to_console=True,
        tool_to_console=False,
        with_resource_overrides=False,
        to_disk=False,
        write_inputs_file=True,
        export_path=ExportPathKeywords.default,
        should_validate=False,
        should_zip=True,
        merge_resources=False,
        hints=None,
        allow_null_if_not_optional=True,
        additional_inputs: Dict = None,
        max_cores=None,
        max_mem=None,
        max_duration=None,
        with_container=True,
        allow_empty_container=False,
        container_override=None,
    ):

        str_tool, tr_tools = None, []

        if tool.type() == ToolType.Workflow:
            tr_tool, tr_tools = self.translate_workflow(
                tool,
                with_container=with_container,
                with_resource_overrides=with_resource_overrides,
                allow_empty_container=allow_empty_container,
                container_override=lowercase_dictkeys(container_override),
            )
            str_tool = self.stringify_translated_workflow(tr_tool)
        elif isinstance(tool, CodeTool):
            tr_tool = self.translate_code_tool_internal(
                tool,
                allow_empty_container=allow_empty_container,
                container_override=lowercase_dictkeys(container_override),
            )
            str_tool = self.stringify_translated_tool(tr_tool)
        else:
            tr_tool = self.translate_tool_internal(
                tool,
                with_container=with_container,
                with_resource_overrides=with_resource_overrides,
                allow_empty_container=allow_empty_container,
                container_override=lowercase_dictkeys(container_override),
            )
            str_tool = self.stringify_translated_tool(tr_tool)

        tr_inp = self.build_inputs_file(
            tool,
            recursive=False,
            merge_resources=merge_resources,
            hints=hints,
            additional_inputs=additional_inputs,
            max_cores=max_cores,
            max_mem=max_mem,
            max_duration=max_duration,
        )
        tr_res = self.build_resources_input(tool, hints)

        str_inp = self.stringify_translated_inputs(tr_inp)
        str_tools = [
            (
                "tools/" + self.tool_filename(t),
                self.stringify_translated_workflow(tr_tools[t]),
            )
            for t in tr_tools
        ]
        str_resources = self.stringify_translated_inputs(tr_res)

        if to_console:
            print("=== WORKFLOW ===")
            print(str_tool)
            if tool_to_console:
                print("\n=== TOOLS ===")
                [print(f":: {t[0]} ::\n" + t[1]) for t in str_tools]
            print("\n=== INPUTS ===")
            print(str_inp)
            if not merge_resources and with_resource_overrides:
                print("\n=== RESOURCES ===")
                print(str_resources)

        d = ExportPathKeywords.resolve(
            export_path, workflow_spec=self.name, workflow_name=tool.versioned_id()
        )

        fn_workflow = self.workflow_filename(tool)
        fn_inputs = self.inputs_filename(tool)
        fn_resources = self.resources_filename(tool)

        if to_disk and write_inputs_file:
            if not os.path.isdir(d):
                os.makedirs(d)

            with open(os.path.join(d, fn_inputs), "w+") as f:
                Logger.log(f"Writing {fn_inputs} to disk")
                f.write(str_inp)
                Logger.log(f"Written {fn_inputs} to disk")
        else:
            Logger.log("Skipping writing input (yaml) job file")

        if to_disk:

            toolsdir = os.path.join(d, "tools")
            if not os.path.isdir(toolsdir):
                os.makedirs(toolsdir)

            Logger.info(f"Exporting tool files to '{d}'")

            with open(os.path.join(d, fn_workflow), "w+") as wf:
                Logger.log(f"Writing {fn_workflow} to disk")
                wf.write(str_tool)
                Logger.log(f"Wrote {fn_workflow}  to disk")

            for (fn_tool, disk_str_tool) in str_tools:
                with open(os.path.join(d, fn_tool), "w+") as toolfp:
                    Logger.log(f"Writing {fn_tool} to disk")
                    toolfp.write(disk_str_tool)
                    Logger.log(f"Written {fn_tool} to disk")

            if not merge_resources and with_resource_overrides:
                print("\n=== RESOURCES ===")
                with open(os.path.join(d, fn_resources), "w+") as wf:
                    Logger.log(f"Writing {fn_resources} to disk")
                    wf.write(str_inp)
                    Logger.log(f"Wrote {fn_resources}  to disk")
                print(str_resources)

            import subprocess

            if should_zip:
                Logger.debug("Zipping tools")
                with Path(d):
                    FNULL = open(os.devnull, "w")
                    zip_result = subprocess.run(
                        ["zip", "-r", "tools.zip", "tools/"], stdout=FNULL
                    )
                    if zip_result.returncode == 0:
                        Logger.debug("Zipped tools")
                    else:
                        Logger.critical(str(zip_result.stderr.decode()))

            if should_validate:
                with Path(d):

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

    def translate_tool(
        self,
        tool,
        to_console=True,
        to_disk=False,
        export_path=None,
        with_container=True,
        with_resource_overrides=False,
        max_cores=None,
        max_mem=None,
        allow_empty_container=False,
        container_override=None,
    ):

        tool_out = self.stringify_translated_tool(
            self.translate_tool_internal(
                tool,
                with_container=with_container,
                with_resource_overrides=with_resource_overrides,
                allow_empty_container=allow_empty_container,
                container_override=lowercase_dictkeys(container_override),
            )
        )

        if to_console:
            print(tool_out)

        if to_disk:
            d = ExportPathKeywords.resolve(
                export_path, workflow_spec=self.name, workflow_name=tool.id()
            )
            if not os.path.exists(d):
                os.makedirs(d)
            fn_tool = self.tool_filename(tool)
            with open(os.path.join(d, fn_tool), "w+") as wf:
                Logger.log(f"Writing {fn_tool} to disk")
                wf.write(tool_out)
                Logger.log(f"Wrote {fn_tool}  to disk")

        return tool_out

    def translate_code_tool(
        self,
        codetool,
        to_console=True,
        to_disk=False,
        export_path=None,
        with_docker=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override=None,
    ):

        tool_out = self.stringify_translated_tool(
            self.translate_code_tool_internal(
                codetool,
                with_docker=with_docker,
                allow_empty_container=allow_empty_container,
                container_override=lowercase_dictkeys(container_override),
            )
        )

        if to_console:
            print(tool_out)

        if to_disk:
            d = ExportPathKeywords.resolve(
                export_path, workflow_spec=self.name, workflow_name=codetool.id()
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
    def validate_inputs(cls, inputs, allow_null_if_optional):
        return True
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

    @staticmethod
    def get_type(t):
        if isinstance(t, list):
            q = set(TranslatorBase.get_type(tt) for tt in t)
            if len(q) == 0:
                return "empty array"
            val = q.pop() if len(q) == 1 else "Union[" + ", ".join(q) + "]"
            return f"Array<{val}>"

        return type(t).__name__

    @classmethod
    @abstractmethod
    def translate_workflow(
        cls,
        workflow,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ) -> Tuple[any, Dict[str, any]]:
        pass

    @classmethod
    @abstractmethod
    def translate_tool_internal(
        cls,
        tool,
        with_container=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        pass

    @classmethod
    @abstractmethod
    def translate_code_tool_internal(
        cls,
        tool,
        with_docker=True,
        allow_empty_container=False,
        container_override: dict = None,
    ):
        pass

    @classmethod
    @abstractmethod
    def unwrap_expression(cls, expression):
        pass

    @classmethod
    def build_inputs_file(
        cls,
        tool,
        recursive=False,
        merge_resources=False,
        hints=None,
        additional_inputs: Dict = None,
        max_cores=None,
        max_mem=None,
        max_duration=None,
    ) -> Dict[str, any]:

        ad = additional_inputs or {}
        values_provided_from_tool = {}
        if tool.type() == ToolType.Workflow:
            values_provided_from_tool = {
                i.id(): i.value or i.default
                for i in tool.input_nodes.values()
                if i.value or (i.default and not isinstance(i.default, Selector))
            }

        inp = {
            i.id(): ad.get(i.id(), values_provided_from_tool.get(i.id()))
            for i in tool.tool_inputs()
            if i.default is not None
            or not i.intype.optional
            or i.id() in ad
            or i.id() in values_provided_from_tool
        }

        if merge_resources:
            for k, v in cls.build_resources_input(
                tool,
                hints,
                max_cores=max_cores,
                max_mem=max_mem,
                max_duration=max_duration,
            ).items():
                inp[k] = ad.get(k, v)

        return inp

    @classmethod
    def build_resources_input(
        cls,
        tool,
        hints,
        max_cores=None,
        max_mem=None,
        max_duration=None,
        inputs=None,
        prefix="",
    ):

        inputs = inputs or {}

        if not tool.type() == ToolType.Workflow:
            cpus = inputs.get(f"{prefix}runtime_cpu", tool.cpus(hints) or 1)
            mem = inputs.get(f"{prefix}runtime_memory", tool.memory(hints))
            disk = inputs.get(f"{prefix}runtime_disks", 20)
            seconds = inputs.get(f"{prefix}runtime_seconds", 86400)

            if max_cores and cpus > max_cores:
                Logger.info(
                    f"Tool '{tool.id()}' exceeded ({cpus}) max number of cores ({max_cores}), "
                    "this was dropped to the new maximum"
                )
                cpus = max_cores
            if mem and max_mem and mem > max_mem:
                Logger.info(
                    f"Tool '{tool.id()}' exceeded ({mem} GB) max amount of memory ({max_mem} GB), "
                    "this was dropped to the new maximum"
                )
                mem = max_mem

            if seconds and max_duration and seconds > max_duration:
                Logger.info(
                    f"Tool '{tool.id()}' exceeded ({seconds} secs) max duration in seconds ({max_duration} secs), "
                    "this was dropped to the new maximum"
                )
                seconds = max_duration

            return {
                prefix + "runtime_memory": mem,
                prefix + "runtime_cpu": cpus,
                prefix + "runtime_disks": disk,
                prefix + "runtime_seconds": seconds,
            }

        new_inputs = {}
        for s in tool.step_nodes.values():
            new_inputs.update(
                cls.build_resources_input(
                    s.tool,
                    hints=hints,
                    max_cores=max_cores,
                    max_mem=max_mem,
                    max_duration=max_duration,
                    prefix=prefix + s.id() + "_",
                    inputs=inputs,
                )
            )

        return new_inputs

    @staticmethod
    def inp_can_be_skipped(inp, override_value=None):
        return (
            inp.default is None
            and override_value is None
            # and not inp.include_in_inputs_file_if_none
            and (inp.intype.optional and inp.default is None)
        )

    # Resource overrides
    @staticmethod
    def get_resource_override_inputs() -> List[ToolInput]:
        return [
            ToolInput("runtime_cpu", Int(optional=True)),  # number of CPUs
            ToolInput("runtime_memory", Int(optional=True)),  # GB of memory
            ToolInput("runtime_seconds", Int(optional=True)),  # seconds of running time
            ToolInput("runtime_disks", Int(optional=True)),  # GB of storage required
        ]

    # STRINGIFY

    @staticmethod
    @abstractmethod
    def stringify_translated_workflow(wf):
        pass

    @staticmethod
    @abstractmethod
    def stringify_translated_tool(tool):
        pass

    @staticmethod
    @abstractmethod
    def stringify_translated_inputs(inputs):
        pass

    # OUTPUTS

    @classmethod
    def filename(cls, tool):

        if tool.type() == ToolType.Workflow:
            return cls.workflow_filename(tool)
        return cls.tool_filename(tool)

    @staticmethod
    @abstractmethod
    def workflow_filename(workflow):
        pass

    @staticmethod
    @abstractmethod
    def inputs_filename(workflow):
        pass

    @staticmethod
    @abstractmethod
    def tool_filename(tool):
        pass

    @staticmethod
    def dependencies_filename(workflow):
        return "tools.zip"

    @staticmethod
    @abstractmethod
    def resources_filename(workflow):
        pass

    # VALIDATION

    @staticmethod
    @abstractmethod
    def validate_command_for(wfpath, inppath, tools_dir_path, tools_zip_path):
        pass

    @staticmethod
    def get_container_override_for_tool(tool, container_override):
        if not container_override:
            return None

        if tool.id().lower() in container_override:
            return container_override.get(tool.id().lower())
        elif tool.versioned_id().lower() in container_override:
            return container_override.get(tool.versioned_id().lower())
        elif "*" in container_override:
            return container_override["*"]
