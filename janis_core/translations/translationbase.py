from abc import ABC, abstractmethod
import os
from path import Path
from typing import Tuple, List, Dict, Optional, Any

from janis_core.tool.tool import ToolInput
from janis_core.types.common_data_types import Int
from janis_core.translations.exportpath import ExportPathKeywords
from janis_core.utils.logger import Logger

# from janis_core.workflow.input import Input, InputNode


class TranslatorBase(ABC):
    """
    So you're thinking about adding a new translation :)

    This class will hopefully give you a pretty good indication
    on what's required to add a new translation, however what I
    can't fully talk about are all the concepts.

    To add a new translation, you should understand at least these
    workflow and janis concepts:
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

    def __init__(self, name):
        self.name = name

    def translate(
        self,
        workflow,
        to_console=True,
        tool_to_console=False,
        with_docker=True,
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
    ):

        # self.validate_inputs(workflow._inputs, allow_null_if_not_optional)

        tr_wf, tr_tools = self.translate_workflow(
            workflow,
            with_docker=with_docker,
            with_resource_overrides=with_resource_overrides,
        )
        tr_inp = self.build_inputs_file(
            workflow,
            recursive=False,
            merge_resources=merge_resources,
            hints=hints,
            additional_inputs=additional_inputs,
            max_cores=max_cores,
            max_mem=max_mem,
        )
        tr_res = self.build_resources_input(workflow, hints)

        str_wf = self.stringify_translated_workflow(tr_wf)
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
            print(str_wf)
            if tool_to_console:
                print("\n=== TOOLS ===")
                [print(f":: {t[0]} ::\n" + t[1]) for t in str_tools]
            print("\n=== INPUTS ===")
            print(str_inp)
            if not merge_resources and with_resource_overrides:
                print("\n=== RESOURCES ===")
                print(str_resources)

        d = ExportPathKeywords.resolve(
            export_path, workflow_spec=self.name, workflow_name=workflow.id()
        )

        fn_workflow = self.workflow_filename(workflow)
        fn_inputs = self.inputs_filename(workflow)
        fn_resources = self.resources_filename(workflow)

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

            Logger.info(f"Exporting workflow files to '{d}'")

            with open(os.path.join(d, fn_workflow), "w+") as wf:
                Logger.log(f"Writing {fn_workflow} to disk")
                wf.write(str_wf)
                Logger.log(f"Wrote {fn_workflow}  to disk")

            for (fn_tool, str_tool) in str_tools:
                with open(os.path.join(d, fn_tool), "w+") as toolfp:
                    Logger.log(f"Writing {fn_tool} to disk")
                    toolfp.write(str_tool)
                    Logger.log(f"Written {fn_tool} to disk")

            if not merge_resources and with_resource_overrides:
                print("\n=== RESOURCES ===")
                with open(os.path.join(d, fn_resources), "w+") as wf:
                    Logger.log(f"Writing {fn_resources} to disk")
                    wf.write(str_wf)
                    Logger.log(f"Wrote {fn_resources}  to disk")
                print(str_resources)

            import subprocess

            if should_zip:
                Logger.info("Zipping tools")
                with Path(d):
                    zip_result = subprocess.run(["zip", "-r", "tools.zip", "tools/"])
                    if zip_result.returncode == 0:
                        Logger.info("Zipped tools")
                    else:
                        Logger.critical(zip_result.stderr)

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
                            "Exported workflow was validated by: " + " ".join(enved_vcs)
                        )
                    else:
                        Logger.critical(cwltool_result.stderr)

        return str_wf, str_inp, str_tools

    def translate_tool(
        self,
        tool,
        to_console=True,
        to_disk=False,
        export_path=None,
        with_docker=True,
        with_resource_overrides=False,
        max_cores=None,
        max_mem=None,
    ):

        tool_out = self.stringify_translated_tool(
            self.translate_tool_internal(
                tool,
                with_docker=with_docker,
                with_resource_overrides=with_resource_overrides,
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
        cls, workflow, with_docker=True, with_resource_overrides=False
    ) -> Tuple[any, Dict[str, any]]:
        pass

    @classmethod
    @abstractmethod
    def translate_tool_internal(
        cls, tool, with_docker=True, with_resource_overrides=False
    ):
        pass

    @classmethod
    @abstractmethod
    def build_inputs_file(
        cls,
        workflow,
        recursive=False,
        merge_resources=False,
        hints=None,
        additional_inputs: Dict = None,
        max_cores=None,
        max_mem=None,
    ) -> Dict[str, any]:
        pass

    @classmethod
    def build_resources_input(
        cls, workflow, hints, max_cores=None, max_mem=None, inputs=None, prefix=""
    ):
        from janis_core.workflow.workflow import Workflow, Tool, CommandTool

        # returns a list of key, value pairs
        steps: Dict[str, Optional[Any]] = {}
        inputs = inputs or {}

        # prefix = ""
        # if not prefix:
        #     prefix = ""
        # else:
        #     prefix += ""

        for s in workflow.step_nodes.values():
            tool: Tool = s.tool

            if isinstance(tool, CommandTool):
                tool_pre = prefix + s.id() + "_"
                cpus = inputs.get(f"{s.id()}_runtime_cpu", tool.cpus(hints) or 1)
                mem = inputs.get(f"{s.id()}_runtime_memory", tool.memory(hints))
                disk = inputs.get(f"{s.id()}_runtime_disks", "local-disk 60 SSD")

                if max_cores and cpus > max_cores:
                    Logger.info(
                        f"Tool '{tool.tool()}' exceeded ({cpus}) max number of cores ({max_cores}), "
                        "this was dropped to the new maximum"
                    )
                    cpus = max_cores
                if max_mem and mem > max_mem:
                    Logger.info(
                        f"Tool '{tool.tool()}' exceeded ({mem} GB) max amount of memory({max_mem} GB), "
                        "this was dropped to the new maximum"
                    )
                    mem = max_mem
                steps.update(
                    {
                        tool_pre + "runtime_memory": mem,
                        tool_pre + "runtime_cpu": cpus,
                        tool_pre + "runtime_disks": disk,
                    }
                )
            elif isinstance(tool, Workflow):
                tool_pre = prefix + s.id() + "_"
                steps.update(
                    cls.build_resources_input(
                        tool,
                        hints,
                        max_cores=max_cores,
                        max_mem=max_mem,
                        prefix=tool_pre,
                        inputs={
                            k[len(s.id()) + 1 :]: v
                            for k, v in inputs.items()
                            if k.startswith(s.id())
                        },
                    )
                )

        return steps

    @staticmethod
    def inp_can_be_skipped(inp, override_value=None):
        return (
            inp.default is None
            and override_value is None
            # and not inp.include_in_inputs_file_if_none
            and (inp.datatype.optional and inp.default is None)
        )

    # Resource overrides
    @staticmethod
    def get_resource_override_inputs() -> List[ToolInput]:
        return [
            ToolInput("runtime_cpu", Int(optional=True), default=1),
            ToolInput("runtime_memory", Int(optional=True), default=4),
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
