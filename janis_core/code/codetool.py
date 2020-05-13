from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional, Union

from janis_core.types import Filename, String

from janis_core.tool.tool import Tool, TOutput, TInput, ToolType, ToolTypes


class CodeTool(Tool, ABC):

    # User should inherit from these blocks

    @abstractmethod
    def inputs(self) -> List[TInput]:
        pass

    @abstractmethod
    def outputs(self) -> List[TOutput]:
        pass

    def memory(self, hints: Dict[str, Any]) -> Optional[float]:
        """
        These values are used to generate a separate runtime.json / runtime.yaml input
        that can be passed to the execution engine to fill in for the specified hints.

        These are now (2019-04-10) to be kept out of the workflow, to leave the workflow
        truly portable.

        This memory must be in GB!
        :param hints: Dict[Key: value] of hints
        :return: Optional[int]
        """
        return None

    def cpus(self, hints: Dict[str, Any]) -> Optional[int]:
        """
        These values are used to generate a separate runtime.json / runtime.yaml input
        that can be passed to the execution engine to fill in for the specified hints.

        These are now (2019-04-10) to be kept out of the workflow, to leave the workflow
        truly portable.

        The CPU must be a whole number. If your tool contains threads
        :return:
        """
        return None

    # Janis developer should inherit these methods

    @abstractmethod
    def base_command(self):
        pass

    @abstractmethod
    def script_name(self):
        pass

    @abstractmethod
    def container(self):
        pass

    @abstractmethod
    def prepared_script(self):
        pass

    # Other internal methods

    @classmethod
    def type(cls) -> ToolType:
        return ToolTypes.CodeTool

    def containers(self) -> Dict[str, str]:
        return {self.versioned_id(): self.container()}

    def tool_inputs(self) -> List[TInput]:
        return self.inputs()

    def tool_outputs(self) -> List[TOutput]:
        return self.outputs()

    def generate_inputs_override(
        self,
        additional_inputs=None,
        with_resource_overrides=False,
        hints=None,
        include_defaults=True,
    ):
        d, ad = {}, additional_inputs or {}
        for i in self.inputs():
            if (
                not i.intype.optional
                or i.id() in ad
                or (include_defaults and i.default)
            ):
                d[i] = ad.get(i.id(), i.default)

        if with_resource_overrides:
            cpus = self.cpus(hints) or 1
            mem = self.memory(hints)
            d.update(
                {
                    "runtime_memory": mem,
                    "runtime_cpu": cpus,
                    "runtime_disks": "local-disk 60 SSD",
                }
            )

    def translate(
        self,
        translation: str,
        to_console=True,
        to_disk=False,
        export_path=None,
        with_docker=True,
        with_resource_overrides=False,
        allow_empty_container=False,
        container_override: Union[str, dict] = None,
    ):
        from janis_core import translations

        if isinstance(container_override, str):
            container_override = {self.id().lower(): container_override}

        return translations.translate_code_tool(
            self,
            translation=translation,
            to_console=to_console,
            to_disk=to_disk,
            with_docker=with_docker,
            export_path=export_path,
            allow_empty_container=allow_empty_container,
            container_override=container_override,
        )

    def wrapped_in_wf(self):
        from copy import copy
        from janis_core.workflow.workflow import WorkflowBuilder

        wf = WorkflowBuilder(self.id() + "Wf")
        inpmap = {}
        for i in self.inputs():

            if isinstance(i.intype, Filename):
                intp = String(optional=True)
            else:
                intp = copy(i.intype)
                if i.default:
                    intp.optional = True

            inpmap[i.id()] = wf.input(i.id(), intp)

        stp = wf.step(self.id().lower(), self(**inpmap))

        for o in self.outputs():
            wf.output(o.id(), source=stp[o.id()])

        return wf
