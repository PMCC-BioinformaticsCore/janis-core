from typing import List, Type, Optional
from inspect import isfunction, ismodule, isabstract, isclass

from janis_core.tool.commandtool import Tool, ToolTypes, CommandTool, CommandToolBuilder
from janis_core.code.pythontool import CodeTool, PythonTool
from janis_core.workflow.workflow import Workflow, WorkflowBuilder
from janis_core.types.data_types import DataType
from janis_core.utils.logger import Logger, LogLevel
import janis_core.toolbox.entrypoints as EP
from janis_core.toolbox.register import TaggedRegistry, Registry


class JanisShed:

    MAX_RECURSION_DEPTH = 4

    _byclassname = Registry()
    _toolshed = TaggedRegistry("latest")
    _typeshed = Registry()

    _has_been_hydrated = False

    should_trace = False
    recognised_types = {ToolTypes.Workflow, ToolTypes.CommandTool, ToolTypes.CodeTool}

    # getters

    @staticmethod
    def get_by_class_name(name: str):
        JanisShed.hydrate()
        return JanisShed._byclassname.get(name)

    @staticmethod
    def get_tool(tool: str, version: str = None):
        JanisShed.hydrate()
        if version:
            version = version.lower()
        return JanisShed._toolshed.get(tool.lower(), version)

    @staticmethod
    def get_datatype(datatype: str):
        JanisShed.hydrate()
        return JanisShed._typeshed.get(datatype.lower())

    @staticmethod
    def get_all_tools() -> List[List[Tool]]:
        JanisShed.hydrate()
        return JanisShed._toolshed.objects()

    @staticmethod
    def get_all_datatypes() -> List[Type[DataType]]:
        JanisShed.hydrate()
        return JanisShed._typeshed.objects()

    # setters

    @staticmethod
    def add_tool(tool: Tool) -> bool:
        v: Optional[str] = tool.version()
        if not v:
            t = f"The tool {tool.id()} did not have a version and will not be registered"
            Logger.critical(t)
            return False
        Logger.log("Adding tool: " + tool.id())

        JanisShed._byclassname.register(tool.__class__.__name__, tool)
        return JanisShed._toolshed.register(tool.id().lower(), v.lower(), tool)

    @staticmethod
    def add_type(datatype: Type[DataType]) -> bool:
        JanisShed._byclassname.register(datatype.__name__, datatype)
        return JanisShed._typeshed.register(datatype.name().lower(), datatype)

    @staticmethod
    def hydrate(force=False, modules: list = None):
        # go get everything
        if JanisShed._has_been_hydrated and not force:
            return

        if not modules:
            modules = []
            modules.extend(JanisShed._get_datatype_entrypoints())
            modules.extend(JanisShed._get_tool_entrypoints())

        level = None
        cl = Logger.CONSOLE_LEVEL
        if JanisShed.should_trace:
            level = cl if cl >= LogLevel.DEBUG else LogLevel.DEBUG
        Logger.log(
            f"Setting CONSOLE_LEVEL to {LogLevel.get_str(level) or 'None'} while traversing modules"
        )
        Logger.set_console_level(level)
        seen_modules = set()
        seen_classes = set()
        for m in modules:
            JanisShed.traverse_module(
                m, seen_modules=seen_modules, seen_classes=seen_classes
            )
        Logger.set_console_level(cl)
        Logger.log(
            f"Restoring CONSOLE_LEVEL to {LogLevel.get_str(cl)} now that Janis shed has been hydrated"
        )

        JanisShed._has_been_hydrated = True

    @staticmethod
    def _get_datatype_entrypoints():
        import importlib_metadata

        ep = []
        eps = importlib_metadata.entry_points().get(EP.DATATYPES, [])
        for entrypoint in eps:
            try:
                m = entrypoint.load()
                ep.append(m)
            except ImportError as e:
                t = f"Couldn't import janis data_type extension '{entrypoint.name}': {e}"
                Logger.critical(t)
                continue
        return ep

    @staticmethod
    def _get_tool_entrypoints():
        import importlib_metadata

        ep = []
        eps = importlib_metadata.entry_points().get(EP.TOOLS, [])
        for entrypoint in eps:
            try:
                m = entrypoint.load()
                ep.append(m)
            except ImportError as e:
                t = f"Couldn't import janis data_type extension '{entrypoint.name}': {e}"
                Logger.critical(t)
                continue
        return ep

    @staticmethod
    def traverse_module(module, seen_modules: set, seen_classes: set, current_layer=1):
        if module.__name__ in seen_modules:
            return
        Logger.log("Traversing module " + str(module.__name__))
        seen_modules.add(module.__name__)

        q = {
            n: cls
            for n, cls in list(module.__dict__.items())
            if not n.startswith("__")
            and type(cls) != type
            and not (ismodule(cls) and cls.__name__ in seen_modules)
            and (not isinstance(cls, list) and cls not in seen_classes)
        }

        for k in q:
            cls = q[k]
            JanisShed.process_cls(cls, seen_modules, seen_classes, current_layer)

    @staticmethod
    def process_cls(cls, seen_modules, seen_classes: set, current_layer: int):
        try:
            if ismodule(cls):
                if current_layer <= JanisShed.MAX_RECURSION_DEPTH:
                    return JanisShed.traverse_module(
                        cls, seen_modules, seen_classes, current_layer=current_layer + 1
                    )
                return Logger.log(
                    f"Skip traversing module '{str(cls)}' as reached maximum depth ({JanisShed.MAX_RECURSION_DEPTH})"
                )
            elif isfunction(cls):
                return

            seen_classes.add(cls)
            if isclass(cls) and issubclass(cls, DataType):
                return JanisShed.add_type(cls)
            elif not hasattr(cls, "type") or not callable(cls.type):
                return

            if (
                cls == Tool
                or cls == Workflow
                or cls == CommandTool
                or cls == CodeTool
                or cls == PythonTool
                or cls == WorkflowBuilder
                or cls == CommandToolBuilder
            ):
                return

            tp = cls.type()
            if isinstance(tp, str) and tp in JanisShed.recognised_types:
                if isabstract(cls):
                    if issubclass(cls, Tool):
                        abstractmethods = list(cls.__abstractmethods__)
                        return Logger.warn(
                            f"The tool '{cls.__name__}' had abstract methods: "
                            + ", ".join(abstractmethods)
                        )
                    return
                ic = cls() if isclass(cls) else cls
                return JanisShed.add_tool(ic)

        except Exception as e:
            Logger.log(f"{str(e)} for type {str(cls)}")
