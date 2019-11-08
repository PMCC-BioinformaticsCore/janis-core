from typing import List, Type, Optional
import pkg_resources, sys
from inspect import isfunction, ismodule, isabstract, isclass

from janis_core.tool.tool import Tool, ToolTypes
from janis_core.types.data_types import DataType
from janis_core.utils.logger import Logger, LogLevel
import janis_core.registry.entrypoints as EP
from janis_core.registry.registry import TaggedRegistry, Registry


class JanisShed:
    _byclassname = Registry()
    _toolshed = TaggedRegistry("latest")
    _typeshed = Registry()

    _has_been_hydrated = False

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

        Logger.log("Setting CONSOLE_LEVEL to None while traversing modules")
        cl = Logger.CONSOLE_LEVEL
        Logger.set_console_level(None)
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
        ep = []
        eps = pkg_resources.iter_entry_points(group=EP.DATATYPES)
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
        ep = []
        eps = pkg_resources.iter_entry_points(group=EP.TOOLS)
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
    def traverse_module(module, seen_modules: set, seen_classes: set):
        if module.__name__ in seen_modules:
            return

        seen_modules.add(module.__name__)

        q = {
            n: cls
            for n, cls in list(module.__dict__.items())
            if not n.startswith("__")
            and type(cls) != type
            and not (ismodule(cls) and cls.__name__ in seen_modules)
            and not (isclass(cls) and cls in seen_classes)
        }

        for k in q:
            cls = q[k]
            JanisShed.process_cls(cls, seen_modules, seen_classes)

    @staticmethod
    def process_cls(cls, seen_modules, seen_classes: set):
        try:
            if ismodule(cls):
                return JanisShed.traverse_module(cls, seen_modules, seen_classes)
            elif isfunction(cls) or isabstract(cls):
                return
            elif not isclass(cls):
                return

            seen_classes.add(cls)

            if issubclass(cls, DataType):
                return JanisShed.add_type(cls)
            elif not hasattr(cls, "type") or not callable(cls.type):
                return
            elif cls.type() == ToolTypes.Workflow:
                return JanisShed.add_tool(cls())
            elif cls.type() == ToolTypes.CommandTool:
                return JanisShed.add_tool(cls())

        except Exception as e:
            Logger.log(f"{str(e)} for type {type(cls)}")
