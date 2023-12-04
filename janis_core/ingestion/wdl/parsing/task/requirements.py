
import WDL 
import janis_core as j
import regex as re 
from typing import Any , Optional

from janis_core import Selector
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory

from ..expressions import parse_expr
from ..EntityParser import TaskParser


DEFAULT_CONTAINER = "ubuntu:latest"

class TaskContainerParser(TaskParser):

    def do_parse(self) -> str:
        container = self.task.runtime.get("container", self.task.runtime.get("docker"))
        if isinstance(container, WDL.Expr.Get):
            # relevant input
            # print(str(container.expr))
            inp = [i.expr for i in self.task.inputs if i.name == str(container.expr)]
            if len(inp) > 0:
                container = inp[0]
            else:
                msg = f"Expression for determining containers was '{container}' but couldn't find input called {str(container.expr)}"
                log_message(self.cmdtool.uuid, msg, category=ErrorCategory.SCRIPTING)
                raise RuntimeError
        if isinstance(container, WDL.Expr.String):
            container = container.literal
        if isinstance(container, WDL.Value.String):
            container = container.value
        if container is None:
            container = DEFAULT_CONTAINER
        if not isinstance(container, str):
            # TODO improve this using parse expr
            msg = f"Expression for determining containers ({container}) are not supported in Janis, using default container"
            log_message(self.cmdtool.uuid, msg, category=ErrorCategory.SCRIPTING)
            raise RuntimeError
        return container
    
    def fallback(self) -> str:
        msg = 'Error parsing container requirement'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return DEFAULT_CONTAINER


class TaskCpusParser(TaskParser):

    def do_parse(self) -> Optional[int]:
        value = self.task.runtime.get("cpu")
        if value is None:
            return None
        cpus, success = parse_expr(value, self.task, self.cmdtool)
        # if cpus is not None and not isinstance(cpus, j.Selector) and not isinstance(cpus, (int, float)):
        if isinstance(cpus, str):
            cpus = int(cpus)
        return cpus
    
    def fallback(self) -> Optional[int]:
        msg = 'Error parsing cpus requirement'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return None


class TaskMemoryParser(TaskParser):

    def do_parse(self) -> Optional[float]:
        value = self.task.runtime.get("memory")
        if value is None:
            return None
        s, success = parse_expr(value, self.task, self.cmdtool)
        if s is None:
            return None
        elif isinstance(s, str):
            if s.lower().endswith("g"):
                return float(s[:-1].strip())
            if s.lower().endswith("gb"):
                return float(s[:-2].strip())
            elif s.lower().endswith("gib"):
                return float(s[:-3].strip()) * 1.074
            elif s.lower().endswith("mb"):
                return float(s[:-2].strip()) / 1000
            elif s.lower().endswith("mib"):
                return float(s[:-3].strip()) / 1024
            msg = f'Error parsing memory requirement from string: {s}'
            log_message(self.cmdtool.uuid, msg, category=ErrorCategory.SCRIPTING)
            raise RuntimeError
        elif isinstance(s, (float, int)):
            # in bytes?
            return s / (1024 ** 3)
        elif isinstance(s, j.Selector):
            return s
        msg = f"Couldn't recognise memory requirement '{value}'"
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.SCRIPTING)
        raise RuntimeError
    
    def fallback(self) -> Optional[float]:
        msg = 'Error parsing memory requirement'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return None


class TaskDiskParser(TaskParser):

    def do_parse(self) -> Optional[float]:
        value = self.task.runtime.get("disks")
        if value is None:
            return None
        s, success = parse_expr(value, self.task, self.cmdtool)
        if s is None:
            return None
        if isinstance(s, str):
            try:
                return int(s)
            except ValueError:
                pass
            pattern_matcher = re.match(r"local-disk (\d+) .*", s)
            if not pattern_matcher:
                msg = f"Couldn't recognise disk type '{value}'"
                log_message(self.cmdtool.uuid, msg, category=ErrorCategory.SCRIPTING)
                raise RuntimeError
            s = pattern_matcher.groups()[0]
            try:
                return int(s)
            except ValueError:
                pass
            if s.lower().endswith("gb"):
                return float(s[:-2].strip())
            elif s.lower().endswith("gib"):
                return float(s[:-3].strip()) * 1.074
            elif s.lower().endswith("mb"):
                return float(s[:-2].strip()) / 1000
            elif s.lower().endswith("mib"):
                return float(s[:-3].strip()) / 1024
            raise Exception(f"Disk type type {s}")
        elif isinstance(s, (float, int)):
            # in GiB
            return s * 1.07374
        elif isinstance(s, j.Selector):
            return s
        elif s is None:
            return 2.14748  # 2 GiB
        msg = f"Couldn't recognise memory requirement '{value}'"
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.SCRIPTING)
        raise RuntimeError
    
    def fallback(self) -> Optional[float]:
        msg = 'Error parsing disk requirement'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return None


class TaskEnvVarsParser(TaskParser):

    def do_parse(self) -> dict[str, Any]:
        # TODO???
        return {}
    
    def fallback(self) -> dict[str, Any]:
        msg = 'Error parsing environment variables'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return {}


class TaskFilesToCreateParser(TaskParser):

    def do_parse(self) -> dict[str, Any]:
        # TODO???
        return {}
    
    def fallback(self) -> dict[str, Any]:
        msg = 'Error parsing files to create'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return {}


class TaskDirsToCreateParser(TaskParser):

    def do_parse(self) -> list[str | Selector]:
        # TODO???
        return []
    
    def fallback(self) -> list[str | Selector]:
        msg = 'Error parsing directories to create'
        log_message(self.cmdtool.uuid, msg, category=ErrorCategory.FALLBACKS)
        return []
