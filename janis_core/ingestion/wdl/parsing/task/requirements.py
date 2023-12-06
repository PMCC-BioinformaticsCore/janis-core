
import WDL 
import janis_core as j
import regex as re 
from typing import Any , Optional

from janis_core import Selector, StringFormatter, InputSelector
from janis_core.messages import log_message
from janis_core.messages import ErrorCategory

from ..expressions import parse_expr
from ..EntityParser import TaskParser


DEFAULT_CONTAINER = "ubuntu:latest"

class TaskContainerParser(TaskParser):

    def do_parse(self) -> Any:
        container = self.task.runtime.get("container", self.task.runtime.get("docker"))
        # from expr (NOTE: Janis doesn't really allow this. unsure if works for CWL / Nextflow)
        if isinstance(container, WDL.Expr.Base):
            res, success = parse_expr(container, self.task, self.cmdtool)
            if not success:
                raise RuntimeError
            if isinstance(res, StringFormatter):
                has_single_kwarg = True if len(res.kwargs) == 1 else False
                has_single_format = True if res._format.startswith('{') and res._format.endswith('}') else False
                if has_single_kwarg and has_single_format:
                    arg = list(res.kwargs.values())[0]
                    str_literal = self.get_string_literal_expr(arg)
                    if str_literal is not None:
                        container = str_literal
                    else:
                        container = str(arg)
            else:
                str_literal = self.get_string_literal_expr(res)
                if str_literal is not None:
                    container = str_literal
                else:
                    container = str(res)
        # from string literal
        elif isinstance(container, str):
            container = container
        # no container
        elif container is None:
            container = DEFAULT_CONTAINER
        else:
            raise NotImplementedError
        return container
    
    def get_string_literal_expr(self, sel: Selector) -> Optional[str]:
        if isinstance(sel, InputSelector):
            assert self.task.inputs is not None
            for inp in self.task.inputs:
                if isinstance(inp.expr, WDL.Expr.String):
                    if isinstance(inp.expr.literal, WDL.Value.String):
                        return inp.expr.literal.value
        return None

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
