
import WDL 
import janis_core as j
from typing import Optional, Any
import regex as re 

from .expressions import parse_expr

DEFAULT_CONTAINER = "ubuntu:latest"


def parse_container_requirement(task: WDL.Tree.Task) -> str:
    container = task.runtime.get("container", task.runtime.get("docker"))
    if isinstance(container, WDL.Expr.Get):
        # relevant input
        inp = [i.expr for i in task.inputs if i.name == str(container.expr)]
        if len(inp) > 0:
            container = inp[0]
        else:
            j.Logger.warn(
                f"Expression for determining containers was '{container}' "
                f"but couldn't find input called {str(container.expr)}"
            )
    if isinstance(container, WDL.Expr.String):
        container = container.literal
    if isinstance(container, WDL.Value.String):
        container = container.value
    if container is None:
        container = DEFAULT_CONTAINER
    if not isinstance(container, str):
        j.Logger.warn(
            f"Expression for determining containers ({container}) are not supported in Janis, using ubuntu:latest"
        )
        container = DEFAULT_CONTAINER
    return container

def parse_cpus_requirement(task: WDL.Tree.Task) -> int:
    value = task.runtime.get("cpu")
    cpus = parse_expr(value)
    # if cpus is not None and not isinstance(cpus, j.Selector) and not isinstance(cpus, (int, float)):
    if isinstance(cpus, str):
        cpus = int(cpus)

def parse_memory_requirement(task: WDL.Tree.Task) -> int:
    value = task.runtime.get("memory")
    s = parse_expr(value)
    if s is None:
        return 1.074
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
        raise Exception(f"Memory type {s}")
    elif isinstance(s, (float, int)):
        # in bytes?
        return s / (1024 ** 3)
    elif isinstance(s, j.Selector):
        return s
    raise Exception(f"Couldn't recognise memory requirement '{value}'")

def parse_disk_requirement(task: WDL.Tree.Task) -> int:
    value = task.runtime.get("disks")
    s = parse_expr(value)
    if s is None:
        return None
    if isinstance(s, str):
        try:
            return int(s)
        except ValueError:
            pass
        pattern_matcher = re.match(r"local-disk (\d+) .*", s)
        if not pattern_matcher:
            raise Exception(f"Couldn't recognise disk type '{value}'")
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
    raise Exception(f"Couldn't recognise memory requirement '{value}'")