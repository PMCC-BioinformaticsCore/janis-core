from typing import Dict, Any, List, Tuple


def first_value(d: Dict):
    return next(iter(d.values()))


def get_value_for_hints_and_ordered_resource_tuple(
    hints: Dict[str, Any], tuples: List[Tuple[str, Dict[str, int]]]
):
    if not hints:
        return None
    for k, d in tuples:
        if k not in hints:
            continue
        v = hints[k]
        if v not in d:
            continue
        return d[v]
    return None


def zip_directory(parent_dir, dir_name):
    import subprocess
    from .logger import Logger
    from os import chdir

    Logger.info("Zipping tools")
    chdir(parent_dir)

    zip_result = subprocess.run(["zip", "-r", f"{dir_name}.zip", f"{dir_name}/"])
    if zip_result.returncode == 0:
        Logger.info("Zipped tools")
    else:
        Logger.critical(zip_result.stderr)


def is_array_prefix(prefix, to):
    if len(prefix) > len(to):
        return False

    for i in range(len(prefix)):
        if prefix[i] != to[i]:
            return False

    return True


def recursive_2param_wrap(methodname, items):
    # It would be awesome if this is generic, but it opens up a new can of worms about how to wrap it
    # 2 might even be confusing, as you could say zip(zip(A, B), zip(C, D)) instead of zip(A, zip(B, zip(C, D)))
    # And more than length=2 you could have cases where the last function call might not have a clean 'n' params

    if len(items) < 1:
        raise Exception(
            "'recursive_2param_wrap' required two parameters to wrap a workflow"
        )
    if len(items) == 2:
        return f"{methodname}({items[0]}, {items[1]})"
    return f"{methodname}({items[0]}, {recursive_2param_wrap(methodname, items[1:])})"
