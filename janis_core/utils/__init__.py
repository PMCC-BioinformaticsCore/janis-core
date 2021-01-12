from typing import Dict, Any, List, Tuple, TypeVar, Optional

T = TypeVar("T")


def first_value(d: Dict[Any, T]) -> T:
    # if isinstance(d, list):
    #     return d[0]
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


def fully_qualify_filename(fn):
    """
    The theory is, if the user types in a relative path (from the cwd), we should fully qualify this path.
    We'd also want to resolve `~` / `.` and other operators too.
    :param fn:
    :return:
    """
    from re import compile
    import os.path

    uri_prefix = compile("^[A-z0-9]{2,}:\/\/")

    if fn is None:
        return None
    if isinstance(fn, list):
        return [fully_qualify_filename(f) for f in fn]
    if uri_prefix.match(fn):
        return fn
    return os.path.abspath(os.path.expanduser(os.path.expandvars(fn)))


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


def is_module_available(module_name):
    import sys

    torch_loader = None
    if sys.version_info < (3, 0):
        # python 2
        import importlib

        torch_loader = importlib.find_loader(module_name)
    elif sys.version_info <= (3, 3):
        # python 3.0 to 3.3
        import pkgutil

        torch_loader = pkgutil.find_loader(module_name)
    elif sys.version_info >= (3, 4):
        # python 3.4 and above
        import importlib

        torch_loader = importlib.util.find_spec(module_name)

    return torch_loader is not None


def find_duplicates(ar) -> List:
    counts = {}
    for x in ar:
        counts[x] = (counts[x] if x in counts else 0) + 1

    return list(k for k, v in counts.items() if v > 1)


def lowercase_dictkeys(d: Optional[Dict]) -> Optional[Dict]:
    if not d:
        return None

    return {k.lower(): v for k, v in d.items()}


def generate_cat_command_from_statements(path, contents):
    wrap_tags = "EOT"
    potential_tags = ["EOF", "ENDOFFILE", "ENDOFTHISFILE"]
    while wrap_tags in contents:
        # generate new END tag
        if len(potential_tags) == 0:
            raise Exception(
                "Couldn't determine UNIQUE start / end tags for CAT <<{tag} >> $PATH {{contents}} {tag}"
            )
        wrap_tags = potential_tags.pop(0)

    return f"""\
cat <<{wrap_tags} >> '{path}'
{contents}
{wrap_tags}\
"""
