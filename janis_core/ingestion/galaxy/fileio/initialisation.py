
import os


### MISC ###

def _make_parent_folders(path: str) -> None:
    heirarchy = _gen_heirarchy(path)
    for folder in heirarchy[:-1]:
        _safe_init_folder(folder)

def _gen_heirarchy(path: str) -> list[str]:
    path = os.path.normpath(path)
    path_list = path.split(os.sep)

    out: list[str] = []
    for i, elem in enumerate(path_list):
        parent = os.sep.join(path_list[:i])
        if parent == '':
            out.append(elem)
        else:
            out.append(f'{parent}{os.sep}{elem}')
    return out


### FILES ###

def init_file(path: str, override: bool=False, contents: str='') -> None:
    _make_parent_folders(path)
    if override:
        _unsafe_init_file(path, contents)
    else:
        _safe_init_file(path, contents)

def _safe_init_file(path: str, contents: str) -> None:
    if not os.path.exists(path):
        with open(path, 'w') as fp:
            fp.write(contents)
    
def _unsafe_init_file(path: str, contents: str) -> None:
    if os.path.exists(path):
        os.remove(path)
    with open(path, 'w') as fp:
        fp.write(contents)


### FOLDERS ###

def init_folder(path: str, override: bool=False) -> None:
    _make_parent_folders(path)
    if override:
        _unsafe_init_folder(path)
    else:
        _safe_init_folder(path)

def _safe_init_folder(path: str) -> None:
    if not os.path.isdir(path):
        os.mkdir(path)

def _unsafe_init_folder(path: str) -> None:
    if os.path.isdir(path):
        os.remove(path)
    os.mkdir(path)

