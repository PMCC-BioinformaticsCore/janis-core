
import os


def safe_init_file(path: str, override: bool=False, contents: str='') -> None:
    dirname = os.path.dirname(path)
    safe_init_folder(dirname)
    with open(path, 'w') as fp:
        fp.write(contents)

def safe_init_folder(path: str, override: bool=False) -> None:
    if not os.path.isdir(path):
        os.makedirs(path)

