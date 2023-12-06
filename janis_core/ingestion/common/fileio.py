

import os 
import shutil
PERMISSIONS=0o777

def safe_init_file(path: str, override: bool=False, contents: str='') -> None:
    dirname = os.path.dirname(path)
    safe_init_folder(dirname)
    with open(path, 'w') as fp:
        fp.write(contents)

def safe_init_folder(path: str, override: bool=False) -> None:
    if override:
        if os.path.isdir(path):
            shutil.rmtree(path)
    os.makedirs(path, PERMISSIONS, exist_ok=True)