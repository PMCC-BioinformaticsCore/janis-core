


from typing import Optional
from dataclasses import dataclass


def get_id_filename(identifier: str) -> str:
    cwl_ref = get_cwl_reference(identifier)
    assert(cwl_ref.filename)
    return cwl_ref.filename

def get_id_path(identifier: str) -> Optional[str]:
    cwl_ref = get_cwl_reference(identifier)
    return cwl_ref.internal_path

def get_id_entity(identifier: str) -> str:
    cwl_ref = get_cwl_reference(identifier)
    assert(cwl_ref.entity)
    return cwl_ref.entity

def remove_output_name_from_output_source(identifier: str) -> str:
    assert('#' in identifier)
    block1, block2 = identifier.split('#')
    block2 = block2.split('/', 1)[1]
    return f'{block1}#{block2}'



@dataclass
class CWLReference:
    filename: Optional[str] # myfile.cwl
    internal_path: Optional[str]   # a step
    entity: Optional[str]   # a step or workflow input or tool input


def get_cwl_reference(identifier: str) -> CWLReference:
    identifier = identifier.replace('-', '_')
    
    if identifier.startswith('_:'):
        filename = None
        internal_path = None
        entity = None
    
    elif '#' in identifier:
        block1, block2 = identifier.split('#')
        filename = block1.split('.')[0].split('/')[-1]
        if '/' in block2:
            internal_path, entity = block2.split('/')[-2:]
            # if block2.count('/') == 1:
            #     internal_path, entity = block2.split('/')
            # elif block2.count('/') == 2:
            #     internal_path, entity = block2.split('/')[1:]
            # else:
            #     raise NotImplementedError
        else:
            internal_path = None
            entity = block2
    elif '.cwl' in identifier:
        filename = identifier.split('.')[0].split('/')[-1]
        internal_path = None
        entity = None
    else:
        raise NotImplementedError

    return CWLReference(filename, internal_path, entity)





