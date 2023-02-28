


from typing import Optional
from dataclasses import dataclass


"""
file reference
'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl'

workflow input reference
'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl#annotation_file'

step reference
'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl#assign_genes'

step output reference
file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl#assign_genes/result_file'

internal clt definition
'_:637c1ad6-a30f-4c92-80fe-b9f1448c02bd'

clt input reference
'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl#assign_genes/annotation_filename'

clt output reference
'file:///home/grace/work/pp/translation/janis-assistant/janis_assistant/tests/data/cwl/super_enhancer_wf.cwl#assign_genes/result_file'
"""


@dataclass
class CWLReference:
    filename: Optional[str] # myfile.cwl
    prefix: Optional[str]   # a step
    suffix: Optional[str]   # a step or workflow input or tool input


def get_cwl_reference(identifier: str) -> CWLReference:
    if identifier.startswith('_:'):
        filename = None
        prefix = None
        suffix = None
    
    elif '#' in identifier:
        block1, block2 = identifier.split('#')
        filename = block1.split('.')[0].split('/')[-1]
        if '/' in block2:
            prefix, suffix = block2.split('/')
        else:
            prefix = None
            suffix = block2
    else:
        filename = identifier.split('.')[0].split('/')[-1]
        prefix = None
        suffix = None
        
    return CWLReference(filename, prefix, suffix)

