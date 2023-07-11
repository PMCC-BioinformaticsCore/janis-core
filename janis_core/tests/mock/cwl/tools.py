import os
from cwl_utils.parser.cwl_v1_0 import load_document

CWL_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/cwl')

rose_tool = load_document(f'{CWL_TESTDATA_PATH}/tools/rose.cwl')