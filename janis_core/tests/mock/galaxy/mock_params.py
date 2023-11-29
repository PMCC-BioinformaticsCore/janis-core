

from janis_core.ingestion.galaxy.gxtool.model import XMLBoolParam
from janis_core.ingestion.galaxy.gxtool.model import XMLFloatParam
from janis_core.ingestion.galaxy.gxtool.model import XMLSelectOption
from janis_core.ingestion.galaxy.gxtool.model import XMLSelectParam
from janis_core.ingestion.galaxy.gxtool.model import XMLDataParam
from janis_core.ingestion.galaxy.gxtool.model import XMLParamRegister
from janis_core.ingestion.galaxy.gxtool.model import XMLDataOutputParam


MOCK_DATAPARAM1 = XMLDataParam('file_input')
MOCK_DATAPARAM1.formats = ['fastq']
MOCK_DATAPARAM1.multiple = False
MOCK_DATAPARAM1.helptext = "input fastq"

MOCK_BOOLPARAM1 = XMLBoolParam('adv.no_header')
MOCK_BOOLPARAM1.argument = '--noheader'
MOCK_BOOLPARAM1.checked = False
MOCK_BOOLPARAM1.truevalue = '--noheader'
MOCK_BOOLPARAM1.falsevalue = ''
MOCK_BOOLPARAM1.helptext = "Suppress output file's column headings"
MOCK_BOOLPARAM1.label = 'Suppress header'

MOCK_FLOATPARAM1 = XMLFloatParam('adv.min_dna_id')
MOCK_FLOATPARAM1.argument = '--minid'
MOCK_FLOATPARAM1.helptext = ''
MOCK_FLOATPARAM1.label = 'Minimum DNA %identity'
MOCK_FLOATPARAM1.value = '80'
MOCK_FLOATPARAM1.min = '0.0'
MOCK_FLOATPARAM1.max = '100.0'

MOCK_FLOATPARAM2 = XMLFloatParam('adv.min_cov')
MOCK_FLOATPARAM2.argument = '--mincov'
MOCK_FLOATPARAM2.helptext = ''
MOCK_FLOATPARAM2.label = 'Minimum DNA %coverage'
MOCK_FLOATPARAM2.value = '80'
MOCK_FLOATPARAM2.min = '0.0'
MOCK_FLOATPARAM2.max = '100.0'

MOCK_SELECTPARAM1 = XMLSelectParam('adv.db')
MOCK_SELECTPARAM1.argument = '--db'
MOCK_SELECTPARAM1.helptext = 'Option to switch to other AMR/other database'
MOCK_SELECTPARAM1.label = "Database to use - default is 'resfinder'"
MOCK_SELECTPARAM1.multiple = False
MOCK_SELECTPARAM1.options = [
    XMLSelectOption('argannot', False, 'ARG-ANNOT'),
    XMLSelectOption('card', False, 'CARD'),
    XMLSelectOption('ecoh', False, 'EcOH'),
    XMLSelectOption('ncbi', False, 'NCBI Bacterial Antimicrobial Resistance Reference Gene Database'),
    XMLSelectOption('resfinder', True, 'Resfinder'),
    XMLSelectOption('plasmidfinder', False, 'PlasmidFinder'),
    XMLSelectOption('vfdb', False, 'VFDB'),
    XMLSelectOption('megares', False, 'megares'),
    XMLSelectOption('ecoli_vf', False, 'Ecoli_VF')
]

MOCK_OUTPARAM1 = XMLDataOutputParam('report')
MOCK_OUTPARAM1.label = 'report file'
MOCK_OUTPARAM1.formats = ['txt']
MOCK_OUTPARAM1.from_work_dir = 'report.txt'

MOCK_PARAM_REGISTER = XMLParamRegister()

