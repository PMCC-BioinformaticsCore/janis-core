

from .tool import XMLTool
from .citations import XMLCitation
from .configfile import XMLConfigfile
from .metadata import XMLMetadata
from .script import XMLScript

from .requirement import XMLRequirement
from .requirement import XMLCondaRequirement
from .requirement import XMLContainerRequirement

from .params.param import XMLParam
from .params.param_register import XMLParamRegister

from .params.input_param import XMLInputParam
from .params.input_param import XMLSelectOption
from .params.input_param import XMLTextParam
from .params.input_param import XMLIntegerParam
from .params.input_param import XMLFloatParam
from .params.input_param import XMLBoolParam
from .params.input_param import XMLSelectParam
from .params.input_param import XMLDataParam
from .params.input_param import XMLDataCollectionParam

from .params.output_param import XMLOutputParam
from .params.output_param import XMLDataOutputParam
from .params.output_param import XMLCollectionOutputParam

from .tests import XMLTest
from .tests import XMLTestRegister

