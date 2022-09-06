

from galaxy2janis.gx.command.components import Positional
from galaxy2janis.gx.command.components import Flag
from galaxy2janis.gx.command.components import Option
from galaxy2janis.gx.command.components import RedirectOutput
from galaxy2janis.gx.command.components.outputs.InputOutput import InputOutput
from galaxy2janis.gx.command.components.outputs.WildcardOutput import WildcardOutput


from .mock_params import MOCK_DATAPARAM1
from .mock_params import MOCK_BOOLPARAM1
from .mock_params import MOCK_FLOATPARAM1
from .mock_params import MOCK_SELECTPARAM1
from .mock_params import MOCK_OUTPARAM1

from .mock_params import (
    MOCK_BOOLPARAM1,
    MOCK_SELECTPARAM1,
    MOCK_FLOATPARAM1,
    MOCK_OUTPARAM1
)

MOCK_POSIT2 = Positional()
MOCK_POSIT2.values.add('abricate')
MOCK_POSIT2.cmd_pos = 0
MOCK_POSIT2.before_opts = True




MOCK_POSITIONAL1 = Positional()
MOCK_POSITIONAL1.gxparam = MOCK_DATAPARAM1
MOCK_POSITIONAL1.values.add('$in_fastq')
MOCK_POSITIONAL1.cmd_pos = 1
MOCK_POSITIONAL1.before_opts = True

MOCK_FLAG1 = Flag(prefix='--noheader')
MOCK_FLAG1.cmd_pos = 2
MOCK_FLAG1.gxparam = MOCK_BOOLPARAM1

MOCK_OPTION1 = Option(prefix='--minid')
MOCK_OPTION1.values.add('$adv.min_dna_id')
MOCK_OPTION1.delim = '='
MOCK_OPTION1.cmd_pos = 2
MOCK_OPTION1.forced_array = True
MOCK_OPTION1.gxparam = MOCK_FLOATPARAM1

MOCK_OPTION2 = Option(prefix='--db')
MOCK_OPTION2.values.add('card')
MOCK_OPTION2.delim = '='
MOCK_OPTION2.cmd_pos = 2
MOCK_OPTION2.gxparam = MOCK_SELECTPARAM1

MOCK_REDIRECT_OUTPUT = RedirectOutput('>', 'report.txt')
MOCK_REDIRECT_OUTPUT.gxparam = MOCK_OUTPARAM1

MOCK_WILDCARD_OUTPUT = WildcardOutput()
MOCK_WILDCARD_OUTPUT.gxparam = MOCK_OUTPARAM1

MOCK_INPUT_OUTPUT = InputOutput(MOCK_POSITIONAL1)

