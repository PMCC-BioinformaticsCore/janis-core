
ALL = r'.*'

WORD = r'(\'.*?(?<!\\)\'[^\s]*)|(".*?(?<!\\)"[^\s]*)|([^\s]+)'
STRING  = r'^([\/\\\w\d-.*`@])[-\w\d\{\}\$.\/\\_:*`@]*$'
SIMPLE_STRING = r'[\w$_-]+'
EMPTY_STRING = r'\'\'|""'

BACKTICK_SECTION = r'`.+?`'
QUOTED_SECTION = r'"([^\"]*?)"|\'([^\']*?)\''
QUOTED_SECTION_W_NEWLINE = r'\'[^\']*?\n[^\']*?\'|"[^"]*?\n[^"]*?"'

INTEGER = r'^-?\d+$'
FLOAT   = r'^(-?(\.\d+)|-?(\d+\.\d+))(e-?\d+)?$'
NUMBER  = r'^-?\.?\d+(\.\d+)?$'
OPERATOR = r'^[-+\\/*=]?=$'
VERSION = r'^(\d+)(\.\d+)*$'

LINUX_LN = r'(?:\s|^)ln(?:\s-[sf]+)* [\'"]?([-\{}$./\\*\w\d]+)[\'"]? [\'"]?([-\{}$./\\*\w\d]+)[\'"]?(?=\s|$)'
LINUX_MV = r'(?:\s|^)mv(?:\s-f+)* ([-\'"{}$./\\*\w\d]+) ([-\'"{}$./\\*\w\d]+)(?=\s|$)'
LINUX_CP = r'(?:\s|^)cp(?:\s-[rpfR]+)* ([-_\'"{}$.*/\\\w\d]+) ([-_\'"{}$.*/\\\w\d]+)(?=\s|$)'
LINUX_STATEMENT_DELIMS = r'(?<!\\)(&&|\|?\|(?! tee |tee ))(?=\s|$)' 
LINUX_REDIRECT = r'((?<=\s)\d|&)?>[>&]?(?![>&]?\d)'
LINUX_TEE = r'(?<![\d&])\| ?tee( -a)?'
LINUX_STREAM_MERGE = r'(?<=\s|^)\d?>&\d'

KEYVAL_PAIR = r'(?<=\s|^)(\S+?)([=:])(\S+?)(?=\s|$)'
COMPOUND_OPT = r'^(-\w)(\d+?)$'

CHEETAH_SET = r'(?:^|\s)#set ([$\w\d]+) = [\'"]?([$\w\d.]+)[\'"]?(?=\s|$)'
CHEETAH_EDGE_CASE_INPUT = r'\${?input([ =.}\'")]|$)'

VARIABLES_FMT1 = r'\$\w[\w._]+'
VARIABLES_FMT2 = r'\$\{\w[\w._]+\}'
FUNCTION_CALL_FMT1 = r'\$\{[^(].+?(\(.*\))[^(]*\}'
FUNCTION_CALL_FMT2 = r'\$[^(){} \n\'"]+(\(.*\))[^(){} \n\'"]*'

SCRIPT = r'(\$__tool_directory__\/)([^\s\$]+)'
GX_DYNAMIC_KEYWORDS = r'\\?\${[\w\d]+\:\-(\d+)}'
GX_STATIC_KEYWORDS = r'\$__new_file_path__|\$__tool_data_path__|\$__root_dir__|\$__datatypes_config__|\$__user_id__|\$__user_email__|\$__app__|\$__target_datatype__'

WILDCARD_GROUP = r'\((\?P<.+?>)(.*?)\)'

