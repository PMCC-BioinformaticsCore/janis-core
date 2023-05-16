

"""
glob            regex	            Description

?	            .	                Any single character
*	            .*	                Zero or more characters
[a-z]	        [a-z]	            Any character from the range
[!a-m]	        [^a-m]	            A character not in the range
[a,b,c]	        [abc]	            One of the given characters
{cat,dog,bat}	(cat|dog|bat)	    One of the given options
{*.tar,*.gz}	(.*\.tar|.*\.gz)	One of the given options, considering nested wildcards

regex                       glob                
\.                  ->      .
.                   ->      ?
.*                  ->      *
.+                  ->      *
[a-z]               ->      [a-z]
[^a-m]              ->      [!a-m]
[abc]               ->      [a,b,c]
(cat|dog|bat)       ->      {cat,dog,bat}
(.*\.tar|.*\.gz)    ->      {*.tar,*.gz}


galaxy specific 
__name_and_ext__    ->      *.*

"""

# LOGICAL_OR      = r"\([^|()]+(\|[^|()]+)+\)"    # (cat|dog|bat) -> {cat,dog,bat}
# CHAR_SET        = r"\[(\^?)[^|()]+?\]"          # [^a-z134]     -> [!a-z,1,2,3]
#                                                 # (capturing group 1 identifies whether negated set)
# ZERO_OR_MORE    = r"((\[[^|()\[\]]+\])|(\([^|()]+(\|[^|()]+)+\))|(\\.))\*\??"   # [^a-z134]* -> *, (cat|dog)* -> *
# ONE_OR_MORE     = r"((\[[^|()\[\]]+\])|(\([^|()]+(\|[^|()]+)+\))|(\\.))\+\??"   # [^a-z134]+ -> *, (cat|dog)+ -> *
# SPECIAL_CHARS   = r"(?<!\\)(\\S|\\s|\.|\\c|\\d|\\D|\\w|\\W|\\x|\\O|\\A|\$|\^|\\Z|\\b|\\B|\\<|\\>)(?![*+])" # \S -> ?, \S+ doesn't match

SPECIAL_CHARS   = r"(?<!\\)(\\S|\\s|\.|\\c|\\d|\\D|\\w|\\W|\\x|\\O|\\A|\$|\^|\\Z|\\b|\\B|\\<|\\>)" # \S -> ?, \S+ doesn't match
LOGICAL_OR      = r"\(([^|()]+(\|[^|()]+)+)\)"    # (cat|dog|bat) -> {cat,dog,bat}
CHAR_SET        = r"\[((\^?)[^|()]+?)\]"          # [^a-z134]     -> [!a-z,1,2,3]
                                                # (capturing group 1 identifies whether negated set)
ZERO_OR_MORE    = fr"(({SPECIAL_CHARS})|({LOGICAL_OR})|({CHAR_SET}))\*\??"   # [^a-z134]* -> *, (cat|dog)* -> *
ONE_OR_MORE     = fr"(({SPECIAL_CHARS})|({LOGICAL_OR})|({CHAR_SET}))\+\??"   # [^a-z134]+ -> *, (cat|dog)+ -> *
BRACKETS        = r"(\(.+?\))|(\[.+?\])|(\{.+?\})"

"""
ORDER:

- ZERO_OR_MORE
- ONE_OR_MORE
- LOGICAL_OR
- CHAR_SET
- SPECIAL_CHARS
- remove escape characters

"""






