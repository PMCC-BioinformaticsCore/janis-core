

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
[a-z]               ->      [a-z]
[^a-m]              ->      [!a-m]
[abc]               ->      [a,b,c]
(cat|dog|bat)       ->      {cat,dog,bat}
(.*\.tar|.*\.gz)    ->      {*.tar,*.gz}


"""



