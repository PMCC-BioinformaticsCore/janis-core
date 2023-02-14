
# refgenconf needs the following aliases to be done manually.
# setting these aliases here is global scope. 
# this allows refgenconf==0.9.3 to be used in python 3.10.
import collections.abc
collections.Iterable = collections.abc.Iterable
collections.Mapping = collections.abc.Mapping
