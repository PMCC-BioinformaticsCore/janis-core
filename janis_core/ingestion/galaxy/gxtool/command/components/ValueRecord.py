

from collections import defaultdict
from typing import Optional

from janis_core.ingestion.galaxy import expressions


class ValueRecord:
    def __init__(self):
        self.record: list[str] = []

    def add(self, value: str) -> None:
        self.record.append(value)
    
    @property
    def unique(self) -> list[str]:
        values = list(set([text for text in self.record]))
        values.sort()
        return values
    
    @property
    def counts(self) -> defaultdict[str, int]:
        counts: defaultdict[str, int] = defaultdict(int) 
        for text in self.record:
            if text != '': # TODO how would this happen ??
                counts[text] += 1
        return counts

    @property
    def most_common_value(self) -> Optional[str]:
        counts_dict = self.counts
        counts_list = list(counts_dict.items())
        counts_list.sort(key=lambda x: x[1], reverse=True)
        if len(counts_list) > 0:
            return counts_list[0][0]
        else:
            return None

    @property
    def env_var(self) -> Optional[str]:
        for text in self.record:
            if text.startswith('$'):
                return text
        return None
    
    # @property
    # def script(self) -> Optional[str]:
    #     for text in self.record:
    #         if expressions.is_script(text):
    #             matches = expressions.get_matches(text, expressions.patterns.SCRIPT)
    #             return matches[0].group(2)
    #     return None

    





