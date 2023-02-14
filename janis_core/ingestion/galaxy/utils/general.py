

from Bio import pairwise2


def global_align(pattern: str, template: str) -> int:
    pattern = pattern.lower()
    template = template.lower()
    outcome = pairwise2.align.globalms(pattern, template, 2, -1, -.5, -.1) # type: ignore
    if len(outcome) > 0: # type: ignore
        score = outcome[0].score # type: ignore
    else:
        score = 0
    return score # type: ignore


