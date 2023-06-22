
from typing import Any

from .blocks import get_next_block 
from .blocks import CheetahBlock
from .. import utils


def sectional_evaluate(text: str, inputs: dict[str, Any]) -> str:
    raw_lines = utils.split_lines_blanklines(text)
    evaluator = PartialCheetahEvaluator(raw_lines, inputs)
    eval_lines = evaluator.evaluate()
    return utils.join_lines(eval_lines)


class PartialCheetahEvaluator:
    def __init__(self, lines: list[str], input_dict: dict[str, Any]):
        self.lines = lines
        self.input_dict = input_dict
        self.ptr: int = 0

    def evaluate(self) -> list[str]:
        try:
            return self.evaluation_worker()
        except Exception as e:
            return self.lines

    def evaluation_worker(self) -> list[str]:
        while self.ptr < len(self.lines):
            block = get_next_block(self.ptr, self.lines)
            block.evaluate(self.input_dict)
            self.update_lines(block)
            self.update_ptr(block)
        return self.lines

    def update_lines(self, block: CheetahBlock) -> None:
        if block.evaluated:
            old_lines_top = self.lines[:block.start]
            old_lines_bottom = self.lines[block.stop + 1:]
            self.lines = old_lines_top + block.lines + old_lines_bottom
    
    def update_ptr(self, block: CheetahBlock) -> None:
        if block.evaluated:
            self.ptr += 1  # success, go to next line
        else:
            self.ptr += block.height 




