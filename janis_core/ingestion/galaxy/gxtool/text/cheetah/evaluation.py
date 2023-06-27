
from typing import Any
from collections import defaultdict

from .blocks import get_next_block 
from .blocks import CheetahBlock
from .. import utils


def sectional_evaluate(text: str, inputs: dict[str, Any]) -> str:
    raw_lines = utils.split_lines(text)
    evaluator = PartialCheetahEvaluator(raw_lines, inputs)
    eval_lines = evaluator.evaluate()
    return utils.join_lines(eval_lines)


class EvaluationMetrics:
    def __init__(self) -> None:
        self.total_blocks: int = 0
        self.success_blocks: int = 0
        self.block_types: dict[str, int] = defaultdict(int)

    def add(self, block: CheetahBlock) -> None:
        # block type metrics
        self.block_types[block.btype.name] += 1

        # evaluation metrics
        if block.will_evaluate:
            self.total_blocks += 1
            if block.evaluated:
                self.success_blocks += 1

    def report(self) -> None:
        print('--- CHEETAH EVAL ---')
        success_percent = self.success_blocks / self.total_blocks * 100
        print(f'success blocks: {success_percent:0.1f}%')
        for btype_name, count in self.block_types.items():
            print(f'{btype_name}: {count}')
        print()


class PartialCheetahEvaluator:
    def __init__(self, lines: list[str], input_dict: dict[str, Any]):
        self.lines = lines
        self.input_dict = input_dict
        self.ptr: int = 0
        self.metrics = EvaluationMetrics()

    def evaluate(self) -> list[str]:
        # do eval
        try:
            eval_lines = self.evaluation_worker()
        except Exception as e:
            print('critical cheetah templating error')
            eval_lines = self.lines

        # report metrics & return
        self.metrics.report()
        return eval_lines

    def evaluation_worker(self) -> list[str]:
        while self.ptr < len(self.lines):
            # do evaluation
            block = get_next_block(self.ptr, self.lines)
            block.evaluate(self.input_dict)
            # update metrics, original lines & line ptr
            self.metrics.add(block)
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




