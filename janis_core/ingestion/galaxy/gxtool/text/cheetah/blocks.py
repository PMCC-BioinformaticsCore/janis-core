

from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional
from uuid import uuid4
from enum import Enum, auto
from Cheetah.Template import Template
from galaxy.util import unicodify

from .. import utils
from janis_core.ingestion.galaxy.gxtool.command.cmdstr import constructs

from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.expressions.patterns import CHEETAH_EDGE_CASE_INPUT


def get_blocks(ptr: int, lines: list[str], indent_level: int) -> list[CheetahBlock]:
    block_lines = LineFactory(lines).get_lines()
    factory = BlockFactory(ptr, block_lines, indent_level)
    return factory.get_blocks()

def get_next_block(ptr: int, lines: list[str]) -> CheetahBlock:
    blocks = get_blocks(ptr=ptr, lines=lines[ptr:], indent_level=0)
    return blocks[0]


@dataclass
class BlockLine:
    line_num: int
    indent: int
    text: str

class LineFactory:
    """creates BlockLines from normal text lines"""

    def __init__(self, lines: list[str]):
        self.text_lines = lines
        self.block_lines: list[BlockLine] = []
        self.construct_tracker: constructs.ConstructTracker = constructs.ConstructTracker()

    def get_lines(self) -> list[BlockLine]:
        if len(self.text_lines) == 0:
            return []
        else:
            for i, line in enumerate(self.text_lines):
                self.add_block_line(i, line)
            return self.block_lines

    def add_block_line(self, i: int, line: str) -> None:
        indent = self.construct_tracker.stack.depth
        self.construct_tracker.update(line) # LEAVE THIS HERE
        block_line = BlockLine(
            line_num=i,
            indent=indent,
            text=line
        )
        self.block_lines.append(block_line)
    
    


IGNORE_CONSTRUCTS = constructs.CH_WITHIN_CONDITIONAL | constructs.CH_CLOSE_CONDITIONAL | constructs.CH_CLOSE_LOOP | constructs.CH_CLOSE_FUNC

class BlockFactory:
    """
    gets all cheetah blocks occurring in self.lines on a particular 
    indent level
    """
    def __init__(self, offset: int, block_lines: list[BlockLine], indent_level: int):
        self.offset = offset
        self.block_lines = block_lines
        self.target_indent = indent_level
        self.active_lines: list[BlockLine] = []
        self.blocks: list[CheetahBlock] = []

    def get_blocks(self) -> list[CheetahBlock]:
        for line in self.block_lines:
            if self.start_block(line): # new block begins
                self.add_block() # add previous block
                self.active_lines = [] # reset active lines
            if self.within_block(line):
                self.active_lines.append(line)
        self.add_block()
        self.filter_blocks()
        return self.blocks

    def start_block(self, line: BlockLine) -> bool:
        if line.indent == self.target_indent:
            return True
        return False
    
    def within_block(self, line: BlockLine) -> bool:
        if line.indent >= self.target_indent:
            return True
        return False

    def add_block(self) -> None:
        if self.active_lines:
            block = CheetahBlock(
                btype=self.select_block_type(),
                start=self.active_lines[0].line_num + self.offset,
                stop=self.active_lines[-1].line_num + self.offset,
                lines=[ln.text for ln in self.active_lines]
            )
            self.blocks.append(block)
    
    def select_block_type(self) -> BlockType:
        line = self.active_lines[0].text
        types = [
            (constructs.CH_OPEN_CONDITIONAL, BlockType.CONDITIONAL),
            (constructs.CH_OPEN_LOOP, BlockType.LOOP),
            (constructs.CH_OPEN_FUNC, BlockType.FUNCTION),
            (constructs.CH_MACRO, BlockType.INLINE_CH),
            (constructs.LINUX_ALIAS, BlockType.INLINE_ALIAS),
        ]
        for construct, block_type in types:
            if any ([line.startswith(text) for text in construct]):
                return block_type
        return BlockType.INLINE

    def filter_blocks(self):
        filtered: list[CheetahBlock] = []
        for block in self.blocks:
            if len(block.lines) == 1:
                line = block.lines[0]
                if not any ([line.startswith(text) for text in IGNORE_CONSTRUCTS]):
                    filtered.append(block)
            else:
                filtered.append(block)
        self.blocks = filtered



class EvaluationStrategy(ABC):
    def __init__(self, lines: list[str], input_dict: dict[str, Any]):
        self.lines = lines
        self.input_dict = input_dict
    
    def eval(self) -> Optional[list[str]]:
        template = self.prepare_template()
        outcome = self.evaluate_template(template)
        if outcome is not None: 
            return self.handle_outcome(outcome)
        return None

    @abstractmethod
    def prepare_template(self) -> list[str]:
        """prepares the template text for evaluation"""
        ...
    
    def evaluate_template(self, source_lines: list[str]) -> Optional[list[str]]:
        """performs cheetah evaluation of template"""
        try:
            source = utils.join_lines(source_lines)
            t = Template(source, searchList=[self.input_dict]) # type: ignore
            evaluation = str(unicodify(t))
            return utils.split_lines(evaluation)
        except:
            return None

    @abstractmethod
    def handle_outcome(self, outcome: list[str]) -> list[str]:
        """handles the evaluated text (if successful) and applies any transformations needed"""
        ...



class InlineEvaluationStrategy(EvaluationStrategy):

    def prepare_template(self) -> list[str]:
        assert(len(self.lines) == 1)
        return self.lines

    def handle_outcome(self, outcome: list[str]) -> list[str]:
        return outcome


class ConditionalEvaluationStrategy(EvaluationStrategy):
    immune_lines = '__JANIS_MAIN__'

    def prepare_template(self) -> list[str]:
        self.mask_children()
        return self.lines

    def handle_outcome(self, outcome: list[str]) -> list[str]:
        output = self.create_blank_output()
        output = self.restore_surviving_children(outcome, output)
        return output

    def mask_children(self) -> None:
        """prepares text ready for evaluation. swaps child blocks with identifiers"""
        self.masked_blocks: dict[str, CheetahBlock] = {}
        for child in self.get_child_blocks():
            self.masked_blocks[child.uuid] = child  # register identifier/block
            self.substitute_identifier(child)

    def get_child_blocks(self) -> list[CheetahBlock]:
        """returns next level of cheetah blocks within this block"""
        return get_blocks(
            ptr=0,
            lines=self.lines, 
            indent_level=1
        )

    def substitute_identifier(self, block: CheetahBlock) -> None:
        """swaps block block with identifier in template"""
        old_lines_top = self.lines[:block.start]
        old_lines_bottom = self.lines[block.stop+1:]
        self.lines = old_lines_top + [block.uuid] * block.height + old_lines_bottom

    def create_blank_output(self) -> list[str]:
        return [''] * len(self.lines)

    def restore_surviving_children(self, evaluation: list[str], output: list[str]) -> list[str]:
        surviving_children = self.get_surviving_children(evaluation)
        for identifier, block in surviving_children.items():
            line_num = self.get_identifier_line_num(self.lines, identifier)
            if line_num is not None:
                block = self.masked_blocks[identifier]
                output = self.substitute_block(output, block)
        return output

    def get_surviving_children(self, evaluation: list[str]) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for identifier, block in self.masked_blocks.items():
            if block.lines[0] in self.immune_lines:
                out[identifier] = block
            else:
                position = self.get_identifier_line_num(evaluation, identifier)
                if position is not None:
                    out[identifier] = block
        return out

    def get_identifier_line_num(self, lines: list[str], identifier: str) -> Optional[int]:
        for i, line in enumerate(lines):
            if line == identifier:
                return i
        return None

    def substitute_block(self, output: list[str], block: CheetahBlock) -> list[str]:
        """swaps block block with identifier in template"""
        old_lines_top = output[:block.start]
        old_lines_bottom = output[block.stop+1:]
        return old_lines_top + block.lines + old_lines_bottom



class BlockType(Enum):
    MAIN            = auto()
    INLINE          = auto()
    INLINE_CH       = auto()
    INLINE_ALIAS    = auto()
    CONDITIONAL     = auto()
    FUNCTION        = auto()
    LOOP            = auto()


class CheetahBlock:
    """
    The smallest unit of evaluatable cheetah logic
    examples: single line statement, conditional block, loop block
    """
    def __init__(self, btype: BlockType, start: int, stop: int, lines: list[str]):
        self.btype = btype
        self.start = start
        self.stop = stop
        self.lines = lines
        self.uuid: str = str(uuid4())
        self.evaluated: bool = False
    
    @property 
    def height(self) -> int:
        return self.stop - self.start + 1
    
    def evaluate(self, input_dict: dict[str, Any]) -> None:
        self.will_evaluate = self.should_evaluate(input_dict)
        if self.will_evaluate:
            strategy = self.get_eval_strategy(input_dict)
            evaluation = strategy.eval()
            if evaluation is not None:
                assert(len(self.lines) == len(evaluation))
                self.evaluated = True
                self.lines = evaluation
        
    def should_evaluate(self, input_dict: dict[str, Any]) -> bool:
        """dictates whether this block should be evaluated or left as original text"""
        permitted_blocks = [BlockType.INLINE, BlockType.INLINE_ALIAS, BlockType.INLINE_CH, BlockType.CONDITIONAL]
        if self.lines == ['']:
            return False
        elif self.btype not in permitted_blocks:
            return False
        elif self.edge_case_input(input_dict):
            return False
        return True

    def edge_case_input(self, input_dict: dict[str, Any]) -> bool:
        for line in self.lines:
            matches = expressions.get_matches(line, CHEETAH_EDGE_CASE_INPUT)
            if matches:
                if 'input' not in input_dict:
                    return True
        return False

    def get_eval_strategy(self, input_dict: dict[str, Any]) -> EvaluationStrategy:
        if self.btype in [BlockType.INLINE, BlockType.INLINE_ALIAS, BlockType.INLINE_CH]:
            return InlineEvaluationStrategy(self.lines, input_dict)
        elif self.btype == BlockType.CONDITIONAL:
            return ConditionalEvaluationStrategy(self.lines, input_dict)
        else:
            raise RuntimeError()


