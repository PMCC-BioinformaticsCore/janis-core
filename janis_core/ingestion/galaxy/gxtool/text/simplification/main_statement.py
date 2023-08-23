

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

from janis_core.ingestion.galaxy.gxtool.model import XMLTool
from janis_core.ingestion.galaxy.utils.general import global_align
from janis_core.ingestion.galaxy import expressions
from janis_core.ingestion.galaxy.expressions.patterns import (
    LINUX_STATEMENT_DELIMS,
    VARIABLES_FMT1,
    VARIABLES_FMT2,
)

def mark_main_statement(text: str, xmltool: XMLTool) -> str:
    cmdstmts = split_text_statements(text)
    index = infer_best_statement(cmdstmts, xmltool)
    cmdstmts = mark_main_stmt(index, cmdstmts)
    text = ''.join(cmdstmts)
    return text

def infer_best_statement(cmdstmts: list[str], xmltool: XMLTool) -> int:
    inferrer = MainStatementInferrer(cmdstmts, xmltool)
    return inferrer.infer()

def split_text_statements(text: str) -> list[str]:
    statements: list[str] = []
    
    current_stmt = ''
    lines = text.split('\n')
    for line in lines:
        delim_matches = expressions.get_matches(line, LINUX_STATEMENT_DELIMS)
        quoted_sections = expressions.get_quoted_sections(line)
        
        if not delim_matches:
            current_stmt += f'{line}\n'
            continue
        
        offset = 0
        for m in delim_matches:
            if quoted_sections[m.start()] == False and quoted_sections[m.end() - 1] == False:
                current_stmt += line[offset:m.end()]
                statements.append(current_stmt)
                current_stmt = ''
                offset = m.end()
                print()
        current_stmt += f'{line[offset:]}\n'

    statements.append(current_stmt)
    return statements

def mark_main_stmt(index: int, cmdstmts: list[str]) -> list[str]:
    cmdstmts[index] = f'\n__JANIS_MAIN__\n{cmdstmts[index]}\n__JANIS_MAIN__\n'
    return cmdstmts


# HELPER FUNCS: IDENTIFYING MAIN STATEMENT

def undefined_future_function(metric_records: list[StatementMetricRecord]) -> Optional[int]:
    # TODO statment with most flags? statement with redirect? statement firstword is not in list of known linux commands (except in cases like where the tool is actually 'awk')?
    raise NotImplementedError()

def perfect_requirement_match(metric_records: list[StatementMetricRecord]) -> Optional[int]:
    metric_records.sort(key=lambda x: x.reqsim, reverse=True)
    if metric_records[0].reqsim == 1 and metric_records[1].reqsim < 1:
        return metric_records[0].index
    return None

def top_all_metrics(metric_records: list[StatementMetricRecord]) -> Optional[int]:
    # same statement at the top of both metrics
    best_reqsim = sorted(metric_records, key=lambda x: x.reqsim, reverse=True)[0]
    best_gxrefs = sorted(metric_records, key=lambda x: x.gxrefs, reverse=True)[0]
    if best_reqsim.gxrefs == best_gxrefs.gxrefs:
        return best_reqsim.index
    return None

def requirement_priority(metric_records: list[StatementMetricRecord]) -> Optional[int]:
    # if 1+ statments firstword is very similar to the tool id (> 80% similar)
    # return the statement in this group with most gxrefs
    candidates = [x for x in metric_records if x.reqsim > 0.8]
    if candidates:
        return sorted(candidates, key=lambda x: x.gxrefs, reverse=True)[0].index
    return None

def gxref_priority(metric_records: list[StatementMetricRecord]) -> Optional[int]:
    # one statment has at least 3 galaxy references, and this is 2x more than others
    metric_records.sort(key=lambda x: x.gxrefs, reverse=True)
    if metric_records[0].gxrefs >= 3 and metric_records[0].gxrefs >= 2 * metric_records[1].gxrefs:
        return metric_records[0].index
    return None


# IDENTIFYING MAIN STATEMENT

@dataclass
class StatementMetricRecord:
    index: int
    statement: str
    gxrefs: int = 0 # the number of references to galaxy params in the statement
    reqsim: float = 0 # the similarity between the main tool requirement, and the first word of the statement


class MainStatementInferrer:
    
    def __init__(self, cmdstmts: list[str], xmltool: XMLTool) -> None:
        self.cmdstmts = cmdstmts
        self.xmltool = xmltool
        self.records: list[StatementMetricRecord] = []

    @property
    def requirement(self) -> str:
        if self.xmltool.metadata.main_requirement is not None:
            requirement = self.xmltool.metadata.main_requirement.name
        else:
            requirement = self.xmltool.metadata.id
        return requirement

    def infer(self) -> int:
        if len(self.cmdstmts) <= 1:
            return 0  # fallback: first statement 
        self.init_metric_records()
        self.set_gxref_counts()
        self.set_mainreq_similarities()
        return self.select_best_stmt()

    def init_metric_records(self) -> None:
        for i, statement in enumerate(self.cmdstmts):
            new_record = StatementMetricRecord(i, statement)
            self.records.append(new_record)
        
    def set_gxref_counts(self) -> None:
        for record in self.records:
            gxrefs = 0
            cmdstmt = record.statement
            var_matches = expressions.get_matches(cmdstmt, VARIABLES_FMT1)
            var_matches += expressions.get_matches(cmdstmt, VARIABLES_FMT2)
            for match in var_matches:
                varname = match.group(1)
                if self.xmltool.inputs.get(varname):
                    gxrefs += 1
            record.gxrefs = gxrefs
        
    def set_mainreq_similarities(self) -> None:
        banned_command_starters = ['cp', 'ln', 'mv']
        max_possible_score = global_align(self.requirement, self.requirement)
        for record in self.records:
            cmdstmt = record.statement
            words = cmdstmt.split()
            if len(words) == 0:
                record.reqsim = 0
            else:
                firstword = words[0]
                if firstword in banned_command_starters:
                    record.reqsim = -99999.9 # HACK
                else:
                    raw_similarity = global_align(firstword, self.requirement)
                    record.reqsim = raw_similarity / max_possible_score

    def select_best_stmt(self) -> int:
        identifiers = [
            top_all_metrics,
            gxref_priority,
            perfect_requirement_match,
            requirement_priority,
        ]
        for identifier_func in identifiers:
            selection = identifier_func(self.records)
            if selection is not None: 
                return selection
        return len(self.cmdstmts) - 1 # fallback: the last statement





