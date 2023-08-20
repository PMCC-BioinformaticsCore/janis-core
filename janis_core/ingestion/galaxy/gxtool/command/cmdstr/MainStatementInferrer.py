

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional

from janis_core.ingestion.galaxy.utils.general import global_align
from .DynamicCommandStatement import DynamicCommandStatement
from .CommandString import CommandStringSource

# functions to identify best statement
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



@dataclass
class StatementMetricRecord:
    index: int
    statement: DynamicCommandStatement
    gxrefs: int = 0 # the number of references to galaxy params in the statement
    reqsim: float = 0 # the similarity between the main tool requirement, and the first word of the statement


class MainStatementInferrer:
    def __init__(self, statements: list[DynamicCommandStatement], source: CommandStringSource, requirement: str):
        self.statements = statements
        self.source = source
        self.requirement = requirement
        self.metric_records: list[StatementMetricRecord] = []

    def infer(self) -> int:
        if len(self.statements) <= 1:
            return 0  # fallback: first statement 

        self.init_metric_record()
        self.set_gxref_counts()
        self.set_mainreq_similarities()
        if self.source == CommandStringSource.XML:
            best = self.choose_best_statement_vanilla()
        elif self.source in [CommandStringSource.TEST, CommandStringSource.TOOL_STATE]:
            best = self.choose_best_statement_templated()
        else:
            raise RuntimeError(f'Unknown source: {self.source}')
        return best
        
    def init_metric_record(self) -> None:
        for i, statement in enumerate(self.statements):
            new_record = StatementMetricRecord(i, statement)
            self.metric_records.append(new_record)
        
    def set_gxref_counts(self) -> None:
        for record in self.metric_records:
            record.gxrefs = record.statement.get_galaxy_reference_count()
        
    def set_mainreq_similarities(self) -> None:
        banned_command_starters = ['cp', 'ln', 'mv']
        max_possible_score = global_align(self.requirement, self.requirement)
        for record in self.metric_records:
            if len(record.statement.realised_tokens) == 0:
                record.reqsim = 0
            elif record.statement.get_first_word() in banned_command_starters:
                record.reqsim = -99999.9 # HACK
            else:
                raw_similarity = global_align(record.statement.get_first_word(), self.requirement)
                record.reqsim = raw_similarity / max_possible_score
    
    def choose_best_statement_vanilla(self) -> int:
        identifiers = [
            top_all_metrics,
            gxref_priority,
            perfect_requirement_match,
            requirement_priority,
        ]
        for identifier_func in identifiers:
            selection = identifier_func(self.metric_records)
            if selection is not None: 
                return selection
        return len(self.statements) - 1 # fallback: the last statement

    def choose_best_statement_templated(self) -> int:
        identifiers = [
            perfect_requirement_match,
            top_all_metrics,
            requirement_priority,
            #gxref_priority
        ]
        for identifier_func in identifiers:
            selection = identifier_func(self.metric_records)
            if selection is not None: 
                return selection
        return len(self.statements) - 1 # fallback: the last statement

      