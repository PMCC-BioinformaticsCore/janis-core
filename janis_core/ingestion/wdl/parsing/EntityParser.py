

from abc import abstractmethod
from typing import Any
import WDL 

from janis_core import settings
from janis_core import CommandToolBuilder, WorkflowBuilder


class EntityParser:

    def parse(self) -> Any:
        # normal mode
        if settings.ingest.SAFE_MODE:
            try:
                j_entity = self.do_parse()
                self.success = True
            except Exception:
                j_entity = self.fallback()
        
        # dev mode
        else:
            j_entity = self.do_parse()
            self.success = True
        
        return j_entity

    @abstractmethod
    def do_parse(self) -> Any:
        ...
    
    @abstractmethod
    def fallback(self) -> Any:
        ...



class TaskParser(EntityParser):

    def __init__(self, task: WDL.Tree.Task, cmdtool: CommandToolBuilder) -> None:
        self.task = task 
        self.cmdtool = cmdtool 
        self.success = False


class WorkflowParser(EntityParser):

    def __init__(self, wdl_wf: WDL.Tree.Workflow, janis_wf: WorkflowBuilder) -> None:
        self.wdl_wf = wdl_wf 
        self.janis_wf = janis_wf 
        self.success = False