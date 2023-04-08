

from janis_core import Workflow, Tool
from .history import TaskInputCollector


class TaskInputsCategoriser:

    def __init__(self, tool: Tool, wf: Workflow) -> None:
        self.tool = tool
        self.wf = wf

        collector = TaskInputCollector(tool)
        collector.collect(wf)
        self.histories = collector.histories

        self.task_inputs: set[str] = set()
        self.param_inputs: set[str] = set()
        self.static_inputs: set[str] = set()
        self.ignored_inputs: set[str] = set()

    def categorise(self) -> None:
        for tinput_id, history in self.histories.items():
            # file types

            # RULE 1: files (mandatory) are always task inputs
            if history.is_file and not history.is_optional:  
                self.task_inputs.add(tinput_id)
            
            elif history.is_file:
                # RULE 2: files (optional) will be task inputs if 2+ values supplied
                if len(history.non_null_unique_values) >= 2:
                    self.task_inputs.add(tinput_id)

                # RULE 3: files (optional) with single invariant value should still be passed as task input
                elif len(history.unique_values) == 1 and len(history.non_null_unique_values) == 1:
                    self.task_inputs.add(tinput_id)
                
                # # RULE 3 (dep): files (optional) with single invariant value can be hard-coded param in task
                # elif len(history.unique_values) == 1 and len(history.non_null_unique_values) == 1:
                #     self.param_inputs.add(tinput_id)
                
                # RULE 4: files (optional) with single value (and some null values) must be task inputs
                elif len(history.unique_values) == 2 and len(history.non_null_unique_values) == 1:
                    self.task_inputs.add(tinput_id)

                # RULE 5: files (optional) will be ignored if unused
                elif len(history.non_null_unique_values) == 0:
                    self.ignored_inputs.add(tinput_id)
                
                else:
                    raise RuntimeError

            # non-file types
            else:
                # RULE 6: non-files will be task inputs if 2+ values supplied
                if len(history.non_null_unique_values) >= 2:
                    self.task_inputs.add(tinput_id)
                
                elif len(history.unique_values) == 1 and len(history.non_null_unique_values) == 1:
                    # RULE 7: non-files with single invariant InputNode value can be hardcoded param in task
                    if list(history.non_null_unique_values)[0].startswith('Input:'):
                        self.param_inputs.add(tinput_id)
                    
                    # RULE 8: non-files with single invariant static value can be hardcoded static in task
                    else:
                        self.static_inputs.add(tinput_id)

                # RULE 9: non-files with null and single InputNode value must be task inputs
                elif len(history.unique_values) == 2 and len(history.non_null_unique_values) == 1:
                    self.task_inputs.add(tinput_id)

                # RULE 10: non-files will be ignored if unused
                elif len(history.non_null_unique_values) == 0:
                    self.ignored_inputs.add(tinput_id)
                
                else:
                    raise RuntimeError



