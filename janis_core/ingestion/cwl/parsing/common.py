

from typing import Any
from abc import ABC, abstractmethod
from janis_core.messages import log_warning
from ..expressions import parse_basic_expression


class EntityParser(ABC):
    def __init__(self, cwl_utils: Any) -> None:
        self.cwl_utils = cwl_utils
        self.error_msgs: list[str] = []
        self.success: bool = False

    @abstractmethod
    def parse(self, entity: Any) -> Any:
        ...

    def log_messages(self, entity_uuid: str) -> None:
        for msg in self.error_msgs:
            log_warning(entity_uuid, msg)



class RequirementsParser(EntityParser):

    def parse(self, entity: Any, janis_uuid: str) -> dict[str, Any]:  # type: ignore I FKN KNOW

        out: dict[str, Any] = {
            'container': 'ubuntu:latest',
            'memory': None,
            'cpus': None,
            'time': None,
            'files_to_create': {},
            'env_vars': {}
        }

        for req in entity.requirements or []:

            if isinstance(req, self.cwl_utils.DockerRequirement):
                out['container'] = str(req.dockerPull)

            elif isinstance(req, self.cwl_utils.EnvVarRequirement):
                for envdef in req.envDef:
                    # entry name
                    name_expr, success = parse_basic_expression(envdef.envName)
                    if not success:
                        msg = 'untranslated javascript expression in environment variable name'
                        self.error_msgs.append(msg)
                    # entry 
                    entry_expr, success = parse_basic_expression(envdef.envValue)
                    if not success:
                        msg = 'untranslated javascript expression in environment variable value'
                        self.error_msgs.append(msg)
                    out['env_vars'][name_expr] = entry_expr

            elif isinstance(req, self.cwl_utils.InitialWorkDirRequirement):
                anonymous_files_count: int = 0
                for dirent in req.listing:
                    if isinstance(dirent, str):
                        anonymous_files_count += 1
                        label = f'unnamed_{anonymous_files_count}'
                        # entry
                        entry_expr, success = parse_basic_expression(dirent)
                        if not success:
                            msg = 'untranslated javascript expression in file / directory localisation value'
                            self.error_msgs.append(msg)
                        out['files_to_create'][label] = entry_expr

                    else:
                        # entry name
                        name_expr, success = parse_basic_expression(dirent.entryname)
                        if not success:
                            msg = 'untranslated javascript expression in file / directory localisation name'
                            self.error_msgs.append(msg)
                        # entry 
                        entry_expr, success = parse_basic_expression(dirent.entry)
                        if not success:
                            msg = 'untranslated javascript expression in file / directory localisation value'
                            self.error_msgs.append(msg)
                        out['files_to_create'][name_expr] = entry_expr

            elif isinstance(req, self.cwl_utils.ResourceRequirement):
                # maybe convert mebibytes to megabytes?
                memory, success = parse_basic_expression(req.ramMin or req.ramMax)
                out['memory'] = memory
                if not success:
                    msg = 'untranslated javascript expression in task MEM requirement'
                    self.error_msgs.append(msg)

                cpus, success = parse_basic_expression(req.coresMin)
                out['cpus'] = cpus
                if not success:
                    msg = 'untranslated javascript expression in task CPU requirement'
                    self.error_msgs.append(msg)

            elif hasattr(req, 'timelimit') and isinstance(req, self.cwl_utils.ToolTimeLimit):
                time, success = parse_basic_expression(req.timelimit)
                out['time'] = time
                if not success:
                    msg = 'untranslated javascript expression in task TIME requirement'
                    self.error_msgs.append(msg)
        
        self.log_messages(janis_uuid)
        return out