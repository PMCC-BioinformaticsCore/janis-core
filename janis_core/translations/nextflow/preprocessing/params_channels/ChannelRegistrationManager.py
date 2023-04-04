


from janis_core.workflow.workflow import InputNode
from janis_core.types import File, Filename, Directory, DataType
from janis_core import translation_utils as utils
from janis_core import settings

from ... import params
from ... import channels
from ...scope import Scope


class ChannelRegistrationManager:
    def __init__(self, inp: InputNode, scope: Scope) -> None:
        self.inp = inp
        self.scope = scope

    @property
    def basetype(self) -> DataType:
        basetype = utils.get_base_type(self.inp.datatype)
        basetype = utils.ensure_single_type(basetype)
        return basetype

    @property
    def is_subworkflow(self) -> bool:
        if self.scope.labels != [settings.translate.nextflow.NF_MAIN_NAME]:
            return True
        return False
    
    @property
    def method(self) -> str:
        # if isinstance(self.basetype, (File, Filename, Directory)) and not self.inp.datatype.optional:
        if isinstance(self.basetype, (File, Filename, Directory)):
            return 'fromPath'
        return 'of'

    @property
    def source(self) -> str:
        if self.is_subworkflow:
            src = ''
        else:
            if not params.exists(self.inp.uuid):
                print()
            param_name = params.get(self.inp.uuid).name
            # @secondaryarrays
            if utils.is_array_secondary_type(self.inp.datatype):
                src = f'params.{param_name}.flatten()'
            else:
                src = f'params.{param_name}'
        return src

    def register(self) -> None:
        operations = self.get_operations()
        channels.add(
            scope=self.scope,
            janis_tag=self.inp.id(),
            method=self.method,
            source=self.source,
            operations=operations,
            janis_dtype=self.inp.datatype,
            janis_uuid=self.inp.uuid,
            define=False if self.is_subworkflow else True
        )

    def get_operations(self) -> str:
        # secondary array
        if utils.is_array_secondary_type(self.inp.datatype):
            ops = self.get_operations_secondary_array()
        
        # secondary
        elif utils.is_secondary_type(self.inp.datatype):
            ops = self.get_operations_secondary()
        
        # file array
        elif isinstance(self.basetype, (File, Filename, Directory)) and self.inp.datatype.is_array():
            ops = self.get_operations_file_array()
        
        # nonfile array
        elif self.inp.datatype.is_array():
            ops = self.get_operations_nonfile_array()
        
        # anything else
        else:
            ops = self.get_operations_generic()
        return ops

    def get_operations_secondary_array(self) -> str:
        exts = utils.get_extensions(self.basetype)
        size = len(exts)
        
        ops: str = ''
        ops += f'.collate( {size} )'
        # if self.inp.datatype.optional:
        #     ops += '.ifEmpty( null )'
        return ops

    def get_operations_secondary(self) -> str:
        ops: str = ''
        ops += '.toList()'
        # if self.inp.datatype.optional:
        #     ops += '.ifEmpty( null )'
        return ops

    def get_operations_file_array(self) -> str:
        ops: str = ''
        ops += '.toList()'
        # if self.inp.datatype.optional:
        #     ops += '.ifEmpty( null )'
        return ops

    def get_operations_nonfile_array(self) -> str:
        ops: str = ''
        ops += '.toList()'
        # if self.inp.datatype.optional:
        #     ops += '.ifEmpty( null )'
        return ops

    def get_operations_generic(self) -> str:
        ops: str = ''
        # if self.inp.datatype.optional:
        #     ops += '.ifEmpty( null )'
        return ops
    
