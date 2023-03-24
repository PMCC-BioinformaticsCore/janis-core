

from typing import Optional
from collections import defaultdict

from janis_core import CommandTool, PythonTool
from . import data_sources
from ..scope import Scope


class VariableManager:
    """
    within a process, for a given TInput, get the variable which currently represents that TInput.
    eg 

    input:
    1 path indexed_bam_flat
    2 path cath_crossmapped, stageAs: 'cath_crossmapped'

    script:
    3 def in_bams = get_primary_files(indexed_bam_flat, 2)
    4 def in_bams_joined = in_bams.join(' ')
    5 def cath_crossmapped = cath_crossmapped ? "-cx ${cath_crossmapped}" : ""
    '''
    6 echo \
    7 ${in_bams_joined} \
    8 ${cath_crossmapped} \
    9 --promoter=${params.promoter_bp} \
    '''

    TInput('inBams')
    line 1 -> indexed_bam_flat
    line 3 -> in_bams
    line 4 -> in_bams_joined
    line 6 -> in_bams_joined

    TInput('cathCrossmapped')
    line 2 -> cath_crossmapped
    line 8 -> cath_crossmapped

    TInput('promoterBp')
    line 9 -> params.promoter_bp

    data_structure description:
    {
        tinput_id: [
            first,  (original)
            second, 
            third,  (current)
            ..
        ]
    }

    """
    def __init__(self, scope: Scope) -> None:
        self.scope = scope
        self.data_structure: dict[str, list[Optional[str | list[str]]]] = defaultdict(list)

    def update(self, tinput_id: str, varname: Optional[str | list[str]]) -> None:
        self.data_structure[tinput_id].append(varname)

    def update_for_tool(self, tool: CommandTool | PythonTool) -> None:
        for tinput in tool.inputs():
            src = data_sources.get_variable(self.scope, tinput)
            self.update(tinput.id(), src)

    def original(self, tinput_id: str, index: Optional[int]=None) -> Optional[str]:
        """gets the original varname for this tinput_id in current scope"""
        varname = self.data_structure[tinput_id][0]
        varname = self._apply_index(varname, index)
        return varname
    
    def current(self, tinput_id: str, index: Optional[int]=None) -> Optional[str]:
        """gets the current varname for this tinput_id in current scope"""
        if not self.data_structure[tinput_id]:
            print()
        varname = self.data_structure[tinput_id][-1]
        varname = self._apply_index(varname, index)
        return varname
    
    def history(self, tinput_id: str, index: Optional[int]=None) -> list[Optional[str]]:
        """gets each varname this tinput_id has had in the current scope"""
        out: list[Optional[str]] = []
        for varname in self.data_structure[tinput_id]:
            out.append(self._apply_index(varname, index))
        return out
    
    # helper methods
    def _apply_index(self, varname: Optional[str | list[str]], index: Optional[int]) -> Optional[str]:
        """
        for same variables, the varname is actually a list of varnames.
        
        this will happen in the case of a ToolInput, where the type is a secondary file eg BamBai. 
        in the nextflow process, this bambai may appear like the following:
            inputs:
            tuple path(bam), path(bai)
        
        in this case, the ToolInput is split into the varname ['bam', 'bai']
        to get the current varname for the bam, we want varname[0].
        to get the current varname for the bai, we want varname[1].

        """
        if isinstance(varname, list):
            if index is not None:
                varname = varname[index]
            else:
                varname = varname[0]
        return varname

