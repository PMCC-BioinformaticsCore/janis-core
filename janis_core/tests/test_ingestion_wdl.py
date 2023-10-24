
import os 
import unittest
from janis_core import CommandToolBuilder
from janis_core.ingestion.wdl import WdlParser
from janis_core.ingestion.wdl.parsing import parse_task
from janis_core.ingestion import ingest
import WDL

WDL_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/wdl')


class TestCommandParsing(unittest.TestCase):
    
    def test_rename_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/rename_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        tool = parse_task(task)
        assert isinstance(tool, CommandToolBuilder)
        assert tool.base_command == ['cp']
        assert len(tool._inputs) == 2
        assert tool._inputs[0].id() == 'sourceFile'
        assert tool._inputs[0].position == 1
        assert tool._inputs[0].prefix == None
        assert tool._inputs[1].id() == 'targetFilename'
        assert tool._inputs[1].position == 2
        assert tool._inputs[1].prefix == None
    
    def test_io_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/io_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)
    
    def test_fastqc_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/fastqc.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)
    
    def test_bwa_mem(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/bwa_mem.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        parse_task(task)



class TestFromWdl(unittest.TestCase):
    parser = WdlParser()
    
    def test_ingest_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/rename_tool.wdl'
        tool = ingest(filepath, 'wdl')
        raise NotImplementedError

    def test_ingest_workflow(self) -> None:
        raise NotImplementedError