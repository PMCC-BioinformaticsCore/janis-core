
import os 
import unittest
from janis_core import CommandToolBuilder
from janis_core.ingestion.wdl import WdlParser
from janis_core.ingestion.wdl.parsing import parse_task
from janis_core.ingestion import ingest
from janis_core import (
    File, 
    String, 
    Int,
    Boolean, 
    Array,
)
import WDL

WDL_TESTDATA_PATH = os.path.join(os.getcwd(), 'janis_core/tests/data/wdl')


class TestCommandParsing(unittest.TestCase):
    
    def test_rename_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/rename_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        tool = parse_task(task)
        self.assertIsInstance(tool, CommandToolBuilder)
        self.assertEqual(tool.base_command, ['cp'])
        self.assertEqual(len(tool._inputs), 2)
        self.assertEqual(tool._inputs[0].id(), 'sourceFile')
        self.assertIsInstance(tool._inputs[0].position, File)
        self.assertEqual(tool._inputs[0].position, 1)
        self.assertEqual(tool._inputs[0].prefix, None)
        self.assertEqual(tool._inputs[1].id(), 'targetFilename')
        self.assertIsInstance(tool._inputs[1].position, String)
        self.assertEqual(tool._inputs[1].position, 2)
        self.assertEqual(tool._inputs[1].prefix, None)
    
    def test_io_tool(self) -> None:
        filepath = f'{WDL_TESTDATA_PATH}/basic/io_tool.wdl'
        d = WDL.load(filepath)
        task = d.tasks[0]
        tool = parse_task(task)
        
        ### basics ###
        self.assertIsInstance(tool, CommandToolBuilder)
        self.assertEqual(tool.base_command, ['echo'])
        self.assertEqual(len(tool._inputs), 7)
        self.assertEqual(len(tool._outputs), 5)
        
        ### inputs ###
        # inInt
        tinp = tool._inputs[0]
        self.assertEqual(tinp.id(), 'inInt')
        self.assertIsInstance(tinp.input_type, Int)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 1)
        self.assertIsNone(tinp.prefix)
        # inStr
        tinp = tool._inputs[1]
        self.assertEqual(tinp.id(), 'inStr')
        self.assertIsInstance(tinp.input_type, String)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 2)
        self.assertIsNone(tinp.prefix)
        # inBool
        tinp = tool._inputs[2]
        self.assertEqual(tinp.id(), 'inBool')
        self.assertIsInstance(tinp.input_type, Boolean)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 3)
        self.assertEqual(tinp.prefix, '--flag1')
        # inFile
        tinp = tool._inputs[3]
        self.assertEqual(tinp.id(), 'inFile')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 4)
        self.assertEqual(tinp.prefix, '--in-file=')
        self.assertFalse(tinp.separate_value_from_prefix)
        # inFileOpt
        tinp = tool._inputs[4]
        self.assertEqual(tinp.id(), 'inFileOpt')
        self.assertIsInstance(tinp.input_type, File)
        self.assertTrue(tinp.input_type.optional)
        self.assertEqual(tinp.position, 5)
        self.assertIsNone(tinp.prefix)
        # inFileArr
        tinp = tool._inputs[5]
        self.assertEqual(tinp.id(), 'inFileArr')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 6)
        self.assertIsNone(tinp.prefix)
        # inSecondary
        tinp = tool._inputs[6]
        self.assertEqual(tinp.id(), 'inSecondary')
        self.assertIsInstance(tinp.input_type, File)
        self.assertNotEqual(tinp.input_type.optional, True)
        self.assertEqual(tinp.position, 7)
        self.assertIsNone(tinp.prefix)

        ### outputs ###
        tout = tool._outputs[0]

        

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