

import unittest

from janis_core import WorkflowBuilder
from janis_core import WorkflowMetadata
from janis_core.workflow.workflow import InputNode
from janis_core.workflow.workflow import OutputNode
from janis_core.workflow.workflow import StepNode
from janis_core.workflow.workflow import ScatterMethod

from janis_core import CommandToolBuilder
from janis_core import CommandTool
from janis_core import ToolInput
from janis_core import ToolOutput
from janis_core import ToolMetadata
from janis_core import (
    InputSelector,
    WildcardSelector
)
from janis_core import (
    Stdout,
    Array,
    Float,
    Boolean,
    File,
    String
)

from galaxy2janis import datatypes

from janis_core.ingestion.janis_mapping.tool import to_janis_datatype
from janis_core.ingestion.janis_mapping.tool import to_janis_selector
from janis_core.ingestion.janis_mapping.tool import to_janis_metadata
from janis_core.ingestion.janis_mapping.tool import to_janis_tool_input
from janis_core.ingestion.janis_mapping.tool import to_janis_tool_output
from janis_core.ingestion.janis_mapping.tool import to_janis_tool

from janis_core.ingestion.janis_mapping.workflow import to_janis_workflow
from janis_core.ingestion.janis_mapping.workflow import to_janis_inputs_dict

from .mock.mock_tool import MOCK_TOOL_ABRICATE
from .mock.mock_workflow import MOCK_WORKFLOW


class TestJanisGeneralMapping(unittest.TestCase):
    """
    tests ability to map or generate janis datatypes / selectors etc
    from internal galaxy2janis tool objects.
    """
    def setUp(self) -> None:
        datatypes.populate()
        self.tool = MOCK_TOOL_ABRICATE

    def test_to_janis_datatype(self) -> None:
        """
        checks that galaxy2janis 'JanisDatatype' objects can be mapped to
        janis 'DataType' objects 
        """
        jtype1 = to_janis_datatype(self.tool.inputs[0])  
        jtype2 = to_janis_datatype(self.tool.inputs[1])
        jtype3 = to_janis_datatype(self.tool.inputs[2])
        jtype4 = to_janis_datatype(self.tool.outputs[0])
        # check janis objects are correct
        self.assertIsInstance(jtype1, File)
        self.assertIsInstance(jtype2, Boolean)
        self.assertIsInstance(jtype3, Array)
        self.assertIsInstance(jtype4, Stdout)
        # check attributes are correct
        self.assertEquals(jtype2.optional, True)
        self.assertIsInstance(jtype3.subtype(), Float)
        self.assertIsInstance(jtype4.subtype, File)

    def test_to_janis_selector(self) -> None:
        """
        checks that galaxy2janis 'OutputComponent' objects can generate 
        janis 'Selector' objects
        """
        jsel1 = to_janis_selector(self.tool.outputs[0])
        jsel2 = to_janis_selector(self.tool.outputs[1])
        jsel3 = to_janis_selector(self.tool.outputs[2])
        # check janis objects are correct
        self.assertIsInstance(jsel2, WildcardSelector)
        self.assertIsInstance(jsel3, InputSelector)
        self.assertIsNone(jsel1)
        # check attributes are correct
        self.assertEquals(jsel2.wildcard, 'report.txt')
        self.assertEquals(jsel3.input_to_select, 'fileInput')


class TestJanisToolMapping(unittest.TestCase):
    """
    tests ability to map or generate janis tool objects 
    from internal galaxy2janis tool objects.
    MOCK_TOOL is used to test mapping functions.
    """
    def setUp(self) -> None:
        datatypes.populate()
        self.tool = MOCK_TOOL_ABRICATE
    
    def test_to_janis_tool(self) -> None:
        """
        checks that galaxy2janis 'Tool' objects can be mapped to
        janis 'CommandToolBuilder' objects 
        """
        jtool = to_janis_tool(self.tool)
        self.assertIsInstance(jtool, CommandToolBuilder)

    def test_to_janis_tool_input(self) -> None:
        """
        checks that galaxy2janis 'InputComponent' objects can be mapped to
        janis 'ToolInput' objects 
        """
        jinp1 = to_janis_tool_input(self.tool.inputs[0])  
        jinp2 = to_janis_tool_input(self.tool.inputs[1])
        jinp3 = to_janis_tool_input(self.tool.inputs[2])
        jinp4 = to_janis_tool_input(self.tool.inputs[3])
        # check janis objects are correct
        self.assertIsInstance(jinp1, ToolInput)
        self.assertIsInstance(jinp2, ToolInput)
        self.assertIsInstance(jinp3, ToolInput)
        self.assertIsInstance(jinp4, ToolInput)
        # check attributes are correct
        self.assertEquals(jinp1.tag, 'fileInput')
        self.assertEquals(jinp1.prefix, None)
        self.assertEquals(jinp1.separate_value_from_prefix, True)
        self.assertIsInstance(jinp1.input_type, File)
        
        self.assertEquals(jinp2.tag, 'noheader')
        self.assertEquals(jinp2.prefix, '--noheader')
        self.assertIsInstance(jinp2.input_type, Boolean)
        self.assertEquals(jinp2.input_type.optional, True)

        self.assertEquals(jinp3.tag, 'minid')
        self.assertEquals(jinp3.prefix, '--minid=')
        self.assertEquals(jinp3.separate_value_from_prefix, False)
        self.assertIsInstance(jinp3.input_type, Array)
        
        self.assertEquals(jinp4.tag, 'db')
        self.assertEquals(jinp4.prefix, '--db=')
        self.assertEquals(jinp4.separate_value_from_prefix, False)
        self.assertEquals(jinp4.default, 'resfinder')
        self.assertIsInstance(jinp4.input_type, String)

    def test_to_janis_tool_output(self) -> None:
        """
        checks that galaxy2janis 'OutputComponent' objects can be mapped to
        janis 'ToolOutput' objects 
        """
        jout1 = to_janis_tool_output(self.tool.outputs[0])
        jout2 = to_janis_tool_output(self.tool.outputs[1])
        jout3 = to_janis_tool_output(self.tool.outputs[2])
        # check janis objects are correct
        self.assertIsInstance(jout1, ToolOutput)
        self.assertIsInstance(jout2, ToolOutput)
        self.assertIsInstance(jout3, ToolOutput)
        # check attributes are correct
        self.assertIsInstance(jout1.output_type, Stdout)
        self.assertEquals(jout1.tag, 'outReport1')
        self.assertIsNone(jout1.selector)
        
        self.assertIsInstance(jout2.output_type, File)
        self.assertEquals(jout2.tag, 'outReport2')
        self.assertIsInstance(jout2.selector, WildcardSelector)
        self.assertEquals(jout2.selector.wildcard, 'report.txt')
        
        self.assertIsInstance(jout3.output_type, File)
        self.assertEquals(jout3.tag, 'outFileInput')
        self.assertIsInstance(jout3.selector, InputSelector)
        self.assertEquals(jout3.selector.input_to_select, 'fileInput')

    def test_to_janis_metadata(self) -> None:
        """
        checks that galaxy2janis 'ToolMetadata' objects can be mapped to
        janis 'ToolXMLMetadata' objects 
        """
        # MOCK_TOOL
        jmeta = to_janis_metadata(self.tool.metadata)
        self.assertIsInstance(jmeta, ToolMetadata)
        self.assertIsNotNone(jmeta.contributors)
        self.assertIsNotNone(jmeta.dateCreated)
        self.assertIsNotNone(jmeta.dateUpdated)
        self.assertIsNotNone(jmeta.citation)
        self.assertIsNotNone(jmeta.documentation)
        self.assertIsNotNone(jmeta.short_documentation)
        self.assertIsNotNone(jmeta.version)


class TestJanisWorkflowMapping(unittest.TestCase):
    """
    tests ability to map or generate janis workflow objects 
    from internal galaxy2janis workflow objects
    """
    def setUp(self) -> None:
        datatypes.populate()
        self.workflow = MOCK_WORKFLOW
        self.jworkflow = to_janis_workflow(self.workflow)
        self.jinputs = to_janis_inputs_dict(self.workflow)

    def test_to_janis_workflow(self) -> None:
        self.assertIsInstance(self.jworkflow, WorkflowBuilder)
    
    def test_to_janis_inputs_dict(self) -> None:
        # single input, no value
        self.assertEquals(self.jinputs['inFasta'], None)
        # single input, provided value
        self.workflow.inputs[0].value = 'path/to/file.fasta'
        jinputs = to_janis_inputs_dict(self.workflow)
        self.assertEquals(jinputs['inFasta'], 'path/to/file.fasta')

    def test_janis_metadata(self) -> None:
        self.assertIsInstance(self.jworkflow.metadata, WorkflowMetadata)

    def test_janis_inputs(self) -> None:
        target_inputs = {
            'inFasta', 
        }
        actual_inputs = set(self.jworkflow.input_nodes.keys())
        self.assertEquals(target_inputs, actual_inputs)

        jinp = self.jworkflow.input_nodes['inFasta']
        self.assertIsInstance(jinp, InputNode)
        self.assertIsInstance(jinp.datatype, File)

    def test_janis_outputs(self) -> None:
        target_outputs = {'abricate_outReport1'}
        actual_outputs = set(self.jworkflow.output_nodes.keys())
        self.assertEquals(target_outputs, actual_outputs)

        jout = self.jworkflow.output_nodes['abricate_outReport1']
        self.assertIsInstance(jout, OutputNode)
        self.assertIsInstance(jout.datatype, Stdout)
        self.assertIsInstance(jout.datatype.subtype, File)
        self.assertEquals(jout.doc.doc, 'report file')

    def test_janis_steps(self) -> None:
        target_steps = {'abricate'}
        actual_steps = set(self.jworkflow.step_nodes.keys())
        self.assertEquals(target_steps, actual_steps)

        # basic object checks
        jstep = self.jworkflow.step_nodes['abricate']
        self.assertIsInstance(jstep, StepNode)
        self.assertIsInstance(jstep.tool, CommandTool)

        # sources (tool input values)
        self.assertIn('fileInput', jstep.sources)
        self.assertNotIn('noheader', jstep.sources)
        self.assertNotIn('minid', jstep.sources)
        self.assertNotIn('db', jstep.sources)
        
        # scatter 
        # self.assertEquals(jstep.scatter.fields, ['minid'])
        # self.assertEquals(jstep.scatter.method, ScatterMethod.dot)


