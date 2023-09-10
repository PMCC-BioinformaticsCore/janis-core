
import os
import tempfile
import regex as re
from typing import Any

from galaxy.tools import Tool as GxTool
from galaxy.tools import create_tool_from_source
from galaxy.tools.parameters.basic import ToolParameter
from galaxy.tool_util.parser import get_tool_source
from galaxy.tool_util.parser.output_objects import ToolOutput
from galaxy.model import History

from ..model import XMLDataParam
from ..model import XMLConfigfile
from ..model import XMLScript
from ..model import XMLTool
from ..model import XMLMetadata
from ..model import XMLCitation
from ..model import XMLParamRegister
from ..model import XMLCondaRequirement
from ..model import XMLContainerRequirement
from ..model import XMLTest
from ..model import XMLTestRegister

from ...expressions.patterns import GX_TOOL_SCRIPT
from ...expressions.matches import get_matches
from ..mock import MockApp, MockObjectStore

from .param_flattener import XMLParamFlattener
from .outputs import parse_output_param
from .inputs import parse_input_param

Requirement = XMLContainerRequirement | XMLCondaRequirement


def load_xmltool(path: str) -> XMLTool:
    gxtool = _load_galaxy_tool(path)
    factory = GalaxyToolFactory(gxtool, path)
    return factory.create()

def _load_galaxy_tool(path: str) -> GxTool:
    app = _get_app()
    tool_source = get_tool_source(path)
    tool = create_tool_from_source(app, tool_source)
    tool.assert_finalized()
    return tool



def _get_app() -> MockApp:
    # basic details
    app = MockApp()
    app.job_search = None
    app.object_store = MockObjectStore()
    # config
    app.config.new_file_path = os.path.join(tempfile.mkdtemp(), "new_files")
    app.config.admin_users = "grace@thebest.com"
    app.config.len_file_path = "moocow"
    # database
    
    app.model.context.add(History())
    app.model.context.flush()
    return app

class GalaxyToolFactory:
    def __init__(self, gxtool: GxTool, xmlpath: str):
        self.gxtool = gxtool
        self.xmlpath = xmlpath
        self.scripts: list[XMLScript] = []
        self.configfiles: list[XMLConfigfile] = []
        self.inputs = self.parse_inputs()
        self.parse_scripts()
        self.parse_configfiles()
        self.outputs = self.parse_outputs()

    # TODO remove script references in command parsing etc
    # scripts / configfiles treated as gxparam XMLInputParams. 
    # still need way to set default value? 

    def create(self) -> XMLTool:
        return XMLTool(
            metadata=self.parse_metadata(),
            raw_command=self.parse_command(),
            configfiles=self.configfiles,
            scripts=self.scripts,
            inputs=self.inputs,
            outputs=self.outputs,
            tests=self.parse_tests()
        ) 
    
    def parse_inputs(self) -> XMLParamRegister:
        """returns a an InputRegister by reformatting the galaxy tool representation's params."""
        register = XMLParamRegister()
        g_in_params = self.flatten_params()

        # removing duplicates
        fingerprints = set()
        unique_params = []
        for g_param in g_in_params:
            items = [g_param.flat_name, g_param.label, g_param.type]
            items = [x for x in items if isinstance(x, str)]
            fingerprint = ''.join(items)
            if fingerprint not in fingerprints:
                unique_params.append(g_param)
                fingerprints.add(fingerprint)

        # parsing individual params
        for g_param in unique_params:
            i_param = parse_input_param(g_param)
            register.add(i_param)
            
        return register

    def flatten_params(self) -> list[ToolParameter]:
        pf = XMLParamFlattener(self.gxtool.inputs)
        return pf.flatten()

    def parse_outputs(self) -> XMLParamRegister:
        """returns a formatted list of outputs using the representation"""
        register = XMLParamRegister()
        g_out_params: list[ToolOutput] = list(self.gxtool.outputs.values())
        for g_param in g_out_params:
            t_params = parse_output_param(g_param, self.inputs)
            for param in t_params:
                register.add(param)
        return register

    def parse_metadata(self) -> XMLMetadata:
        """returns a formatted Metadata using the representation"""
        requirements: list[Requirement] = self.get_requirements()
        citations: list[XMLCitation] = self.get_citations()
        return XMLMetadata(
            name=str(self.gxtool.name),  #type: ignore
            id=str(self.gxtool.id),  #type: ignore
            version=str(self.gxtool.version).split('+galaxy')[0],  #type: ignore
            description=str(self.gxtool.description),  #type: ignore
            help=str(self.gxtool.raw_help),  #type: ignore
            requirements=requirements,
            citations=citations,
            creator=self.gxtool.creator  #type: ignore
        )
    
    def get_requirements(self) -> list[Requirement]:
        """returns a formatted list of Requirements using the representation"""
        reqs: list[Requirement] = []
        reqs += self.get_conda_requirements()
        reqs += self.get_container_requirements()
        return reqs
    
    def get_conda_requirements(self) -> list[XMLCondaRequirement]:
        packages: list[dict[str, str]] = self.gxtool.requirements.to_list() # type: ignore
        return [XMLCondaRequirement(p['name'], p['version']) for p in packages]
    
    def get_container_requirements(self) -> list[XMLContainerRequirement]:
        containers: list[ContainerDescription] = self.gxtool.containers # type: ignore
        return [XMLContainerRequirement(c.identifier) for c in containers] # type: ignore

    def get_citations(self) -> list[XMLCitation]:
        citations: list[XMLCitation] = []
        citations += self.get_biotools_citations()
        citations += self.get_doi_citations()
        #citations += self.get_bibtex_citations()
        return citations

    def get_biotools_citations(self) -> list[XMLCitation]: 
        out: list[XMLCitation] = []
        for ref in self.gxtool.xrefs:
            citation = XMLCitation(
                type='biotools',
                text=f"https://bio.tools/{ref['value']}"
            )
            out.append(citation)
        return out

    def get_doi_citations(self) -> list[XMLCitation]:
        out: list[XMLCitation] = []
        for elem in self.gxtool.tool_source.xml_tree.findall('citations'):
            for ref in elem.findall('citation'):
                if ref.attrib['type'] == 'doi':
                    citation = XMLCitation(
                        type='doi',
                        text=f"https://doi.org/{ref.text}"
                    )
                    out.append(citation)
        return out
        
    def get_bibtex_citations(self) -> list[XMLCitation]:
        out: list[XMLCitation] = []
        for elem in self.gxtool.tool_source.xml_tree.findall('citations'):
            for ref in elem.findall('citation'):
                if ref.attrib['type'] == 'bibtex':
                    citation = XMLCitation(
                        type='bibtex',
                        text=self.parse_bibtex(ref.text)
                    )
                    out.append(citation)
        return out
    
    def parse_bibtex(self, bibtex_citation: dict[str, str]) -> str:
        # define and parse using biblib
        bp = bib.Parser()
        data = bp.parse(bibtex_citation).get_entries() # type: ignore
        # get each citation key: value pair
        entry = list(data.values())[0]  # type: ignore
        if 'url' in entry:
            return f"{entry['url']}" # type: ignore
        elif 'author' in entry and 'title' in entry:
            return f"{entry['author']}.  {entry['title']}"
        return ''

    def parse_command(self) -> str:
        """returns the tool xml command"""
        return str(self.gxtool.command) # type: ignore
    
    def parse_configfiles(self) -> None:
        """
        parses tool configfiles & adds to XMLTool.configfiles. 
        Additionally:
            Adds a XMLDataParam for the configfile. 
            Do not need to replace the configfile in the command, as is already a reference. 
            This way, configfile will be treated as ToolInputs which is desired behaviour. 
            
            Need to add a Workflow InputNode for the configfile. 
            Need to link this InputNode to the extracted ToolInput in the step call. 
            These happen later, in ingestion.galaxy.gx.gxworkflow.values.scripts 
            (called from janis_core.ingestion.galaxy.ingest)
        """
        for name, _, contents in self.gxtool.config_files:  # type: ignore
            if isinstance(contents, str):
                # parse into XMLConfigfile & add to XMLTool
                configfile = XMLConfigfile(name, contents)  # type: ignore
                self.configfiles.append(configfile)

                # add param for XMLConfigfile
                param = XMLDataParam(name=configfile.varname)
                param.formats = ['file']
                param.helptext = 'galaxy script needed to run tool'
                self.inputs.add(param)

    def parse_scripts(self) -> None:
        """
        parses local tool scripts & adds to XMLTool.scripts.  
        (ie $__tool_directory__/my_script.py)
        Additionally:
            Adds a XMLDataParam for the script. 
            Replaces the script in the command with a reference to the param.
            This way, scripts will be treated as ToolInputs which is desired behaviour. 
            
            Need to add a Workflow InputNode for the script file. 
            Need to link this InputNode to the extracted ToolInput in the step call. 
            These happen later, in ingestion.galaxy.gx.gxworkflow.values.scripts 
        (called from janis_core.ingestion.galaxy.ingest)
        """
        while get_matches(self.gxtool.command, GX_TOOL_SCRIPT):
            # get the match
            match = get_matches(self.gxtool.command, GX_TOOL_SCRIPT)[0]
            
            # parse match into script & add to XMLTool
            script = self.parse_script(match)
            self.scripts.append(script)

            # add param for script
            param = XMLDataParam(name=script.varname)
            param.formats = ['file']
            param.helptext = 'galaxy script needed to run tool'
            self.inputs.add(param)
            
            # modify command to replace script with param
            print(self.gxtool.command)
            self.replace_script_in_command(match, param)
            print(self.gxtool.command)
            print()

    def parse_script(self, match: re.Match[str]) -> XMLScript:
        filename = match.group(1)
        varname = filename.rsplit('.', 1)[0]
        if not varname.endswith('_script'):
            varname = f'{varname}_script'
        script_contents = self.read_script_contents(filename)
        return XMLScript(varname, filename, script_contents)  
    
    def read_script_contents(self, filename: str) -> str:
        xmldir = os.path.dirname(self.xmlpath)
        filepath = os.path.join(xmldir, filename)
        with open(filepath, 'r') as f:
            script_contents = f.read()
        return script_contents

    def replace_script_in_command(self, match: re.Match[str], param: XMLDataParam) -> None:
        self.gxtool.command = self.gxtool.command.replace(match.group(0), f'${param.name}')
    
    def parse_tests(self) -> XMLTestRegister:
        """
        returns a formatted list of tests using the representation
        needs to be properly fleshed out later!
        """
        xmltests: list[XMLTest] = []
        if self.gxtool.tests:
            for test in self.gxtool.tests:
                inputs_dict = to_inputs_dict(test.inputs)
                xmltest = XMLTest(test.name, inputs_dict)
                xmltests.append(xmltest)
        return XMLTestRegister(xmltests)
    


### HELPER FUNCTIONS ###

class RepeatSection(dict):
    def __init__(self):
        super().__init__()

def to_inputs_dict(gxinputs: dict[str, Any]) -> dict[str, Any]:
    inputs_dict: dict[str, Any] = {}
    for longname, value in gxinputs.items():
        inputs_dict = generate_structure(inputs_dict, longname)
        inputs_dict = inject_value(inputs_dict, longname, value)
    inputs_dict = collapse_repeats(inputs_dict)
    return inputs_dict

def generate_structure(inputs_dict: dict[str, Any], longname: str) -> dict[str, Any]:
    REPEAT_MATCHER = r'_\d+$'
    longname = longname.replace('|', '.')
    heirarchy = longname.split('.')

    if len(heirarchy) >= 1:
        heirarchy = heirarchy[:-1]  # don't care about leaf nodes
        node = inputs_dict
        # prepare dict structure & inject value at leaf
        for label in heirarchy:
            # if internal repeat section, handle
            if re.search(REPEAT_MATCHER, label):
                label, elem = label.rsplit('_', 1)
                if label not in node:
                    node[label] = RepeatSection()
                if elem not in node[label]:
                    node[label][elem] = {}
                node = node[label][elem]

            # if internal generic node, handle
            elif label not in node:
                node[label] = {}
                node = node[label]

            else:
                node = node[label]
            
    return inputs_dict

def inject_value(inputs_dict: dict[str, Any], longname: str, value: Any) -> dict[str, Any]:
    REPEAT_MATCHER = r'_\d+$'
    longname = longname.replace('|', '.')
    heirarchy = longname.split('.')
    node = inputs_dict
    
    for i, label in enumerate(heirarchy):
        # inject value if leaf node
        if i == len(heirarchy) - 1:
            # TODO get gxparam, check if it expects list or single
            # value rather than the logic below
            if isinstance(value, list) and len(value) == 1:
                value = value[0]
            node[label] = value
            break
        
        # if internal repeat section, go to next level
        elif re.search(REPEAT_MATCHER, label):
            label, elem = label.rsplit('_', 1)
            node = node[label][elem]
        
        # if internal generic node, go to next level
        else:
            node = node[label]

    return inputs_dict

def collapse_repeats(node: dict[str, Any]) -> dict[str, Any]:
    for key, value in node.items():
        # if repeat, do collapse
        if isinstance(value, RepeatSection):
            repeat = value
            num_elems = len(repeat)
            the_list = [None for _ in range(num_elems)]
            for elem, subdict in repeat.items():
                the_list[int(elem)] = collapse_repeats(subdict)  # type: ignore
            node[key] = the_list
        # if normal dict, descend
        elif isinstance(value, dict):
            node[key] = collapse_repeats(value)
    return node
            
