
import os
import tempfile
import regex as re

from galaxy.tools import Tool as GxTool
from galaxy.tool_util.parser import get_tool_source
from galaxy.tools import create_tool_from_source
from galaxy.model import History

from galaxy.tools import Tool as GxTool
from galaxy.tools.parameters.basic import ToolParameter
from galaxy.tool_util.parser.output_objects import ToolOutput
from .param import DataParam

from ...expressions.patterns import GX_TOOL_SCRIPT
from ...expressions.matches import get_matches
from ..configfiles.Configfile import Configfile
from ..scripts import Script
from .mock import MockApp, MockObjectStore
from .tool import XMLToolDefinition
from .metadata import ToolXMLMetadata
from .citations import Citation
from .parsing.ParamFlattener import ParamFlattener
from .param.ParamRegister import ParamRegister
from .TestRegister import TestRegister
from .requirements.model import CondaRequirement, ContainerRequirement
from .parsing.outputs import parse_output_param
from .parsing.inputs import parse_input_param


Requirement = ContainerRequirement | CondaRequirement


def load_xmltool(path: str) -> XMLToolDefinition:
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
        self.scripts: list[Script] = []
        self.configfiles: list[Configfile] = []
        self.inputs = self.parse_inputs()
        self.parse_scripts()
        self.parse_configfiles()
        self.outputs = self.parse_outputs()

    # TODO remove script references in command parsing etc
    # scripts / configfiles treated as gxparam InputParams. 
    # still need way to set default value? 

    def create(self) -> XMLToolDefinition:
        return XMLToolDefinition(
            metadata=self.parse_metadata(),
            raw_command=self.parse_command(),
            configfiles=self.configfiles,
            scripts=self.scripts,
            inputs=self.inputs,
            outputs=self.outputs,
            tests=self.parse_tests()
        ) 
    
    def parse_inputs(self) -> ParamRegister:
        """returns a an InputRegister by reformatting the galaxy tool representation's params."""
        register = ParamRegister()
        g_in_params = self.flatten_params()
        for g_param in g_in_params:
            t_param = parse_input_param(g_param)
            register.add(t_param)
        return register

    def flatten_params(self) -> list[ToolParameter]:
        pf = ParamFlattener(self.gxtool.inputs)
        return pf.flatten()

    def parse_outputs(self) -> ParamRegister:
        """returns a formatted list of outputs using the representation"""
        register = ParamRegister()
        g_out_params: list[ToolOutput] = list(self.gxtool.outputs.values())
        for g_param in g_out_params:
            t_params = parse_output_param(g_param, self.inputs)
            for param in t_params:
                register.add(param)
        return register

    def parse_metadata(self) -> ToolXMLMetadata:
        """returns a formatted Metadata using the representation"""
        requirements: list[Requirement] = self.get_requirements()
        citations: list[Citation] = self.get_citations()
        return ToolXMLMetadata(
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
    
    def get_conda_requirements(self) -> list[CondaRequirement]:
        packages: list[dict[str, str]] = self.gxtool.requirements.to_list() # type: ignore
        return [CondaRequirement(p['name'], p['version']) for p in packages]
    
    def get_container_requirements(self) -> list[ContainerRequirement]:
        containers: list[ContainerDescription] = self.gxtool.containers # type: ignore
        return [ContainerRequirement(c.identifier) for c in containers] # type: ignore

    def get_citations(self) -> list[Citation]:
        citations: list[Citation] = []
        citations += self.get_biotools_citations()
        citations += self.get_doi_citations()
        #citations += self.get_bibtex_citations()
        return citations

    def get_biotools_citations(self) -> list[Citation]: 
        out: list[Citation] = []
        for ref in self.gxtool.xrefs:
            citation = Citation(
                type='biotools',
                text=f"https://bio.tools/{ref['value']}"
            )
            out.append(citation)
        return out

    def get_doi_citations(self) -> list[Citation]:
        out: list[Citation] = []
        for elem in self.gxtool.tool_source.xml_tree.findall('citations'):
            for ref in elem.findall('citation'):
                if ref.attrib['type'] == 'doi':
                    citation = Citation(
                        type='doi',
                        text=f"https://doi.org/{ref.text}"
                    )
                    out.append(citation)
        return out
        
    def get_bibtex_citations(self) -> list[Citation]:
        out: list[Citation] = []
        for elem in self.gxtool.tool_source.xml_tree.findall('citations'):
            for ref in elem.findall('citation'):
                if ref.attrib['type'] == 'bibtex':
                    citation = Citation(
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
    
    def parse_configfiles(self) -> list[Configfile]:
        # TODO HERE SAME AS SCRIPT
        """returns the tool configfiles"""
        out: list[Configfile] = []
        for name, _, contents in self.gxtool.config_files:  # type: ignore
            if isinstance(contents, str):
                new_config = Configfile(name, contents)  # type: ignore
                out.append(new_config)
        if out:
            pass
            # logging.has_configfile()
        return out
    
    def parse_scripts(self) -> None:
        """returns the tools script it calls in the command"""
        while get_matches(self.gxtool.command, GX_TOOL_SCRIPT):
            # get the match
            match = get_matches(self.gxtool.command, GX_TOOL_SCRIPT)[0]
            
            # parse match into script
            script = self.parse_script(match)
            
            # add script to list of scripts
            self.scripts.append(script)

            # add param for script
            param = DataParam(name=script.varname)
            param.formats = ['file']
            param.helptext = 'galaxy script needed to run tool'

            # add param to input register
            self.inputs.add(param)
            
            # modify command to replace script with param
            self.replace_script_in_command(match, param)
            print(self.gxtool.command)
            print()

    def parse_script(self, match: re.Match[str]) -> Script:
        filename = match.group(2)
        varname = filename.rsplit('.', 1)[0]
        if not varname.endswith('_script'):
            varname = f'{varname}_script'
        script_contents = self.read_script_contents(filename)
        return Script(varname, filename, script_contents)  
    
    def read_script_contents(self, filename: str) -> str:
        xmldir = os.path.dirname(self.xmlpath)
        filepath = os.path.join(xmldir, filename)
        with open(filepath, 'r') as f:
            script_contents = f.read()
        return script_contents

    def replace_script_in_command(self, match: re.Match[str], param: DataParam) -> None:
        self.gxtool.command = self.gxtool.command.replace(match.group(0), f'${param.name}')
    
    def parse_tests(self) -> TestRegister:
        """
        returns a formatted list of tests using the representation
        needs to be properly fleshed out later!
        """
        return TestRegister([])
