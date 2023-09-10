
import itertools
import math
import os
import re
import tarfile
import tempfile
import threading
from pathlib import Path
from typing import (
    Any,
    cast,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    Union,
)
from urllib.parse import unquote_plus

import packaging.version
import webob.exc
from mako.template import Template
from webob.compat import cgi_FieldStorage

from galaxy import (
    exceptions,
    model
)
from galaxy.exceptions import ToolInputsNotReadyException
from galaxy.job_execution import output_collect
from galaxy.tool_shed.util.repository_util import get_installed_repository
from galaxy.tool_shed.util.shed_util_common import set_image_paths
from galaxy.tool_util.deps import (
    CachedDependencyManager,
)
from galaxy.tool_util.loader import (
    imported_macro_paths,
    raw_tool_xml_tree,
    template_macro_params
)
from galaxy.tool_util.parser import (
    RequiredFiles,
    ToolOutputCollectionPart
)
from galaxy.tool_util.parser.xml import XmlPageSource
from galaxy.tool_util.provided_metadata import parse_tool_provided_metadata
from galaxy.tools import expressions
from galaxy.tools.actions import DefaultToolAction, ToolAction
from galaxy.tools.evaluation import global_tool_errors
from galaxy.tools.parameters import (
    check_param,
    params_from_strings,
    params_to_incoming,
    params_to_strings,
    populate_state,
    visit_input_values
)
from galaxy.tools.parameters.basic import (
    BaseURLToolParameter,
    DataCollectionToolParameter,
    DataToolParameter,
    HiddenToolParameter,
    ImplicitConversionRequired,
    SelectToolParameter,
    ToolParameter,
    workflow_building_modes,
)
from galaxy.tools.parameters.dataset_matcher import (
    set_dataset_matcher_factory,
    unset_dataset_matcher_factory,
)
from galaxy.tools.parameters.grouping import Conditional, ConditionalWhen, Repeat, Section, UploadDataset
from galaxy.tools.parameters.input_translation import ToolInputTranslator
from galaxy.tools.parameters.meta import expand_meta_parameters
from galaxy.tools.test import parse_tests
from galaxy.util import (
    in_directory,
    Params,
    parse_xml_string,
    rst_to_html,
    string_as_bool,
    unicodify,
    XML,
)
from galaxy.util.expressions import ExpressionContext
from galaxy.util.template import (
    fill_template,
    refactoring_tool,
)
from galaxy.util.tool_shed.common_util import (
    get_tool_shed_repository_url,
    get_tool_shed_url_from_tool_shed_registry,
)
from galaxy.version import VERSION_MAJOR
from galaxy.work.context import proxy_work_context_for_history


HELP_UNINITIALIZED = threading.Lock()

REQUIRE_FULL_DIRECTORY = {
    "includes": [{"path": "**", "path_type": "glob"}],
}
IMPLICITLY_REQUIRED_TOOL_FILES: Dict[str, Dict] = {
    "deseq2": {"version": packaging.version.parse("2.11.40.6"), "required": {"includes": [{"path": "*.R", "path_type": "glob"}]}},
    # minimum example:
    # "foobar": {"required": REQUIRE_FULL_DIRECTORY}
    # if no version is specified, all versions without explicit RequiredFiles will be selected
    "circos": {"required": REQUIRE_FULL_DIRECTORY},
    "cp_image_math": {"required": {"includes": [{"path": "*.py", "path_type": "glob"}]}},
    "enumerate_charges": {"required": REQUIRE_FULL_DIRECTORY},
    "fasta_compute_length": {"required": {"includes": [{"path": "utils/*", "path_type": "glob"}]}},
    "fasta_concatenate0": {"required": {"includes": [{"path": "utils/*", "path_type": "glob"}]}},
    "filter_tabular": {"required": {"includes": [{"path": "*.py", "path_type": "glob"}]}},
    "flanking_features_1": {"required": {"includes": [{"path": "utils/*", "path_type": "glob"}]}},
    "gops_intersect_1": {"required": {"includes": [{"path": "utils/*", "path_type": "glob"}]}},
    "gops_subtract_1": {"required": {"includes": [{"path": "utils/*", "path_type": "glob"}]}},
    "maxquant": {"required": {"includes": [{"path": "*.py", "path_type": "glob"}]}},
    "maxquant_mqpar": {"required": {"includes": [{"path": "*.py", "path_type": "glob"}]}},
    "query_tabular": {"required": {"includes": [{"path": "*.py", "path_type": "glob"}]}},
    "shasta": {"required": {"includes": [{"path": "configs/*", "path_type": "glob"}]}},
    "sqlite_to_tabular": {"required": {"includes": [{"path": "*.py", "path_type": "glob"}]}},
    "sucos_max_score": {"required": {"includes": [{"path": "*.py", "path_type": "glob"}]}},
}

from galaxy.datatypes.registry import Registry


class App:

    def __init__(self) -> None:
        self.datatypes_registry = Registry()
        self.datatypes_registry.load_datatypes()


class GXTool:
    """
    Represents a computational tool that can be executed through Galaxy.
    """

    tool_type = 'default'
    requires_setting_metadata = True
    produces_entry_points = False
    default_tool_action = DefaultToolAction
    tool_action: ToolAction
    tool_type_local = False
    dict_collection_visible_keys = ['id', 'name', 'version', 'description', 'labels']
    __help: Optional[threading.Lock]
    __help_by_page: Union[threading.Lock, List[str]]
    job_search: 'JobSearch'
    version: str

    def __init__(self, config_file, tool_source, guid=None, repository_id=None, tool_shed_repository=None, allow_code_files=True, dynamic=False, tool_dir=None):
        """Load a tool from the config named by `config_file`"""
        # Determine the full path of the directory where the tool config is
        if config_file is not None:
            self.config_file = config_file
            self.tool_dir = tool_dir or os.path.dirname(config_file)
        else:
            self.config_file = None
            self.tool_dir = tool_dir

        # self.app = app
        self.repository_id = repository_id
        self._allow_code_files = allow_code_files
        # setup initial attribute values
        self.stdio_exit_codes = list()
        self.stdio_regexes = list()
        self.inputs_by_page = list()
        self.display_by_page = list()
        self.action: Union[str, Tuple[str, str]] = "/tool_runner/index"
        self.target = "galaxy_main"
        self.method = "post"
        self.labels = []
        self.check_values = True
        self.nginx_upload = False
        self.input_required = False
        self.display_interface = True
        self.require_login = False
        self.rerun = False
        # This will be non-None for tools loaded from the database (DynamicTool objects).
        self.dynamic_tool = None
        # Define a place to keep track of all input   These
        # differ from the inputs dictionary in that inputs can be page
        # elements like conditionals, but input_params are basic form
        # parameters like SelectField objects.  This enables us to more
        # easily ensure that parameter dependencies like index files or
        # tool_data_table_conf.xml entries exist.
        self.input_params = []
        # Attributes of tools installed from Galaxy tool sheds.
        self.tool_shed = None
        self.repository_name = None
        self.repository_owner = None
        self.changeset_revision = None
        self.installed_changeset_revision = None
        self.sharable_url = None
        # The tool.id value will be the value of guid, but we'll keep the
        # guid attribute since it is useful to have.
        self.guid = guid
        self.old_id = None
        self.python_template_version = None
        self._lineage = None
        self.dependencies = []
        # populate toolshed repository info, if available
        self.populate_tool_shed_info(tool_shed_repository)
        # add tool resource parameters
        # self.populate_resource_parameters(tool_source)
        self.tool_errors = None
        # Parse XML element containing configuration
        self.tool_source = tool_source
        self._is_workflow_compatible = None
        self.finalized = False
        try:
            self.parse(tool_source, guid=guid, dynamic=dynamic)
        except Exception as e:
            global_tool_errors.add_error(config_file, "Tool Loading", e)
            raise e
        # # The job search is only relevant in a galaxy context, and breaks
        # # loading tools into the toolshed for validation.
        # if self.app.name == 'galaxy':
        #     self.job_search = self.app.job_search

    def __getattr__(self, name):
        lazy_attributes = {
            'action',
            'check_values',
            'display_by_page',
            'enctype',
            'has_multiple_pages',
            'inputs',
            'inputs_by_page',
            'last_page',
            'method',
            'npages',
            'nginx_upload',
            'target',
            'template_macro_params',
            'outputs',
            'output_collections'
        }
        if name in lazy_attributes:
            self.assert_finalized()
            return getattr(self, name)
        raise AttributeError(name)

    def assert_finalized(self, raise_if_invalid=False):
        if self.finalized is False:
            try:
                self.parse_inputs(self.tool_source)
                self.parse_outputs(self.tool_source)
                self.finalized = True
            except Exception as e:
                toolbox = getattr(self.app, 'toolbox', None)
                if toolbox:
                    toolbox.remove_tool_by_id(self.id)
                if raise_if_invalid:
                    raise
                else:
                    print("An error occured while parsing the tool wrapper xml, the tool is not functional")

    def remove_from_cache(self):
        source_path = self.tool_source._source_path
        if source_path:
            for region in self.app.toolbox.cache_regions.values():
                region.delete(source_path)

    @property
    def history_manager(self):
        return self.app.history_manager

    @property
    def _view(self):
        return self.app.dependency_resolvers_view

    @property
    def version_object(self):
        return packaging.version.parse(self.version)

    @property
    def sa_session(self):
        """Returns a SQLAlchemy session"""
        return self.app.model.context

    @property
    def lineage(self):
        """Return ToolLineage for this tool."""
        return self._lineage

    @property
    def tool_versions(self):
        # If we have versions, return them.
        if self.lineage:
            return list(self.lineage.tool_versions)
        else:
            return []

    @property
    def is_latest_version(self):
        tool_versions = self.tool_versions
        return not tool_versions or self.version == self.tool_versions[-1]

    @property
    def latest_version(self):
        if self.is_latest_version:
            return self
        else:
            return self.app.tool_cache.get_tool_by_id(self.lineage.get_versions()[-1].id)

    @property
    def is_datatype_converter(self):
        return self in self.app.datatypes_registry.converter_tools

    @property
    def tool_shed_repository(self):
        # If this tool is included in an installed tool shed repository, return it.
        if self.tool_shed:
            return get_installed_repository(self.app,
                                            tool_shed=self.tool_shed,
                                            name=self.repository_name,
                                            owner=self.repository_owner,
                                            installed_changeset_revision=self.installed_changeset_revision,
                                            from_cache=True)

    @property
    def produces_collections_with_unknown_structure(self):

        def output_is_dynamic(output):
            if not output.collection:
                return False
            return output.dynamic_structure

        return any(map(output_is_dynamic, self.outputs.values()))

    @property
    def valid_input_states(self):
        return model.Dataset.valid_input_states

    @property
    def requires_galaxy_python_environment(self):
        """Indicates this tool's runtime requires Galaxy's Python environment."""
        # All special tool types (data source, history import/export, etc...)
        # seem to require Galaxy's Python.
        # FIXME: the (instantiated) tool class should emit this behavior, and not
        #        use inspection by string check
        if self.tool_type not in ["default", "manage_data", "interactive", "data_source"]:
            return True

        if self.tool_type == "manage_data" and self.profile < 18.09:
            return True

        if self.tool_type == "data_source" and self.profile < 21.09:
            return True

        config = self.app.config
        preserve_python_environment = config.preserve_python_environment
        if preserve_python_environment == "always":
            return True
        elif preserve_python_environment == "legacy_and_local" and self.tool_shed is None:
            return True
        else:
            unversioned_legacy_tool = self.old_id in GALAXY_LIB_TOOLS_UNVERSIONED
            versioned_legacy_tool = self.old_id in GALAXY_LIB_TOOLS_VERSIONED
            legacy_tool = unversioned_legacy_tool or (
                versioned_legacy_tool
                and self.old_id
                and self.version_object < GALAXY_LIB_TOOLS_VERSIONED[self.old_id]
            )
            return legacy_tool

    def __get_job_tool_configuration(self, job_params=None):
        """Generalized method for getting this tool's job configuration.

        :type job_params: dict or None
        :returns: `galaxy.jobs.JobToolConfiguration` -- JobToolConfiguration that matches this `Tool` and the given `job_params`
        """
        rval = None
        if len(self.job_tool_configurations) == 1:
            # If there's only one config, use it rather than wasting time on comparisons
            rval = self.job_tool_configurations[0]
        elif job_params is None:
            for job_tool_config in self.job_tool_configurations:
                if not job_tool_config.params:
                    rval = job_tool_config
                    break
        else:
            for job_tool_config in self.job_tool_configurations:
                if job_tool_config.params:
                    # There are job params and this config has params defined
                    for param, value in job_params.items():
                        if param not in job_tool_config.params or job_tool_config.params[param] != value:
                            break
                    else:
                        # All params match, use this config
                        rval = job_tool_config
                        break
                else:
                    rval = job_tool_config
        assert rval is not None, f'Could not get a job tool configuration for Tool {self.id} with job_params {job_params}, this is a bug'
        return rval

    def get_configured_job_handler(self, job_params=None):
        """Get the configured job handler for this `Tool` given the provided `job_params`.

        Unlike the former ``get_job_handler()`` method, this does not perform "preassignment" (random selection of
        a configured handler ID from a tag).

        :param job_params: Any params specific to this job (e.g. the job source)
        :type job_params: dict or None

        :returns: str or None -- The configured handler for a job run of this `Tool`
        """
        return self.__get_job_tool_configuration(job_params=job_params).handler

    def get_job_destination(self, job_params=None):
        """
        :returns: galaxy.jobs.JobDestination -- The destination definition and runner parameters.
        """
        return self.app.job_config.get_destination(self.__get_job_tool_configuration(job_params=job_params).destination)

    def get_panel_section(self):
        return self.app.toolbox.get_section_for_tool(self)

    def allow_user_access(self, user, attempting_access=True):
        """
        :returns: bool -- Whether the user is allowed to access the tool.
        """
        if self.require_login and user is None:
            return False
        return True

    def parse(self, tool_source, guid=None, dynamic=False):
        """
        Read tool configuration from the element `root` and fill in `self`.
        """
        self.profile = float(tool_source.parse_profile())
        # Get the UNIQUE id for the tool
        self.old_id = tool_source.parse_id()
        self.id = self.old_id  # added
        # if guid is None:
        #     self.id = self.old_id
        # else:
        #     self.id = guid

        # if not dynamic and not self.id:
        #     raise Exception(f"Missing tool 'id' for tool at '{tool_source}'")

        # profile = packaging.version.parse(str(self.profile))
        # if self.app.name == 'galaxy' and profile >= packaging.version.parse("16.04") and packaging.version.parse(VERSION_MAJOR) < profile:
        #     message = f"The tool [{self.id}] targets version {self.profile} of Galaxy, you should upgrade Galaxy to ensure proper functioning of this tool."
        #     raise Exception(message)

        # self.python_template_version = tool_source.parse_python_template_version()
        # if self.python_template_version is None:
        #     # If python_template_version not specified we assume tools with profile versions >= 19.05 are python 3 ready
        #     if self.profile >= 19.05:
        #         self.python_template_version = packaging.version.parse('3.5')
        #     else:
        #         self.python_template_version = packaging.version.parse('2.7')

        # Get the (user visible) name of the tool
        self.name = tool_source.parse_name()
        if not self.name:
            raise RuntimeError
        # if not self.name and dynamic:
        #     self.name = self.id
        # if not dynamic and not self.name:
        #     raise Exception(f"Missing tool 'name' for tool with id '{self.id}' at '{tool_source}'")

        self.version = tool_source.parse_version()
        if not self.version:
            raise RuntimeError
            # if self.profile < 16.04:
            #     # For backward compatibility, some tools may not have versions yet.
            #     self.version = "1.0.0"
            # else:
            #     raise Exception(f"Missing tool 'version' for tool with id '{self.id}' at '{tool_source}'")

        # Legacy feature, ignored by UI.
        self.force_history_refresh = False

        self.display_interface = False # added
        # self.display_interface = tool_source.parse_display_interface(default=self.display_interface)

        self.require_login = False # added
        # self.require_login = tool_source.parse_require_login(self.require_login)

        self.input_translator = None # added
        # request_param_translation_elem = tool_source.parse_request_param_translation_elem()
        # if request_param_translation_elem is not None:
        #     # Load input translator, used by datasource tools to change names/values of incoming parameters
        #     self.input_translator = ToolInputTranslator.from_element(request_param_translation_elem)
        # else:
        #     self.input_translator = None

        self.parse_command(tool_source)
        self.environment_variables = self.parse_environment_variables(tool_source)
        self.tmp_directory_vars = tool_source.parse_tmp_directory_vars()

        self.home_target = None # added
        self.tmp_target = None # added

        # home_target = tool_source.parse_home_target()
        # tmp_target = tool_source.parse_tmp_target()
        # # If a tool explicitly sets one of these variables just respect that and turn off
        # # explicit processing by Galaxy.
        # for environment_variable in self.environment_variables:
        #     if environment_variable.get("name") == "HOME":
        #         home_target = None
        #         continue
        #     for tmp_directory_var in self.tmp_directory_vars:
        #         if environment_variable.get("name") == tmp_directory_var:
        #             tmp_target = None
        #             break
        # self.home_target = home_target
        # self.tmp_target = tmp_target
        self.docker_env_pass_through = [] # added
        # self.docker_env_pass_through = tool_source.parse_docker_env_pass_through()
        # if self.environment_variables:
        #     if not self.docker_env_pass_through:
        #         self.docker_env_pass_through = []
        #     self.docker_env_pass_through.extend(map(lambda x: x['name'], self.environment_variables))

        # Parameters used to build URL for redirection to external app
        self.redirect_url_params = '' # added
        # redirect_url_params = tool_source.parse_redirect_url_params_elem()
        # if redirect_url_params is not None and redirect_url_params.text is not None:
        #     # get rid of leading / trailing white space
        #     redirect_url_params = redirect_url_params.text.strip()
        #     # Replace remaining white space with something we can safely split on later
        #     # when we are building the params
        #     self.redirect_url_params = redirect_url_params.replace(' ', '**^**')
        # else:
        #     self.redirect_url_params = ''

        # Short description of the tool
        self.description = tool_source.parse_description()

        # Versioning for tools
        self.version_string_cmd = None
        # version_command = tool_source.parse_version_command()
        # if version_command is not None:
        #     self.version_string_cmd = version_command.strip()

        #     version_cmd_interpreter = tool_source.parse_version_command_interpreter()
        #     if version_cmd_interpreter:
        #         executable = self.version_string_cmd.split()[0]
        #         abs_executable = os.path.abspath(os.path.join(self.tool_dir, executable))
        #         command_line = self.version_string_cmd.replace(executable, abs_executable, 1)
        #         self.version_string_cmd = f"{version_cmd_interpreter} {command_line}"

        # Parallelism for tasks, read from tool config.
        self.parallelism = None # added
        # self.parallelism = tool_source.parse_parallelism()

        self.all_ids = [] # added
        # # Get JobToolConfiguration(s) valid for this particular Tool.  At least
        # # a 'default' will be provided that uses the 'default' handler and
        # # 'default' destination.  I thought about moving this to the
        # # job_config, but it makes more sense to store here. -nate
        # if self.id:
        #     self_ids = [self.id.lower()]
        #     if self.old_id != self.id:
        #         # Handle toolshed guids
        #         self_ids = [self.id.lower(), self.id.lower().rsplit('/', 1)[0], self.old_id.lower()]
        # else:
        #     self_ids = []
        # self.all_ids = self_ids

        # In the toolshed context, there is no job config.
        self.job_tool_configurations = None # added
        # if hasattr(self.app, 'job_config'):
        #     # Order of this list must match documentation in job_conf.sample_advanced.yml
        #     tool_classes = []
        #     if self.tool_type_local:
        #         tool_classes.append("local")
        #     elif self.old_id in ['upload1', '__DATA_FETCH__']:
        #         tool_classes.append("local")
        #     if self.requires_galaxy_python_environment:
        #         tool_classes.append("requires_galaxy")

        #     self.job_tool_configurations = self.app.job_config.get_job_tool_configurations(self_ids, tool_classes)

        # Is this a 'hidden' tool (hidden in tool menu)
        self.hidden = tool_source.parse_hidden()
        self.license = tool_source.parse_license()
        self.creator = tool_source.parse_creator()

        self.__parse_legacy_features(tool_source)

        # Load any tool specific options (optional)
        self.options = None # added
        # self.options = _Options(
        #     **dict(
        #         sanitize=tool_source.parse_sanitize(),
        #         refresh=tool_source.parse_refresh(),
        #     )
        # )

        # Read in name of galaxy.json metadata file and how to parse it.
        self.provided_metadata_file = tool_source.parse_provided_metadata_file()
        self.provided_metadata_style = tool_source.parse_provided_metadata_style()

        # Parse tool help
        self.parse_help(tool_source)

        # Parse result handling for tool exit codes and stdout/stderr messages:
        self.parse_stdio(tool_source)

        self.strict_shell = tool_source.parse_strict_shell()

        # Any extra generated config files for the tool
        self.__parse_config_files(tool_source)
        # Action
        self.tool_action = self.default_tool_action() # added
        # action = tool_source.parse_action_module()
        # if action is None:
        #     self.tool_action = self.default_tool_action()
        # else:
        #     module, cls = action
        #     mod = __import__(module, globals(), locals(), [cls])
        #     self.tool_action = getattr(mod, cls)()
        #     if getattr(self.tool_action, "requires_js_runtime", False):
        #         try:
        #             expressions.find_engine(self.app.config)
        #         except Exception:
        #             message = REQUIRES_JS_RUNTIME_MESSAGE % self.tool_id or self.tool_uuid
        #             raise Exception(message)
        # Tests
        self.__parse_tests(tool_source)

        # Requirements (dependencies)
        requirements, containers = tool_source.parse_requirements_and_containers()
        self.requirements = requirements
        self.containers = containers

        required_files = tool_source.parse_required_files()
        if required_files is None:
            old_id = self.old_id
            if old_id in IMPLICITLY_REQUIRED_TOOL_FILES:
                lineage_requirement = IMPLICITLY_REQUIRED_TOOL_FILES[old_id]
                lineage_requirement_until = lineage_requirement.get("version")
                if lineage_requirement_until is None or self.version_object < lineage_requirement_until:
                    required_files = RequiredFiles.from_dict(lineage_requirement["required"])
        self.required_files = required_files

        self.citations = [] # added
        # self.citations = self._parse_citations(tool_source)
        xrefs = tool_source.parse_xrefs()
        has_biotools_reference = any(x["reftype"] == "bio.tools" for x in xrefs)
        if not has_biotools_reference:
            legacy_biotools_ref = self.legacy_biotools_external_reference
            if legacy_biotools_ref is not None:
                xrefs.append({"value": legacy_biotools_ref, "reftype": "bio.tools"})
        self.xrefs = xrefs


        self.edam_operations = []  # added
        self.edam_topics = []  # added
        # edam_operations = tool_source.parse_edam_operations()
        # edam_topics = tool_source.parse_edam_topics()

        # has_missing_data = len(edam_operations) == 0 or len(edam_topics) == 0
        # if has_missing_data:
        #     biotools_reference = self.biotools_reference
        #     metadata_source = self.app.biotools_metadata_source
        #     if biotools_reference and metadata_source:
        #         biotools_entry = metadata_source.get_biotools_metadata(biotools_reference)
        #         if biotools_entry:
        #             edam_info = biotools_entry.edam_info
        #             if len(edam_operations) == 0:
        #                 edam_operations = edam_info.edam_operations
        #             if len(edam_topics) == 0:
        #                 edam_topics = edam_info.edam_topics

        # self.edam_operations = edam_operations
        # self.edam_topics = edam_topics

        self.__parse_trackster_conf(tool_source)
        # Record macro paths so we can reload a tool if any of its macro has changes
        self._macro_paths = tool_source.macro_paths
        self.ports = tool_source.parse_interactivetool()

    def __parse_legacy_features(self, tool_source):
        self.code_namespace: Dict[str, str] = {}
        self.hook_map: Dict[str, str] = {}
        self.uihints: Dict[str, str] = {}

        if not hasattr(tool_source, 'root'):
            return

        # TODO: Move following logic into XmlToolSource.
        root = tool_source.root
        # Load any tool specific code (optional) Edit: INS 5/29/2007,
        # allow code files to have access to the individual tool's
        # "module" if it has one.  Allows us to reuse code files, etc.
        for code_elem in root.findall("code"):
            for hook_elem in code_elem.findall("hook"):
                for key, value in hook_elem.items():
                    # map hook to function
                    self.hook_map[key] = value
            file_name = code_elem.get("file")
            code_path = os.path.join(self.tool_dir, file_name)
            if self._allow_code_files:
                with open(code_path) as f:
                    code_string = f.read()
                try:
                    compiled_code = compile(code_string, code_path, 'exec')
                    exec(compiled_code, self.code_namespace)
                except Exception:
                    if (
                        refactoring_tool
                        and self.python_template_version
                        and self.python_template_version.release[0] < 3
                    ):
                        # Could be a code file that uses python 2 syntax
                        translated_code = str(refactoring_tool.refactor_string(code_string, name='auto_translated_code_file'))
                        compiled_code = compile(translated_code, f"futurized_{code_path}", 'exec')
                        exec(compiled_code, self.code_namespace)
                    else:
                        raise

        # User interface hints
        uihints_elem = root.find("uihints")
        if uihints_elem is not None:
            for key, value in uihints_elem.attrib.items():
                self.uihints[key] = value

    def __parse_tests(self, tool_source):
        self.__tests_source = tool_source
        self.__tests_populated = False

    def __parse_config_files(self, tool_source):
        self.config_files = []
        if not hasattr(tool_source, 'root'):
            return

        root = tool_source.root
        conf_parent_elem = root.find("configfiles")
        if conf_parent_elem is not None:
            inputs_elem = conf_parent_elem.find("inputs")
            if inputs_elem is not None:
                name = inputs_elem.get("name")
                filename = inputs_elem.get("filename", None)
                format = inputs_elem.get("format", "json")
                data_style = inputs_elem.get("data_style", "skip")
                content = dict(format=format, handle_files=data_style, type="inputs")
                self.config_files.append((name, filename, content))
            file_sources_elem = conf_parent_elem.find("file_sources")
            if file_sources_elem is not None:
                name = file_sources_elem.get("name")
                filename = file_sources_elem.get("filename", None)
                content = dict(type="files")
                self.config_files.append((name, filename, content))
            for conf_elem in conf_parent_elem.findall("configfile"):
                name = conf_elem.get("name")
                filename = conf_elem.get("filename", None)
                content = conf_elem.text
                self.config_files.append((name, filename, content))

    def __parse_trackster_conf(self, tool_source):
        self.trackster_conf = None
        if not hasattr(tool_source, 'root'):
            return

        # Trackster configuration.
        trackster_conf = tool_source.root.find("trackster_conf")
        if trackster_conf is not None:
            self.trackster_conf = TracksterConfig.parse(trackster_conf)

    @property
    def tests(self):
        self.assert_finalized()
        if not self.__tests_populated:
            tests_source = self.__tests_source
            if tests_source:
                try:
                    self.__tests = parse_tests(self, tests_source)
                except Exception:
                    self.__tests = None
                    log.exception("Failed to parse tool tests for tool '%s'", self.id)
            else:
                self.__tests = None
            self.__tests_populated = True
        return self.__tests

    @property
    def _repository_dir(self):
        """If tool shed installed tool, the base directory of the repository installed."""
        repository_base_dir = None

        if getattr(self, 'tool_shed', None):
            tool_dir = Path(self.tool_dir)
            for repo_dir in itertools.chain([tool_dir], tool_dir.parents):
                if repo_dir.name == self.repository_name:
                    return str(repo_dir)
            else:
                log.error(f"Problem finding repository dir for tool '{self.id}'")

        return repository_base_dir

    def test_data_path(self, filename):
        repository_dir = self._repository_dir
        test_data = None
        if repository_dir:
            test_data = self.__walk_test_data(dir=repository_dir, filename=filename)
        else:
            if self.tool_dir:
                tool_dir = self.tool_dir
                if isinstance(self, DataManagerTool):
                    tool_dir = os.path.dirname(self.tool_dir)
                test_data = self.__walk_test_data(tool_dir, filename=filename)
        if not test_data:
            # Fallback to Galaxy test data directory for builtin tools, tools
            # under development, and some older ToolShed published tools that
            # used stock test data.
            test_data = self.app.test_data_resolver.get_filename(filename)
        return test_data

    def __walk_test_data(self, dir, filename):
        for root, dirs, _ in os.walk(dir):
            if '.hg' in dirs:
                dirs.remove('.hg')
            if 'test-data' in dirs:
                test_data_dir = os.path.join(root, 'test-data')
                result = os.path.abspath(os.path.join(test_data_dir, filename))
                if not in_directory(result, test_data_dir):
                    # Don't raise an explicit exception and reveal details about what
                    # files are or are not on the path, simply return None and let the
                    # API raise a 404.
                    return None
                else:
                    if os.path.exists(result):
                        return result

    def tool_provided_metadata(self, job_wrapper):
        meta_file = os.path.join(job_wrapper.tool_working_directory, self.provided_metadata_file)
        return parse_tool_provided_metadata(meta_file, provided_metadata_style=self.provided_metadata_style, job_wrapper=job_wrapper)

    def parse_command(self, tool_source):
        """
        """
        # Command line (template). Optional for tools that do not invoke a local program
        command = tool_source.parse_command()
        if command is not None:
            self.command = command.lstrip()  # get rid of leading whitespace
            # Must pre-pend this AFTER processing the cheetah command template
            self.interpreter = tool_source.parse_interpreter()
        else:
            self.command = ''
            self.interpreter = None

    def parse_environment_variables(self, tool_source):
        return tool_source.parse_environment_variables()

    def parse_inputs(self, tool_source):
        """
        Parse the "<inputs>" element and create appropriate `ToolParameter` s.
        This implementation supports multiple pages and grouping constructs.
        """
        # Load parameters (optional)
        self.inputs = {}
        pages = tool_source.parse_input_pages()
        enctypes: Set[str] = set()
        if pages.inputs_defined:
            if hasattr(pages, "input_elem"):
                input_elem = pages.input_elem
                # Handle properties of the input form
                self.check_values = string_as_bool(input_elem.get("check_values", self.check_values))
                self.nginx_upload = string_as_bool(input_elem.get("nginx_upload", self.nginx_upload))
                self.action = input_elem.get('action', self.action)
                # If we have an nginx upload, save the action as a tuple instead of
                # a string. The actual action needs to get url_for run to add any
                # prefixes, and we want to avoid adding the prefix to the
                # nginx_upload_path.
                # if (
                #     self.nginx_upload
                #     and self.app.config.nginx_upload_path
                #     and not isinstance(self.action, tuple)
                # ):
                #     if "?" in unquote_plus(self.action):
                #         raise Exception(
                #             "URL parameters in a non-default tool action can not be used "
                #             "in conjunction with nginx upload.  Please convert them to "
                #             "hidden POST parameters"
                #         )
                #     self.action = (
                #         f"{self.app.config.nginx_upload_path}?nginx_redir=",
                #         unquote_plus(self.action),
                #     )
                self.target = input_elem.get("target", self.target)
                self.method = input_elem.get("method", self.method)
                # Parse the actual parameters
                # Handle multiple page case
            for page_source in pages.page_sources:
                inputs = self.parse_input_elem(page_source, enctypes)
                display = page_source.parse_display()
                self.inputs_by_page.append(inputs)
                self.inputs.update(inputs)
                self.display_by_page.append(display)
        else:
            self.inputs_by_page.append(self.inputs)
            self.display_by_page.append(None)
        self.display = self.display_by_page[0]
        self.npages = len(self.inputs_by_page)
        self.last_page = len(self.inputs_by_page) - 1
        self.has_multiple_pages = bool(self.last_page)
        # Determine the needed enctype for the form
        if len(enctypes) == 0:
            self.enctype = "application/x-www-form-urlencoded"
        elif len(enctypes) == 1:
            self.enctype = enctypes.pop()
        else:
            raise Exception(f"Conflicting required enctypes: {str(enctypes)}")
        # Check if the tool either has no parameters or only hidden (and
        # thus hardcoded)  FIXME: hidden parameters aren't
        # parameters at all really, and should be passed in a different
        # way, making this check easier.
        template_macros = {}
        if hasattr(tool_source, 'root'):
            template_macros = template_macro_params(tool_source.root)
        self.template_macro_params = template_macros
        for param in self.inputs.values():
            if not isinstance(param, (HiddenToolParameter, BaseURLToolParameter)):
                self.input_required = True
                break

    def parse_help(self, tool_source):
        """
        Parse the help text for the tool. Formatted in reStructuredText, but
        stored as Mako to allow for dynamic image paths.
        This implementation supports multiple pages.
        """
        # TODO: Allow raw HTML or an external link.
        self.__help = HELP_UNINITIALIZED
        self.__help_by_page = HELP_UNINITIALIZED
        self.__help_source = tool_source

    def parse_outputs(self, tool_source):
        """
        Parse <outputs> elements and fill in self.outputs (keyed by name)
        """
        self.outputs, self.output_collections = tool_source.parse_outputs(self)

    # TODO: Include the tool's name in any parsing warnings.
    def parse_stdio(self, tool_source):
        """
        Parse <stdio> element(s) and fill in self.return_codes,
        self.stderr_rules, and self.stdout_rules. Return codes have a range
        and an error type (fault or warning).  Stderr and stdout rules have
        a regular expression and an error level (fault or warning).
        """
        exit_codes, regexes = tool_source.parse_stdio()
        self.stdio_exit_codes = exit_codes
        self.stdio_regexes = regexes

    def _parse_citations(self, tool_source):
        # TODO: Move following logic into ToolSource abstraction.
        if not hasattr(tool_source, 'root'):
            return []

        root = tool_source.root
        citations: List[str] = []
        citations_elem = root.find("citations")
        if citations_elem is None:
            return citations

        for citation_elem in citations_elem:
            if citation_elem.tag != "citation":
                pass
            if hasattr(self.app, 'citations_manager'):
                citation = self.app.citations_manager.parse_citation(citation_elem)
                if citation:
                    citations.append(citation)
        return citations

    def parse_input_elem(self, page_source, enctypes, context=None):
        """
        Parse a parent element whose children are inputs -- these could be
        groups (repeat, conditional) or param elements. Groups will be parsed
        recursively.
        """
        rval: Dict[str, Any] = {}
        context = ExpressionContext(rval, context)
        for input_source in page_source.parse_input_sources():
            # Repeat group
            input_type = input_source.parse_input_type()
            if input_type == "repeat":
                group_r = Repeat()
                group_r.name = input_source.get("name")
                group_r.title = input_source.get("title")
                group_r.help = input_source.get("help", None)
                page_source = input_source.parse_nested_inputs_source()
                group_r.inputs = self.parse_input_elem(page_source, enctypes, context)
                group_r.default = int(input_source.get("default", 0))
                group_r.min = int(input_source.get("min", 0))
                # Use float instead of int so that math.inf can be used for no max
                group_r.max = float(input_source.get("max", math.inf))
                assert group_r.min <= group_r.max, ValueError(
                    f"Tool with id '{self.id}': min repeat count must be less-than-or-equal to the max."
                )
                # Force default to be within min-max range
                group_r.default = cast(
                    int, min(max(group_r.default, group_r.min), group_r.max)
                )
                rval[group_r.name] = group_r
            elif input_type == "conditional":
                group_c = Conditional()
                group_c.name = input_source.get("name")
                group_c.value_ref = input_source.get("value_ref", None)
                group_c.value_ref_in_group = input_source.get_bool(
                    "value_ref_in_group", True
                )
                value_from = input_source.get("value_from", None)
                if value_from:
                    value_from = value_from.split(':')
                    temp_value_from = locals().get(value_from[0])
                    group_c.test_param = rval[group_c.value_ref]
                    group_c.test_param.refresh_on_change = True
                    for attr in value_from[1].split('.'):
                        temp_value_from = getattr(temp_value_from, attr)
                    group_c.value_from = temp_value_from  # type: ignore[assignment]
                    # ^^ due to https://github.com/python/mypy/issues/2427
                    assert group_c.value_from
                    for case_value, case_inputs in group_c.value_from(
                        context, group_c, self
                    ).items():
                        # TODO move this to attribute check to remove galaxy class imports
                        case = ConditionalWhen()
                        case.value = case_value
                        if case_inputs:
                            page_source = XmlPageSource(XML(f"<when>{case_inputs}</when>"))
                            case.inputs = self.parse_input_elem(page_source, enctypes, context)
                        else:
                            case.inputs = {}
                        group_c.cases.append(case)
                else:
                    # Should have one child "input" which determines the case
                    test_param_input_source = input_source.parse_test_input_source()
                    group_c.test_param = self.parse_param_elem(
                        test_param_input_source, enctypes, context
                    )
                    if group_c.test_param.optional:
                        log.debug(f"Tool with id '{self.id}': declares a conditional test parameter as optional, this is invalid and will be ignored.")
                        group_c.test_param.optional = False
                    possible_cases = list(
                        group_c.test_param.legal_values
                    )  # store possible cases, undefined whens will have no inputs
                    # Must refresh when test_param changes
                    group_c.test_param.refresh_on_change = True
                    # And a set of possible cases
                    for (value, case_inputs_source) in input_source.parse_when_input_sources():
                        # TODO move this to attribute check to remove galaxy class imports
                        case = ConditionalWhen()
                        case.value = value
                        case.inputs = self.parse_input_elem(case_inputs_source, enctypes, context)
                        group_c.cases.append(case)
                        try:
                            possible_cases.remove(case.value)
                        except Exception:
                            log.debug(
                                "Tool with id '%s': a when tag has been defined for '%s (%s) --> %s', but does not appear to be selectable."
                                % (
                                    self.id,
                                    group_c.name,
                                    group_c.test_param.name,
                                    case.value,
                                )
                            )
                    for unspecified_case in possible_cases:
                        log.warning(
                            "Tool with id '%s': a when tag has not been defined for '%s (%s) --> %s', assuming empty inputs."
                            % (
                                self.id,
                                group_c.name,
                                group_c.test_param.name,
                                unspecified_case,
                            )
                        )
                        # TODO move this to attribute check to remove galaxy class imports
                        case = ConditionalWhen()
                        case.value = unspecified_case
                        case.inputs = {}
                        group_c.cases.append(case)
                rval[group_c.name] = group_c
            elif input_type == "section":
                group_s = Section()
                group_s.name = input_source.get("name")
                group_s.title = input_source.get("title")
                group_s.help = input_source.get("help", None)
                group_s.expanded = input_source.get_bool("expanded", False)
                page_source = input_source.parse_nested_inputs_source()
                group_s.inputs = self.parse_input_elem(page_source, enctypes, context)
                rval[group_s.name] = group_s
            elif input_type == "upload_dataset":
                elem = input_source.elem()
                group_u = UploadDataset()
                group_u.name = elem.get("name")
                group_u.title = elem.get("title")
                group_u.file_type_name = elem.get(
                    "file_type_name", group_u.file_type_name
                )
                group_u.default_file_type = elem.get(
                    "default_file_type", group_u.default_file_type
                )
                group_u.metadata_ref = elem.get("metadata_ref", group_u.metadata_ref)
                try:
                    rval[group_u.file_type_name].refresh_on_change = True
                except KeyError:
                    pass
                group_page_source = XmlPageSource(elem)
                group_u.inputs = self.parse_input_elem(
                    group_page_source, enctypes, context
                )
                rval[group_u.name] = group_u
            elif input_type == "param":
                param = self.parse_param_elem(input_source, enctypes, context)
                rval[param.name] = param
                if hasattr(param, 'data_ref'):
                    param.ref_input = context[param.data_ref]
                self.input_params.append(param)
        return rval

    def parse_param_elem(self, input_source, enctypes, context):
        """
        Parse a single "<param>" element and return a ToolParameter instance.
        Also, if the parameter has a 'required_enctype' add it to the set
        enctypes.
        """
        param = ToolParameter.build(self, input_source)
        param_enctype = param.get_required_enctype()
        if param_enctype:
            enctypes.add(param_enctype)
        # If parameter depends on any other paramters, we must refresh the
        # form when it changes
        for name in param.get_dependencies():
            # Let it throw exception, but give some hint what the problem might be
            if name not in context:
                log.error(f"Tool with id '{self.id}': Could not find dependency '{name}' of parameter '{param.name}'")
            context[name].refresh_on_change = True
        return param

    def populate_resource_parameters(self, tool_source):
        root = getattr(tool_source, 'root', None)
        if root is not None and hasattr(self.app, 'job_config') and hasattr(self.app.job_config, 'get_tool_resource_xml'):
            resource_xml = self.app.job_config.get_tool_resource_xml(root.get('id', '').lower(), self.tool_type)
            if resource_xml is not None:
                inputs = root.find('inputs')
                if inputs is None:
                    inputs = parse_xml_string('<inputs/>')
                    root.append(inputs)
                inputs.append(resource_xml)

    def populate_tool_shed_info(self, tool_shed_repository):
        if tool_shed_repository:
            self.tool_shed = tool_shed_repository.tool_shed
            self.repository_name = tool_shed_repository.name
            self.repository_owner = tool_shed_repository.owner
            self.changeset_revision = tool_shed_repository.changeset_revision
            self.installed_changeset_revision = tool_shed_repository.installed_changeset_revision
            self.sharable_url = get_tool_shed_repository_url(
                self.app, self.tool_shed, self.repository_owner, self.repository_name
            )

    @property
    def legacy_biotools_external_reference(self) -> Optional[str]:
        """Return a bio.tools ID if any of tool's IDs are BIOTOOLS_MAPPING."""
        for tool_id in self.all_ids:
            if tool_id in BIOTOOLS_MAPPING:
                return BIOTOOLS_MAPPING[tool_id]
        return None

    @property
    def biotools_reference(self) -> Optional[str]:
        """Return a bio.tools ID if external reference to it is found.

        If multiple bio.tools references are found, return just the first one.
        """
        for xref in self.xrefs:
            if xref["reftype"] == "bio.tools":
                return xref["value"]
        return None

    @property
    def help(self):
        if self.__help is HELP_UNINITIALIZED:
            self.__ensure_help()
        return self.__help

    @property
    def help_by_page(self):
        if self.__help_by_page is HELP_UNINITIALIZED:
            self.__ensure_help()
        return self.__help_by_page

    @property
    def raw_help(self):
        # may return rst (or Markdown in the future)
        tool_source = self.__help_source
        help_text = tool_source.parse_help()
        return help_text

    def __ensure_help(self):
        with HELP_UNINITIALIZED:
            if self.__help is HELP_UNINITIALIZED:
                self.__inititalize_help()

    def __inititalize_help(self):
        tool_source = self.__help_source
        self.__help = None
        __help_by_page = []
        help_footer = ""
        help_text = tool_source.parse_help()
        if help_text is not None:
            try:
                if help_text.find('.. image:: ') >= 0 and (self.tool_shed_repository or self.repository_id):
                    help_text = set_image_paths(
                        self.app, help_text, encoded_repository_id=self.repository_id, tool_shed_repository=self.tool_shed_repository, tool_id=self.old_id, tool_version=self.version
                    )
            except Exception:
                log.exception("Exception in parse_help, so images may not be properly displayed for tool with id '%s'", self.id)
            try:
                self.__help = Template(rst_to_html(help_text), input_encoding='utf-8',
                                       default_filters=['decode.utf8'],
                                       encoding_errors='replace')
            except Exception:
                log.exception("Exception while parsing help for tool with id '%s'", self.id)

            # Handle deprecated multi-page help text in XML case.
            if hasattr(tool_source, "root"):
                help_elem = tool_source.root.find("help")
                help_header = help_text
                help_pages = help_elem.findall("page")
                # Multiple help page case
                if help_pages:
                    for help_page in help_pages:
                        __help_by_page.append(help_page.text)
                        help_footer = help_footer + help_page.tail
                # Each page has to rendered all-together because of backreferences allowed by rst
                try:
                    __help_by_page = [
                        Template(
                            rst_to_html(help_header + x + help_footer),
                            input_encoding="utf-8",
                            default_filters=["decode.utf8"],
                            encoding_errors="replace",
                        )
                        for x in __help_by_page
                    ]
                except Exception:
                    log.exception("Exception while parsing multi-page help for tool with id '%s'", self.id)
        # Pad out help pages to match npages ... could this be done better?
        while len(__help_by_page) < self.npages:
            __help_by_page.append(self.__help)
        self.__help_by_page = __help_by_page

    def find_output_def(self, name):
        # name is JobToOutputDatasetAssociation name.
        # TODO: to defensive, just throw IndexError and catch somewhere
        # up that stack.
        if ToolOutputCollectionPart.is_named_collection_part_name(name):
            collection_name, part = ToolOutputCollectionPart.split_output_name(name)
            collection_def = self.output_collections.get(collection_name, None)
            if not collection_def:
                return None
            return collection_def.outputs.get(part, None)
        else:
            return self.outputs.get(name, None)

    @property
    def is_workflow_compatible(self):
        is_workflow_compatible = self._is_workflow_compatible
        if is_workflow_compatible is None:
            is_workflow_compatible = self.check_workflow_compatible(self.tool_source)
            if self.finalized:
                self._is_workflow_compatible = is_workflow_compatible
        return is_workflow_compatible

    def check_workflow_compatible(self, tool_source):
        """
        Determine if a tool can be used in workflows. External tools and the
        upload tool are currently not supported by workflows.
        """
        # Multiple page tools are not supported -- we're eliminating most
        # of these anyway
        if self.finalized and self.has_multiple_pages:
            return False
        # This is probably the best bet for detecting external web tools
        # right now
        if self.tool_type.startswith('data_source'):
            return False

        if hasattr(tool_source, "root"):
            root = tool_source.root
            if not string_as_bool(root.get("workflow_compatible", "True")):
                return False

        # TODO: Anyway to capture tools that dynamically change their own
        #       outputs?
        return True

    def new_state(self, trans):
        """
        Create a new `DefaultToolState` for this tool. It will be initialized
        with default values for inputs. Grouping elements are filled in recursively.
        """
        state = DefaultToolState()
        state.initialize(trans, self)
        return state

    def get_param(self, key):
        """
        Returns the parameter named `key` or None if there is no such
        parameter.
        """
        return self.inputs.get(key, None)

    def get_hook(self, name):
        """
        Returns an object from the code file referenced by `code_namespace`
        (this will normally be a callable object)
        """
        if self.code_namespace:
            # Try to look up hook in self.hook_map, otherwise resort to default
            if name in self.hook_map and self.hook_map[name] in self.code_namespace:
                return self.code_namespace[self.hook_map[name]]
            elif name in self.code_namespace:
                return self.code_namespace[name]
        return None

    def visit_inputs(self, values, callback):
        """
        Call the function `callback` on each parameter of this tool. Visits
        grouping parameters recursively and constructs unique prefixes for
        each nested set of  The callback method is then called as:

        `callback( level_prefix, parameter, parameter_value )`
        """
        # HACK: Yet another hack around check_values -- WHY HERE?
        if self.check_values:
            visit_input_values(self.inputs, values, callback)

    def expand_incoming(self, trans, incoming, request_context, input_format='legacy'):
        rerun_remap_job_id = None
        if 'rerun_remap_job_id' in incoming:
            try:
                rerun_remap_job_id = trans.app.security.decode_id(incoming['rerun_remap_job_id'])
            except Exception as exception:
                log.error(str(exception))
                raise exceptions.MessageException("Failure executing tool with id '%s' (attempting to rerun invalid job).", self.id)

        set_dataset_matcher_factory(request_context, self)

        # Fixed set of input parameters may correspond to any number of jobs.
        # Expand these out to individual parameters for given jobs (tool executions).
        expanded_incomings, collection_info = expand_meta_parameters(trans, self, incoming)

        # Remapping a single job to many jobs doesn't make sense, so disable
        # remap if multi-runs of tools are being used.
        if rerun_remap_job_id and len(expanded_incomings) > 1:
            raise exceptions.MessageException(
                "Failure executing tool with id '%s' (cannot create multiple jobs when remapping existing job).", self.id)

        # Process incoming data
        validation_timer = self.app.execution_timer_factory.get_timer(
            'internals.galaxy.tools.validation',
            'Validated and populated state for tool request',
        )
        all_errors = []
        all_params = []
        for expanded_incoming in expanded_incomings:
            params = {}
            errors: Dict[str, str] = {}
            if self.input_translator:
                self.input_translator.translate(expanded_incoming)
            if not self.check_values:
                # If `self.check_values` is false we don't do any checking or
                # processing on input  This is used to pass raw values
                # through to/from external sites.
                params = expanded_incoming
            else:
                # Update state for all inputs on the current page taking new
                # values from `incoming`.
                populate_state(request_context, self.inputs, expanded_incoming, params, errors, simple_errors=False, input_format=input_format)
                # If the tool provides a `validate_input` hook, call it.
                validate_input = self.get_hook('validate_input')
                if validate_input:
                    validate_input(request_context, errors, params, self.inputs)
            all_errors.append(errors)
            all_params.append(params)
        unset_dataset_matcher_factory(request_context)

        log.info(validation_timer)
        return all_params, all_errors, rerun_remap_job_id, collection_info

    def handle_input(self, trans, incoming, history=None, use_cached_job=False, input_format='legacy'):
        """
        Process incoming parameters for this tool from the dict `incoming`,
        update the tool state (or create if none existed), and either return
        to the form or execute the tool (only if 'execute' was clicked and
        there were no errors).
        """
        request_context = proxy_work_context_for_history(trans, history=history)
        all_params, all_errors, rerun_remap_job_id, collection_info = self.expand_incoming(trans=trans, incoming=incoming, request_context=request_context, input_format=input_format)
        # If there were errors, we stay on the same page and display them
        if any(all_errors):
            # simple param_key -> message string for tool form.
            err_data = {key: unicodify(value) for d in all_errors for (key, value) in d.items()}
            param_errors = {}
            for d in all_errors:
                for key, value in d.items():
                    if hasattr(value, 'to_dict'):
                        value_obj = value.to_dict()
                    else:
                        value_obj = {"message": unicodify(value)}
                    param_errors[key] = value_obj
            raise exceptions.RequestParameterInvalidException(', '.join(msg for msg in err_data.values()), err_data=err_data, param_errors=param_errors)
        else:
            mapping_params = MappingParameters(incoming, all_params)
            completed_jobs = {}
            for i, param in enumerate(all_params):
                if use_cached_job:
                    completed_jobs[i] = self.job_search.by_tool_input(
                        trans=trans,
                        tool_id=self.id,
                        tool_version=self.version,
                        param=param,
                        param_dump=self.params_to_strings(param, self.app, nested=True),
                        job_state=None,
                    )
                else:
                    completed_jobs[i] = None
            execution_tracker = execute_job(trans, self, mapping_params, history=request_context.history, rerun_remap_job_id=rerun_remap_job_id, collection_info=collection_info, completed_jobs=completed_jobs)
            # Raise an exception if there were jobs to execute and none of them were submitted,
            # if at least one is submitted or there are no jobs to execute - return aggregate
            # information including per-job errors. Arguably we should just always return the
            # aggregate information - we just haven't done that historically.
            raise_execution_exception = not execution_tracker.successful_jobs and len(all_params) > 0

            if raise_execution_exception:
                raise exceptions.MessageException(execution_tracker.execution_errors[0])

            return dict(out_data=execution_tracker.output_datasets,
                        num_jobs=len(execution_tracker.successful_jobs),
                        job_errors=execution_tracker.execution_errors,
                        jobs=execution_tracker.successful_jobs,
                        output_collections=execution_tracker.output_collections,
                        implicit_collections=execution_tracker.implicit_collections)

    def handle_single_execution(self, trans, rerun_remap_job_id, execution_slice, history, execution_cache=None, completed_job=None, collection_info=None, job_callback=None, flush_job=True):
        """
        Return a pair with whether execution is successful as well as either
        resulting output data or an error message indicating the problem.
        """
        try:
            rval = self.execute(
                trans,
                incoming=execution_slice.param_combination,
                history=history,
                rerun_remap_job_id=rerun_remap_job_id,
                execution_cache=execution_cache,
                dataset_collection_elements=execution_slice.dataset_collection_elements,
                completed_job=completed_job,
                collection_info=collection_info,
                job_callback=job_callback,
                flush_job=flush_job,
            )
            job = rval[0]
            out_data = rval[1]
            if len(rval) > 2:
                execution_slice.history = rval[2]
        except (webob.exc.HTTPFound, exceptions.MessageException) as e:
            # if it's a webob redirect exception, pass it up the stack
            raise e
        except ToolInputsNotReadyException as e:
            return False, e
        except Exception as e:
            log.exception("Exception caught while attempting to execute tool with id '%s':", self.id)
            message = f"Error executing tool with id '{self.id}': {unicodify(e)}"
            return False, message
        if isinstance(out_data, dict):
            return job, list(out_data.items())
        else:
            if isinstance(out_data, str):
                message = out_data
            else:
                message = f"Failure executing tool with id '{self.id}' (invalid data returned from tool execution)"
            return False, message

    def find_fieldstorage(self, x):
        if isinstance(x, cgi_FieldStorage):
            raise InterruptedUpload(None)
        elif isinstance(x, dict):
            [self.find_fieldstorage(y) for y in x.values()]
        elif isinstance(x, list):
            [self.find_fieldstorage(y) for y in x]

    @property
    def params_with_missing_data_table_entry(self):
        """
        Return all parameters that are dynamically generated select lists whose
        options require an entry not currently in the tool_data_table_conf.xml file.
        """
        params = []
        for input_param in self.input_params:
            if isinstance(input_param, SelectToolParameter) and input_param.is_dynamic:
                options = input_param.options
                if options and options.missing_tool_data_table_name and input_param not in params:
                    params.append(input_param)
        return params

    @property
    def params_with_missing_index_file(self):
        """
        Return all parameters that are dynamically generated
        select lists whose options refer to a  missing .loc file.
        """
        params = []
        for input_param in self.input_params:
            if isinstance(input_param, SelectToolParameter) and input_param.is_dynamic:
                options = input_param.options
                if options and options.tool_data_table and options.tool_data_table.missing_index_file and input_param not in params:
                    params.append(input_param)
        return params

    def get_static_param_values(self, trans):
        """
        Returns a map of parameter names and values if the tool does not
        require any user input. Will raise an exception if any parameter
        does require input.
        """
        args = dict()
        for key, param in self.inputs.items():
            # BaseURLToolParameter is now a subclass of HiddenToolParameter, so
            # we must check if param is a BaseURLToolParameter first
            if isinstance(param, BaseURLToolParameter):
                args[key] = param.get_initial_value(trans, None)
            elif isinstance(param, HiddenToolParameter):
                args[key] = model.User.expand_user_properties(trans.user, param.value)
            else:
                args[key] = param.get_initial_value(trans, None)
        return args

    def execute(self, trans, incoming=None, set_output_hid=True, history=None, **kwargs):
        """
        Execute the tool using parameter values in `incoming`. This just
        dispatches to the `ToolAction` instance specified by
        `self.tool_action`. In general this will create a `Job` that
        when run will build the tool's outputs, e.g. `DefaultToolAction`.
        """
        if incoming is None:
            incoming = {}
        try:
            return self.tool_action.execute(self, trans, incoming=incoming, set_output_hid=set_output_hid, history=history, **kwargs)
        except exceptions.ToolExecutionError as exc:
            job = exc.job
            job_id = 'unknown'
            if job is not None:
                job.mark_failed(info=exc.err_msg, blurb=exc.err_code.default_error_message)
                job_id = job.id
            log.error("Tool execution failed for job: %s", job_id)
            raise

    def params_to_strings(self, params, app, nested=False):
        return params_to_strings(self.inputs, params, app, nested)

    def params_from_strings(self, params, app, ignore_errors=False):
        return params_from_strings(self.inputs, params, app, ignore_errors)

    def check_and_update_param_values(self, values, trans, update_values=True, workflow_building_mode=False):
        """
        Check that all parameters have values, and fill in with default
        values where necessary. This could be called after loading values
        from a database in case new parameters have been added.
        """
        messages = {}
        request_context = proxy_work_context_for_history(trans, workflow_building_mode=workflow_building_mode)

        def validate_inputs(input, value, error, parent, context, prefixed_name, prefixed_label, **kwargs):
            if not error:
                value, error = check_param(request_context, input, value, context)
            if error:
                if update_values and not hasattr(input, 'data_ref'):
                    try:
                        previous_value = value
                        value = input.get_initial_value(request_context, context)
                        if not prefixed_name.startswith('__'):
                            messages[prefixed_name] = error if previous_value == value else f'{error} Using default: \'{value}\'.'
                        parent[input.name] = value
                    except Exception:
                        messages[prefixed_name] = 'Attempt to replace invalid value for \'%s\' failed.' % (prefixed_label)
                else:
                    messages[prefixed_name] = error

        visit_input_values(self.inputs, values, validate_inputs)
        return messages

    def build_dependency_cache(self, **kwds):
        if isinstance(self.app.toolbox.dependency_manager, CachedDependencyManager):
            self.app.toolbox.dependency_manager.build_cache(
                requirements=self.requirements,
                installed_tool_dependencies=self.installed_tool_dependencies,
                tool_dir=self.tool_dir,
                job_directory=None,
                metadata=False,
                tool_instance=self,
                **kwds
            )

    def build_dependency_shell_commands(self, job_directory=None, metadata=False):
        """
        Return a list of commands to be run to populate the current environment to include this tools requirements.
        """
        return self.app.toolbox.dependency_manager.dependency_shell_commands(
            requirements=self.requirements,
            installed_tool_dependencies=self.installed_tool_dependencies,
            tool_dir=self.tool_dir,
            job_directory=job_directory,
            preserve_python_environment=self.requires_galaxy_python_environment,
            metadata=metadata,
            tool_instance=self
        )

    @property
    def installed_tool_dependencies(self):
        if self.tool_shed_repository:
            installed_tool_dependencies = self.tool_shed_repository.tool_dependencies_installed_or_in_error
        else:
            installed_tool_dependencies = None
        return installed_tool_dependencies

    @property
    def tool_requirements(self):
        """
        Return all requiremens of type package
        """
        return self.requirements.packages

    @property
    def tool_requirements_status(self):
        """
        Return a list of dictionaries for all tool dependencies with their associated status
        """
        return self._view.get_requirements_status({self.id: self.tool_requirements}, self.installed_tool_dependencies)

    @property
    def output_discover_patterns(self):
        # patterns to collect for remote job execution
        patterns = []
        for output in self.outputs.values():
            patterns.extend(output.output_discover_patterns)
        return patterns

    def build_redirect_url_params(self, param_dict):
        """
        Substitute parameter values into self.redirect_url_params
        """
        if not self.redirect_url_params:
            return
        redirect_url_params = None
        # Substituting parameter values into the url params
        redirect_url_params = fill_template(self.redirect_url_params, context=param_dict)
        # Remove newlines
        redirect_url_params = redirect_url_params.replace("\n", " ").replace("\r", " ")
        return redirect_url_params

    def parse_redirect_url(self, data, param_dict):
        """
        Parse the REDIRECT_URL tool param. Tools that send data to an external
        application via a redirect must include the following 3 tool params:

        1) REDIRECT_URL - the url to which the data is being sent

        2) DATA_URL - the url to which the receiving application will send an
           http post to retrieve the Galaxy data

        3) GALAXY_URL - the url to which the external application may post
           data as a response
        """
        redirect_url = param_dict.get('REDIRECT_URL')
        redirect_url_params = self.build_redirect_url_params(param_dict)
        # Add the parameters to the redirect url.  We're splitting the param
        # string on '**^**' because the self.parse() method replaced white
        # space with that separator.
        params = redirect_url_params.split('**^**')
        rup_dict = {}
        for param in params:
            p_list = param.split('=')
            p_name = p_list[0]
            p_val = p_list[1]
            rup_dict[p_name] = p_val
        DATA_URL = param_dict.get('DATA_URL', None)
        assert DATA_URL is not None, "DATA_URL parameter missing in tool config."
        DATA_URL += f"/{str(data.id)}/display"
        redirect_url += f"?DATA_URL={DATA_URL}"
        # Add the redirect_url_params to redirect_url
        for p_name in rup_dict:
            redirect_url += f"&{p_name}={rup_dict[p_name]}"
        # Add the current user email to redirect_url
        if data.history.user:
            USERNAME = str(data.history.user.email)
        else:
            USERNAME = 'Anonymous'
        redirect_url += f"&USERNAME={USERNAME}"
        return redirect_url

    def call_hook(self, hook_name, *args, **kwargs):
        """
        Call the custom code hook function identified by 'hook_name' if any,
        and return the results
        """
        try:
            code = self.get_hook(hook_name)
            if code:
                return code(*args, **kwargs)
        except Exception as e:
            original_message = ''
            if len(e.args):
                original_message = e.args[0]
            e.args = (f"Error in '{self.name}' hook '{hook_name}', original message: {original_message}", )
            raise

    def exec_before_job(self, app, inp_data, out_data, param_dict=None):
        pass

    def exec_after_process(self, app, inp_data, out_data, param_dict, job=None, **kwds):
        pass

    def job_failed(self, job_wrapper, message, exception=False):
        """
        Called when a job has failed
        """

    def discover_outputs(self, out_data, out_collections, tool_provided_metadata, tool_working_directory, job, input_ext, input_dbkey, inp_data=None, final_job_state='ok'):
        """
        Find any additional datasets generated by a tool and attach (for
        cases where number of outputs is not known in advance).
        """
        # given the job_execution import is the only one, probably makes sense to refactor this out
        # into job_wrapper.
        tool = self
        permission_provider = output_collect.PermissionProvider(inp_data, tool.app.security_agent, job)
        metadata_source_provider = output_collect.MetadataSourceProvider(inp_data)
        job_context = output_collect.JobContext(
            tool,
            tool_provided_metadata,
            job,
            tool_working_directory,
            permission_provider,
            metadata_source_provider,
            input_dbkey,
            object_store=tool.app.object_store,
            final_job_state=final_job_state,
            flush_per_n_datasets=tool.app.config.flush_per_n_datasets,
            max_discovered_files=tool.app.config.max_discovered_files,
        )
        collected = output_collect.collect_primary_datasets(
            job_context,
            out_data,
            input_ext,
        )
        output_collect.collect_dynamic_outputs(
            job_context,
            out_collections,
        )
        # Return value only used in unit tests. Probably should be returning number of collected
        # bytes instead?
        return collected

    def to_archive(self):
        tool = self
        tarball_files = []
        temp_files = []
        with open(os.path.abspath(tool.config_file)) as fh1:
            tool_xml = fh1.read()
        # Retrieve tool help images and rewrite the tool's xml into a temporary file with the path
        # modified to be relative to the repository root.
        image_found = False
        if tool.help is not None:
            tool_help = tool.help._source
            # Check each line of the rendered tool help for an image tag that points to a location under static/
            for help_line in tool_help.split('\n'):
                image_regex = re.compile(r'img alt="[^"]+" src="\${static_path}/([^"]+)"')
                matches = re.search(image_regex, help_line)
                if matches is not None:
                    tool_help_image = matches.group(1)
                    tarball_path = tool_help_image
                    filesystem_path = os.path.abspath(os.path.join(self.app.config.root, 'static', tool_help_image))
                    if os.path.exists(filesystem_path):
                        tarball_files.append((filesystem_path, tarball_path))
                        image_found = True
                        tool_xml = tool_xml.replace('${static_path}/%s' % tarball_path, tarball_path)
        # If one or more tool help images were found, add the modified tool XML to the tarball instead of the original.
        if image_found:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".xml", delete=False
            ) as fh2:
                new_tool_config = fh2.name
                fh2.write(tool_xml)
            tool_tup = (new_tool_config, os.path.split(tool.config_file)[-1])
            temp_files.append(new_tool_config)
        else:
            tool_tup = (os.path.abspath(tool.config_file), os.path.split(tool.config_file)[-1])
        tarball_files.append(tool_tup)
        # TODO: This feels hacky.
        tool_command = tool.command.strip().split()[0]
        tool_path = os.path.dirname(os.path.abspath(tool.config_file))
        # Add the tool XML to the tuple that will be used to populate the tarball.
        if os.path.exists(os.path.join(tool_path, tool_command)):
            tarball_files.append((os.path.join(tool_path, tool_command), tool_command))
        # Find and add macros and code files.
        for external_file in tool.get_externally_referenced_paths(os.path.abspath(tool.config_file)):
            external_file_abspath = os.path.abspath(os.path.join(tool_path, external_file))
            tarball_files.append((external_file_abspath, external_file))
        if os.path.exists(os.path.join(tool_path, "Dockerfile")):
            tarball_files.append((os.path.join(tool_path, "Dockerfile"), "Dockerfile"))
        # Find tests, and check them for test data.
        tests = tool.tests
        if tests is not None:
            for test in tests:
                # Add input file tuples to the list.
                for input in test.inputs:
                    for input_value in test.inputs[input]:
                        input_filename = str(input_value)
                        input_path = os.path.abspath(os.path.join('test-data', input_filename))
                        if os.path.exists(input_path):
                            td_tup = (input_path, os.path.join('test-data', input_filename))
                            tarball_files.append(td_tup)
                # And add output file tuples to the list.
                for _, filename, _ in test.outputs:
                    output_filepath = os.path.abspath(os.path.join('test-data', filename))
                    if os.path.exists(output_filepath):
                        td_tup = (output_filepath, os.path.join('test-data', filename))
                        tarball_files.append(td_tup)
        for param in tool.input_params:
            # Check for tool data table definitions.
            if hasattr(param, 'options'):
                if hasattr(param.options, 'tool_data_table'):
                    data_table = param.options.tool_data_table
                    if hasattr(data_table, 'filenames'):
                        data_table_definitions = []
                        for data_table_filename in data_table.filenames:
                            # FIXME: from_shed_config seems to always be False.
                            if not data_table.filenames[data_table_filename]['from_shed_config']:
                                tar_file = f"{data_table.filenames[data_table_filename]['filename']}.sample"
                                sample_file = os.path.join(data_table.filenames[data_table_filename]['tool_data_path'],
                                                           tar_file)
                                # Use the .sample file, if one exists. If not, skip this data table.
                                if os.path.exists(sample_file):
                                    tarfile_path, tarfile_name = os.path.split(tar_file)
                                    tarfile_path = os.path.join('tool-data', tarfile_name)
                                    tarball_files.append((sample_file, tarfile_path))
                                data_table_definitions.append(data_table.xml_string)
                        if len(data_table_definitions) > 0:
                            # Put the data table definition XML in a temporary file.
                            table_definition = '<?xml version="1.0" encoding="utf-8"?>\n<tables>\n    %s</tables>'
                            table_definition = table_definition % "\n".join(
                                data_table_definitions
                            )
                            with tempfile.NamedTemporaryFile(
                                mode="w", delete=False
                            ) as fh3:
                                table_conf = fh3.name
                                fh3.write(table_definition)
                            tarball_files.append((table_conf, os.path.join('tool-data', 'tool_data_table_conf.xml.sample')))
                            temp_files.append(table_conf)
        # Create the tarball.
        with tempfile.NamedTemporaryFile(suffix=".tgz", delete=False) as fh4:
            tarball_archive = fh4.name
        tarball = tarfile.open(name=tarball_archive, mode='w:gz')
        # Add the files from the previously generated list.
        for fspath, tarpath in tarball_files:
            tarball.add(fspath, arcname=tarpath)
        tarball.close()
        # Delete any temporary files that were generated.
        for temp_file in temp_files:
            os.remove(temp_file)
        return tarball_archive

    def to_dict(self, trans, link_details=False, io_details=False, tool_help=False):
        """ Returns dict of tool. """

        # Basic information
        tool_dict = super().to_dict()

        tool_dict["edam_operations"] = self.edam_operations
        tool_dict["edam_topics"] = self.edam_topics
        tool_dict["hidden"] = self.hidden
        tool_dict["is_workflow_compatible"] = self.is_workflow_compatible
        tool_dict["xrefs"] = self.xrefs

        # Fill in ToolShedRepository info
        if hasattr(self, 'tool_shed') and self.tool_shed:
            tool_dict['tool_shed_repository'] = {
                'name': self.repository_name,
                'owner': self.repository_owner,
                'changeset_revision': self.changeset_revision,
                'tool_shed': self.tool_shed
            }

        # If an admin user, expose the path to the actual tool config XML file.
        if trans.user_is_admin:
            config_file = None if not self.config_file else os.path.abspath(self.config_file)
            tool_dict['config_file'] = config_file

        # Add link details.
        if link_details:
            # Add details for creating a hyperlink to the tool.
            if not isinstance(self, DataSourceTool):
                link = self.app.url_for(controller='tool_runner', tool_id=self.id)
            else:
                link = self.app.url_for(controller='tool_runner', action='data_source_redirect', tool_id=self.id)

            # Basic information
            tool_dict.update({'link': link,
                              'min_width': self.uihints.get('minwidth', -1),
                              'target': self.target})

        # Add input and output details.
        if io_details:
            tool_dict['inputs'] = [input.to_dict(trans) for input in self.inputs.values()]
            tool_dict['outputs'] = [output.to_dict(app=self.app) for output in self.outputs.values()]

        tool_dict['panel_section_id'], tool_dict['panel_section_name'] = self.get_panel_section()

        tool_class = self.__class__
        # FIXME: the Tool class should declare directly, instead of ad hoc inspection
        regular_form = tool_class == Tool or isinstance(self, (DatabaseOperationTool, InteractiveTool))
        tool_dict["form_style"] = "regular" if regular_form else "special"
        if tool_help:
            # create tool help
            help_txt = ''
            if self.help:
                help_txt = self.help.render(static_path=self.app.url_for('/static'), host_url=self.app.url_for('/', qualified=True))
                help_txt = unicodify(help_txt)
            tool_dict['help'] = help_txt

        return tool_dict

    def to_json(self, trans, kwd=None, job=None, workflow_building_mode=False, history=None):
        """
        Recursively creates a tool dictionary containing repeats, dynamic options and updated states.
        """
        if kwd is None:
            kwd = {}
        if workflow_building_mode is workflow_building_modes.USE_HISTORY or workflow_building_mode is workflow_building_modes.DISABLED:
            # We don't need a history when exporting a workflow for the workflow editor or when downloading a workflow
            history = history or trans.get_history()
            if history is None and job is not None:
                history = self.history_manager.get_owned(job.history.id, trans.user, current_history=trans.history)
            if history is None:
                raise exceptions.MessageException('History unavailable. Please specify a valid history id')

        # build request context
        request_context = proxy_work_context_for_history(trans, history, workflow_building_mode=workflow_building_mode)

        # load job parameters into incoming
        tool_message = ''
        tool_warnings = ''
        if job:
            try:
                job_params = job.get_param_values(self.app, ignore_errors=True)
                tool_warnings = self.check_and_update_param_values(job_params, request_context, update_values=True)
                self._map_source_to_history(request_context, self.inputs, job_params)
                tool_message = self._compare_tool_version(job)
                params_to_incoming(kwd, self.inputs, job_params, self.app)
            except Exception as e:
                raise exceptions.MessageException(unicodify(e))

        # create parameter object
        params = Params(kwd, sanitize=False)

        # expand incoming parameters (parameters might trigger multiple tool executions,
        # here we select the first execution only in order to resolve dynamic parameters)
        expanded_incomings, _ = expand_meta_parameters(trans, self, params.__dict__)
        if expanded_incomings:
            params.__dict__ = expanded_incomings[0]

        # do param translation here, used by datasource tools
        if self.input_translator:
            self.input_translator.translate(params)

        set_dataset_matcher_factory(request_context, self)
        # create tool state
        state_inputs: Dict[str, str] = {}
        state_errors: Dict[str, str] = {}
        populate_state(request_context, self.inputs, params.__dict__, state_inputs, state_errors)

        # create tool model
        tool_model = self.to_dict(request_context)
        tool_model['inputs'] = []
        self.populate_model(request_context, self.inputs, state_inputs, tool_model['inputs'])
        unset_dataset_matcher_factory(request_context)

        # create tool help
        tool_help = ''
        if self.help:
            tool_help = self.help.render(static_path=self.app.url_for('/static'), host_url=self.app.url_for('/', qualified=True))
            tool_help = unicodify(tool_help, 'utf-8')

        if isinstance(self.action, tuple):
            action = self.action[0] + self.app.url_for(self.action[1])
        else:
            action = self.app.url_for(self.action)

        # update tool model
        tool_model.update({
            'id': self.id,
            'help': tool_help,
            'citations': bool(self.citations),
            'sharable_url': self.sharable_url,
            'message': tool_message,
            'warnings': tool_warnings,
            'versions': self.tool_versions,
            'requirements': [{'name': r.name, 'version': r.version} for r in self.requirements],
            'errors': state_errors,
            'tool_errors': self.tool_errors,
            'state_inputs': params_to_strings(self.inputs, state_inputs, self.app, use_security=True, nested=True),
            'job_id': trans.security.encode_id(job.id) if job else None,
            'job_remap': job.remappable() if job else None,
            'history_id': trans.security.encode_id(history.id) if history else None,
            'display': self.display_interface,
            'action': action,
            'license': self.license,
            'creator': self.creator,
            'method': self.method,
            'enctype': self.enctype
        })
        return tool_model

    def populate_model(self, request_context, inputs, state_inputs, group_inputs, other_values=None):
        """
        Populates the tool model consumed by the client form builder.
        """
        other_values = ExpressionContext(state_inputs, other_values)
        for input_index, input in enumerate(inputs.values()):
            tool_dict = None
            group_state = state_inputs.get(input.name, {})
            if input.type == 'repeat':
                tool_dict = input.to_dict(request_context)
                group_size = len(group_state)
                tool_dict["cache"] = [None] * group_size
                group_cache: List[List[str]] = tool_dict["cache"]
                for i in range(group_size):
                    group_cache[i] = []
                    self.populate_model(request_context, input.inputs, group_state[i], group_cache[i], other_values)
            elif input.type == 'conditional':
                tool_dict = input.to_dict(request_context)
                if 'test_param' in tool_dict:
                    test_param = tool_dict['test_param']
                    test_param['value'] = input.test_param.value_to_basic(group_state.get(test_param['name'], input.test_param.get_initial_value(request_context, other_values)), self.app)
                    test_param['text_value'] = input.test_param.value_to_display_text(test_param['value'])
                    for i in range(len(tool_dict['cases'])):
                        current_state = {}
                        if i == group_state.get('__current_case__'):
                            current_state = group_state
                        self.populate_model(request_context, input.cases[i].inputs, current_state, tool_dict['cases'][i]['inputs'], other_values)
            elif input.type == 'section':
                tool_dict = input.to_dict(request_context)
                self.populate_model(request_context, input.inputs, group_state, tool_dict['inputs'], other_values)
            else:
                try:
                    initial_value = input.get_initial_value(request_context, other_values)
                    tool_dict = input.to_dict(request_context, other_values=other_values)
                    tool_dict['value'] = input.value_to_basic(state_inputs.get(input.name, initial_value), self.app, use_security=True)
                    tool_dict['default_value'] = input.value_to_basic(initial_value, self.app, use_security=True)
                    tool_dict['text_value'] = input.value_to_display_text(tool_dict['value'])
                except ImplicitConversionRequired:
                    tool_dict = input.to_dict(request_context)
                    # This hack leads client to display a text field
                    tool_dict['textable'] = True
                except Exception:
                    tool_dict = input.to_dict(request_context)
                    log.exception("tools::to_json() - Skipping parameter expansion '%s'", input.name)
            if input_index >= len(group_inputs):
                group_inputs.append(tool_dict)
            else:
                group_inputs[input_index] = tool_dict

    def _map_source_to_history(self, trans, tool_inputs, params):
        # Need to remap dataset parameters. Job parameters point to original
        # dataset used; parameter should be the analygous dataset in the
        # current history.
        history = trans.history

        # Create index for hdas.
        hda_source_dict = {}
        for hda in history.datasets:
            key = f'{hda.hid}_{hda.dataset.id}'
            hda_source_dict[hda.dataset.id] = hda_source_dict[key] = hda

        # Ditto for dataset collections.
        hdca_source_dict = {}
        for hdca in history.dataset_collections:
            key = f'{hdca.hid}_{hdca.collection.id}'
            hdca_source_dict[hdca.collection.id] = hdca_source_dict[key] = hdca

        # Map dataset or collection to current history
        def map_to_history(value):
            id = None
            source = None
            if isinstance(value, self.app.model.HistoryDatasetAssociation):
                id = value.dataset.id
                source = hda_source_dict
            elif isinstance(value, self.app.model.HistoryDatasetCollectionAssociation):
                id = value.collection.id
                source = hdca_source_dict
            else:
                return None
            key = f'{value.hid}_{id}'
            if key in source:
                return source[key]
            elif id in source:
                return source[id]
            else:
                return None

        def mapping_callback(input, value, **kwargs):
            if isinstance(input, DataToolParameter):
                if isinstance(value, list):
                    values = []
                    for val in value:
                        new_val = map_to_history(val)
                        if new_val:
                            values.append(new_val)
                        else:
                            values.append(val)
                    return values
                else:
                    return map_to_history(value)
            elif isinstance(input, DataCollectionToolParameter):
                return map_to_history(value)
        visit_input_values(tool_inputs, params, mapping_callback)

    def _compare_tool_version(self, job):
        """
        Compares a tool version with the tool version from a job (from ToolRunner).
        """
        tool_id = job.tool_id
        tool_version = job.tool_version
        message = ''
        try:
            select_field, tools, tool = self.app.toolbox.get_tool_components(tool_id, tool_version=tool_version, get_loaded_tools_by_lineage=False, set_selected=True)
            if tool is None:
                raise exceptions.MessageException('This dataset was created by an obsolete tool (%s). Can\'t re-run.' % tool_id)
            if (self.id != tool_id and self.old_id != tool_id) or self.version != tool_version:
                if self.id == tool_id:
                    if tool_version:
                        message = f'This job was run with tool version "{tool_version}", which is not available. '
                        if len(tools) > 1:
                            message += 'You can re-run the job with the selected tool or choose another version of the tool. '
                        else:
                            message += 'You can re-run the job with this tool version, which is a different version of the original tool. '
                else:
                    new_tool_shed_url = f'{tool.sharable_url}/{tool.changeset_revision}/'
                    old_tool_shed_url = get_tool_shed_url_from_tool_shed_registry(self.app, tool_id.split('/repos/')[0])
                    old_tool_shed_url = f'{old_tool_shed_url}/view/{tool.repository_owner}/{tool.repository_name}/'
                    message = f'This job was run with <a href=\"{old_tool_shed_url}\" target=\"_blank\">tool id \"{tool_id}\"</a>, version "{tool_version}", which is not available. '
                    if len(tools) > 1:
                        message += f'You can re-run the job with the selected <a href=\"{new_tool_shed_url}\" target=\"_blank\">tool id \"{self.id}\"</a> or choose another derivation of the tool. '
                    else:
                        message += f'You can re-run the job with <a href=\"{new_tool_shed_url}\" target=\"_blank\">tool id \"{self.id}\"</a>, which is a derivation of the original tool. '
            if not self.is_latest_version:
                message += 'There is a newer version of this tool available.'
        except Exception as e:
            raise exceptions.MessageException(unicodify(e))
        return message

    def get_default_history_by_trans(self, trans, create=False):
        return trans.get_history(create=create)

    @classmethod
    def get_externally_referenced_paths(self, path):
        """ Return relative paths to externally referenced files by the tool
        described by file at `path`. External components should not assume things
        about the structure of tool xml files (this is the tool's responsibility).
        """
        tree = raw_tool_xml_tree(path)
        root = tree.getroot()
        external_paths = []
        for code_elem in root.findall('code'):
            external_path = code_elem.get('file')
            if external_path:
                external_paths.append(external_path)
        external_paths.extend(imported_macro_paths(root))
        # May also need to load external citation files as well at some point.
        return external_paths

