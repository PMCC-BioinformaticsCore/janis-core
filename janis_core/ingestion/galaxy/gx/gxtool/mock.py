
import os
import shutil
import tempfile
import yaml
from typing import Any

from janis_core.ingestion.galaxy import runtime
from janis_core.ingestion.galaxy.runtime.paths import GALAXY_CONFIG

from sqlalchemy.orm.scoping import scoped_session

from galaxy import (
    di,
    model,
    objectstore,
    quota,
)
from galaxy.tools.parameters.basic import (
    IntegerToolParameter,
)
from galaxy.auth import AuthManager
from galaxy.datatypes.registry import Registry
from galaxy.jobs.manager import NoopManager
from galaxy.jobs import SharedComputeEnvironment
from galaxy.managers.users import UserManager
from galaxy.model import mapping, tags
from galaxy.model.base import SharedModelMapping
from galaxy.model.mapping import GalaxyModelMapping
from galaxy.security import idencoding
from galaxy.structured_app import BasicApp, MinimalManagerApp, StructuredApp
from galaxy.tool_util.deps.containers import NullContainerFinder
from galaxy.util import StructuredExecutionTimer
from galaxy.util import galaxy_directory
from galaxy.util.bunch import Bunch
from galaxy.util.dbkeys import GenomeBuilds
from galaxy.web_stack import ApplicationStack
from galaxy.util import XML
from galaxy.tools.parameters import params_from_strings # for MockTool
from galaxy.tool_util.parser.output_objects import ToolOutput
from galaxy.tool_util.biotools import BiotoolsMetadataSource
from galaxy.tools.data import ToolDataTableManager
#from galaxy.datatypes.registry import example_datatype_registry_for_sample
# terrible stuff from galaxy. not planning on a refactor.


class MockApp(di.Container):
    def __init__(self, config=None, **kwargs):
        super().__init__()
        self[BasicApp] = self
        self[MinimalManagerApp] = self
        self[StructuredApp] = self
        self.config = config or MockAppConfig(**kwargs)
        self.security = self.config.security
        self[idencoding.IdEncodingHelper] = self.security
        self.name = kwargs.get('name', 'galaxy')
        self.object_store = objectstore.build_object_store_from_config(self.config)
        self.model = mapping.init("/tmp", self.config.database_connection, create_tables=True, object_store=self.object_store)
        self[SharedModelMapping] = self.model
        self[GalaxyModelMapping] = self.model
        self[scoped_session] = self.model.context
        self.security_agent = self.model.security_agent
        self.visualizations_registry = MockVisualizationsRegistry()
        self.tag_handler = tags.GalaxyTagHandler(self.model.context)
        self[tags.GalaxyTagHandler] = self.tag_handler
        self.quota_agent = quota.DatabaseQuotaAgent(self.model)
        self.init_datatypes()
        # self.job_config = Bunch(
        #     dynamic_params=None,
        #     destinations={},
        #     use_messaging=False,
        #     assign_handler=lambda *args, **kwargs: None
        # )
        self.tool_data_tables = self.grace_init_data_tables()
        self.dataset_collections_service = None
        self.container_finder = NullContainerFinder()
        self._toolbox_lock = MockLock()
        self.tool_shed_registry = Bunch(tool_sheds={})
        self.genome_builds = GenomeBuilds(self)
        self.job_manager = NoopManager()
        self.application_stack = ApplicationStack()
        self.auth_manager = AuthManager(self.config)
        self.user_manager = UserManager(self)
        self.execution_timer_factory = Bunch(get_timer=StructuredExecutionTimer)
        self.is_job_handler = False
        self.biotools_metadata_source = BiotoolsMetadataSource()
        # new param not native to the aApp class
        self.dataset_counter: int = 1

        def url_for(*args, **kwds):
            return "/mock/url"
        self.url_for = url_for

    def grace_init_data_tables(self) -> ToolDataTableManager | dict[str, str]:
        # Initialize tool data tables using the config defined by self.config.tool_data_table_config_path.
        if self.config.tool_data_path and self.config.tool_data_table_config_path:
            return ToolDataTableManager(
                tool_data_path=self.config.tool_data_path,
                config_filename=self.config.tool_data_table_config_path,
                other_config_dict=self.config
            )
        return {}

    def init_datatypes(self):
        # config = Bunch(sniff_compressed_dynamic_datatypes_default=True)
        # datatypes_registry = registry.Registry(config=config)
        # datatypes_registry.load_datatypes(root_dir='galaxy', config='config/datatypes_conf.xml.sample')
        # datatypes_registry = example_datatype_registry_for_sample()
        # model.set_datatypes_registry(datatypes_registry)
        # self.datatypes_registry = datatypes_registry
        galaxy_dir = galaxy_directory()
        sample_conf = os.path.join(galaxy_dir, "site-packages", "galaxy", "config", "sample", "datatypes_conf.xml.sample")
        config = Bunch(sniff_compressed_dynamic_datatypes_default=True)
        datatypes_registry = Registry(config)
        datatypes_registry.load_datatypes(root_dir=galaxy_dir, config=sample_conf)
        model.set_datatypes_registry(datatypes_registry)
        self.datatypes_registry = datatypes_registry

    def wait_for_toolbox_reload(self, toolbox):
        return True


def grace_get_config() -> dict[str, Any]:
    config_path = GALAXY_CONFIG
    with open(config_path, "r") as fp:
        return yaml.safe_load(fp)

class MockAppConfig(Bunch):

    class MockSchema(Bunch):
        pass

    def __init__(self, root=None, **kwargs):
        Bunch.__init__(self, **kwargs)
        if not root:
            root = tempfile.mkdtemp()
            self._remove_root = True
        else:
            self._remove_root = False   
        self.schema = self.MockSchema()
        self.security = idencoding.IdEncodingHelper(id_secret='6e46ed6483a833c100e68cc3f1d0dd76')
        self.database_connection = kwargs.get('database_connection', "sqlite:///:memory:")
        self.use_remote_user = kwargs.get('use_remote_user', False)
        self.enable_celery_tasks = False
        self.data_dir = os.path.join(root, 'database')
        self.file_path = os.path.join(self.data_dir, 'files')
        self.jobs_directory = os.path.join(self.data_dir, 'jobs_directory')
        self.new_file_path = os.path.join(self.data_dir, 'tmp')

        # GRACE TOOL DATA
        self.tool_data_path = None
        self.tool_data_table_config_path = None
        self.tool_dependency_dir = None
        self.set_tool_data_attrs()

        # GRACE LOGGING
        self.logging = grace_get_config()

        self.metadata_strategy = 'legacy'

        self.object_store_config_file = ''
        self.object_store = 'disk'
        self.object_store_check_old_style = False
        self.object_store_cache_path = '/tmp/cache'
        self.umask = os.umask(0o77)
        self.gid = os.getgid()

        self.user_activation_on = False
        self.new_user_dataset_access_role_default_private = False

        self.expose_dataset_path = True
        self.allow_user_dataset_purge = True
        self.enable_old_display_applications = True
        self.redact_username_in_logs = False
        self.auth_config_file = "config/auth_conf.xml.sample"
        self.error_email_to = "admin@email.to"
        self.password_expiration_period = 0

        self.umask = 0o77
        self.flush_per_n_datasets = 0

        # Compliance related config
        self.redact_email_in_job_name = False

        # Follow two required by GenomeBuilds
        self.len_file_path = os.path.join('tool-data', 'shared', 'ucsc', 'chrom')
        self.builds_file_path = os.path.join('tool-data', 'shared', 'ucsc', 'builds.txt.sample')

        self.shed_tool_config_file = "config/shed_tool_conf.xml"
        self.shed_tool_config_file_set = False
        self.enable_beta_edam_toolbox = False
        self.preserve_python_environment = "always"
        self.enable_beta_gdpr = False

        self.version_major = "19.09"

        # set by MockDir
        self.root = root
        self.enable_tool_document_cache = False
        self.tool_cache_data_dir = os.path.join(root, 'tool_cache')
        self.delay_tool_initialization = True
        self.external_chown_script = None

        self.default_panel_view = "default"
        self.panel_views_dir = ''
        self.panel_views = {}
        self.edam_panel_views = ''

        self.config_file = None

    def set_tool_data_attrs(self) -> None:
        tool_data = f'{runtime.tool.xml_dir()}/tool-data'
        tool_data_table_conf = f'{runtime.tool.xml_dir()}/tool_data_table_conf.xml.sample'
        if os.path.exists(tool_data):
            self.tool_data_path = tool_data
        if os.path.exists(tool_data_table_conf):
            self.tool_data_table_config_path = tool_data_table_conf

    @property
    def config_dict(self):
        return self.dict()

    def __getattr__(self, name):
        # Handle the automatic [option]_set options: for tests, assume none are set
        if name == 'is_set':
            return lambda x: False
        # Handle the automatic config file _set options
        if name.endswith('_file_set'):
            return False
        raise AttributeError(name)

    def __del__(self):
        if self._remove_root:
            shutil.rmtree(self.root)



class MockVisualizationsRegistry:
    BUILT_IN_VISUALIZATIONS = ['trackster']

    def get_visualizations(self, trans, target):
        return []



class MockLock:
    def __enter__(self):
        pass

    def __exit__(self, type, value, traceback):
        pass



class MockObjectStore:
    def __init__(self):
        self.created_datasets = []
        self.first_create = True
        self.object_store_id = "mycoolid"

    def exists(self, *args, **kwargs):
        return True

    def create(self, dataset):
        self.created_datasets.append(dataset)
        if self.first_create:
            self.first_create = False
            assert dataset.object_store_id is None
            dataset.object_store_id = self.object_store_id
        else:
            assert dataset.object_store_id == self.object_store_id


class ComputeEnvironment(SharedComputeEnvironment):
    def __init__(
        self,
        new_file_path,
        working_directory
    ):
        self._tool_dir = ''
        self._new_file_path = new_file_path
        self._working_directory = working_directory
        self._unstructured_path_rewrites = {}

    def input_path_rewrite(self, dataset):
        return dataset.file_name.rsplit('/', 1)[-1]

    def output_path_rewrite(self, dataset):
        return dataset.file_name.rsplit('/', 1)[-1]

    def output_paths(self):
        return self._output_paths

    def working_directory(self):
        return self._working_directory

    def home_directory(self):
        return self._working_directory

    def tmp_directory(self):
        return self._working_directory

    def new_file_path(self):
        return self._new_file_path

    def unstructured_path_rewrite(self, path):
        for key, val in self._unstructured_path_rewrites.items():
            if path.startswith(key):
                return path.replace(key, val)
        return None

    def tool_directory(self):
        return self._tool_dir

    def galaxy_url(self):
        return 'http://localhost:9090/'


# not used?
class MockTool:
    def __init__(self, app: MockApp):
        self.profile = 16.01
        self.python_template_version = '2.7'
        self.app = app
        self.hooks_called = []
        self.environment_variables = []
        self._config_files = []
        self._command_line = "bwa --thresh=$thresh --in1=$input1 --in2=$input2 --out=$output1"
        self._params = {"thresh": self.test_thresh_param()}
        self.options = Bunch(sanitize=False)
        self.check_values = True

    def test_thresh_param(self):
        elem = XML('<param name="thresh" type="integer" value="5" />')
        return IntegerToolParameter(self, elem)

    def params_from_strings(self, params, app, ignore_errors=False):
        return params_from_strings(self.inputs, params, app, ignore_errors)

    @property
    def config_file(self):
        return "<fake tool>"

    @property
    def template_macro_params(self):
        return {}

    @property
    def inputs(self):
        return self._params

    def set_params(self, params):
        self._params = params

    @property
    def outputs(self):
        return dict(
            output1=ToolOutput("output1"),
        )

    @property
    def tmp_directory_vars(self):
        return ["TMP"]

    @property
    def config_files(self):
        return self._config_files

    @property
    def command(self):
        return self._command_line

    @property
    def interpreter(self):
        return None

    def handle_unvalidated_param_values(self, input_values, app):
        pass

    def build_param_dict(self, incoming, *args, **kwds):
        return incoming

    def call_hook(self, hook_name, *args, **kwargs):
        self.hooks_called.append(hook_name)

    def exec_before_job(self, *args, **kwargs):
        pass