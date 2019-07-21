from core import Workflow, Input, File, Output, Step
from core.unix.data_types.tarfile import TarFile
from core.unix.tools.compile import Compile
from core.unix.tools.untar import Untar

w = Workflow("user_guide")

tarball = Input("tarball", TarFile(), value="path/to/tar")
compiled_class = Output("compiledOutput")

untar = Step("untar", Untar())
compile_step = Step("compile", Compile())

w.add_edges(
    [
        (tarball, untar.tarFile),
        (untar.files, compile_step.file),
        (compile_step.compiled, compiled_class),
    ]
)

w.translate("cwl", to_disk=False)  # or "wdl"
