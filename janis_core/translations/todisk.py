
import os 
import shutil
from typing import Tuple, Any,Optional

from janis_core import Logger
from janis_core import settings
from janis_core.translation_deps.exportpath import ExportPathKeywords


def write_workflow_to_console(
    tup_workflow: Tuple[str, str], 
    tup_tools: list[Tuple[str, str]], 
    tup_inp: Tuple[str, str], 
    tup_resources: Optional[Tuple[str, str]]
    ) -> None:
    """Each tuple is (filename, string). Tells us where to write each file, followed by the file contents."""

    print(tup_workflow[1])
    if settings.translate.TOOL_TO_CONSOLE:
        print("\n=== TOOLS ===")
        [print(f":: {t[0]} ::\n" + t[1]) for t in tup_tools]
    print("\n=== INPUTS ===")
    print(tup_inp)
    if tup_resources is not None:
        if not settings.translate.MERGE_RESOURCES and settings.translate.WITH_RESOURCE_OVERRIDES:
            print("\n=== RESOURCES ===")
            print(tup_resources[1])


def write_workflow_to_disk(
    tup_main: Tuple[str, str], 
    tup_subworkflows: list[Tuple[str, str]], 
    tup_tools: list[Tuple[str, str]], 
    tup_inputs: Tuple[str, str], 
    tup_helpers: list[Tuple[str, str]], 
    tup_resources: Optional[Tuple[str, str]],
    outdir_structure: dict[str, Any]
    ) -> None:
    """Each tuple is (filename, string). Tells us where to write each file, followed by the file contents."""
    
    # prepare base output directory
    basedir = ExportPathKeywords.resolve(
        settings.translate.EXPORT_PATH, 
        workflow_spec=settings.translate.DEST,
        workflow_name=None
    )
    if os.path.isdir(basedir):
        shutil.rmtree(basedir)
    os.makedirs(basedir, exist_ok=True)

    # writing main workflow
    _write_file(tup_main, basedir, 'main', outdir_structure['main'])
    
    # writing tools
    for tup_tool in tup_tools:
        _write_file(tup_tool, basedir, 'tools', outdir_structure['tools'])
    
    # writing subworkflows
    for tup_subworkflow in tup_subworkflows:
        _write_file(tup_subworkflow, basedir, 'subworkflows', outdir_structure['subworkflows'])

    # writing inputs file
    if settings.translate.WRITE_INPUTS_FILE:
        _write_file(tup_inputs, basedir, 'inputs', outdir_structure['inputs'])
    else:
        Logger.log("Skipping writing inputs config file")

    # writing resources file
    if tup_resources is not None:
        if settings.translate.WITH_RESOURCE_OVERRIDES and not settings.translate.MERGE_RESOURCES:
            _write_file(tup_resources, basedir, 'resources', outdir_structure['resources'])
        else:
            Logger.log("Skipping writing resources config file")

    # writing helper files
    for tup_helper in tup_helpers:
        _write_file(tup_helper, basedir, 'helpers', outdir_structure['helpers'])

    # copying source files - this one is a bit weird & specific to galaxy.  
    if settings.general.SOURCE_FILES is not None:
        # create source folder in basedir
        source_dir = os.path.join(basedir, 'source')
        if not os.path.isdir(source_dir):
            os.mkdir(source_dir)
        # copy files
        for src, dest in settings.general.SOURCE_FILES:
            dest = os.path.join(source_dir, dest)
            if not os.path.isdir(os.path.dirname(dest)):
                os.mkdir(os.path.dirname(dest))
            shutil.copy2(src, dest)


def _write_file(tup_file: Tuple[str, str], basedir: str, ftype: str, fsubdir: str | None) -> None:
    filename, contents = tup_file
    
    # format outdir using basedir and subdir if provided
    outdir = os.path.join(basedir, fsubdir) if fsubdir is not None else basedir

    # write to disk
    Logger.info(f"Writing {ftype} to '{outdir}'")
    with open(os.path.join(outdir, filename), "w+") as f:
        Logger.log(f"Writing {filename} to disk")
        f.write(contents)
        Logger.log(f"Written {filename} to disk")


def write_tool_to_console(str_tool: str) -> None:
    print(str_tool)


def write_tool_to_disk(str_tool: str, filename: str, helpers: dict[str, str]) -> None:
    # set output folder
    basedir = ExportPathKeywords.resolve(
        settings.translate.EXPORT_PATH, 
        workflow_spec=settings.translate.DEST, 
        workflow_name=None
    )
    # create output folder
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    # write tool file
    with open(os.path.join(basedir, filename), "w+") as wf:
        Logger.log(f"Writing {filename} to disk")
        wf.write(str_tool)
        Logger.log(f"Wrote {filename}  to disk")

    # write helper files (files_to_create scripts)
    helpers = {fn.split('/')[-1]: fc for fn, fc in helpers.items()}
    for (filename, filecontents) in helpers.items():
        with open(os.path.join(basedir, filename), "w+") as helperfp:
            Logger.log(f"Writing {filename} to disk")
            helperfp.write(filecontents)
            Logger.log(f"Written {filename} to disk")

    # writing source files to output folder (specifically galaxy tool wrappers)
    # copying source files 
    if settings.general.SOURCE_FILES is not None:
        # create source folder in basedir
        source_dir = os.path.join(basedir, 'source')
        if not os.path.isdir(source_dir):
            os.mkdir(source_dir)
        
        # copy files
        for src, dest in settings.general.SOURCE_FILES:
            dest = os.path.join(source_dir, dest)
            if not os.path.isdir(os.path.dirname(dest)):
                os.mkdir(os.path.dirname(dest))
            shutil.copy2(src, dest)

    

# DEPRECATED from write_workflow_to_disk()

    # subfolders: list[str] = []
    # subfolders.append(self.DIR_TOOLS)
    # subfolders += self.SUBDIRS_TO_CREATE
    # for subfolder in subfolders:
    #     path = os.path.join(basedir, subfolder)
    #     if not os.path.isdir(path):
    #         os.makedirs(path, exist_ok=True)


        # if settings.translate.SHOULD_VALIDATE:
    #     with Path(basedir):

    #         Logger.info(f"Validating outputted {self.name}")

    #         enved_vcs = [
    #             (os.getenv(x[1:]) if x.startswith("$") else x)
    #             for x in self.validate_command_for(
    #                 fn_workflow, fn_inputs, "tools/", "tools.zip"
    #             )
    #         ]

    #         cwltool_result = subprocess.run(enved_vcs)
    #         if cwltool_result.returncode == 0:
    #             Logger.info(
    #                 "Exported tool was validated by: " + " ".join(enved_vcs)
    #             )
    #         else:
    #             Logger.critical(str(cwltool_result.stderr))
    
    # # zipping tools file
    # import subprocess

    # if settings.translate.SHOULD_ZIP:
    #     Logger.debug("Zipping tools")
    #     with Path(basedir):
    #         FNULL = open(os.devnull, "w")
    #         zip_result = subprocess.run(
    #             ["zip", "-r", "tools.zip", "tools/"], stdout=FNULL
    #         )
    #         if zip_result.returncode == 0:
    #             Logger.debug("Zipped tools")
    #         else:
    #             Logger.critical(str(zip_result.stderr.decode()))