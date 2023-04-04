from setuptools import setup, find_packages

# Version information is found in the __init__ file of `janiscore/`
DESCRIPTION = "Contains classes and helpers to build a workflow, and provide options to convert to CWL / WDL"

######## SHOULDN'T NEED EDITS BELOW THIS LINE ########

with open("./README.md") as readme:
    long_description = readme.read()


#   Wondering why we have to get the version like this,
#   instead of simply from janis import __version__?
#
#   Well, well well, let me tell you...
#
#   This cost me days, but basically in `pip v19.0.0`,
#   they changed the way that modules are installed with PEP517.
#
#   Module installation is more isolated, and as this setup.py file called
#   for janis.__version__ AND we use a pyproject.toml file, pip / setuptools
#   gets confused and "From a PEP 517 perspective ..., it means it needs itself
#   to build. Which is of course a bit weird." [1] This means it is essentially
#   tries to import the project without resolving any of the setup.install_requires
#   dependencies. Installing the project with the `--no-use-pep517` flag works as expected.
#
#   I've elected to go with option (3) from [2]:
#
#   Links:
#       [1] https://github.com/pypa/pip/issues/6163#issuecomment-456738585
#       [2] https://packaging.python.org/guides/single-sourcing-package-version/
#       [3] https://github.com/pyinstaller/pyinstaller/issues/2730
#       [4] https://github.com/pypa/pip/issues/6163
#
vsn = {}
with open("./janis_core/__meta__.py") as fp:
    exec(fp.read(), vsn)
__version__ = vsn["__version__"]
githuburl = vsn["GITHUB_URL"]


setup(
    name="janis-pipelines.core",
    version=__version__,
    description=DESCRIPTION,
    url=githuburl,
    author="Michael Franklin, Evan Thomas, Mohammad Bhuyan",
    author_email="michael.franklin@petermac.org",
    license="GNU",
    keywords=["janis", "workflows", "pipelines", "cwl", "wdl", "bioinformatics"],
    packages=["janis_core"]
    + ["janis_core." + p for p in sorted(find_packages("./janis_core"))],
    install_requires=[
        "importlib-metadata",
        "illusional.wdlgen >= 0.3.0",
        "ruamel.yaml >= 0.12.4, <= 0.16.5",
        "tabulate",
        "path",
        "cwlformat == 2020.5.19",
        "cwl-utils == 0.11",
        "graphviz",
        "nose",
        "miniwdl",
    ],
    # entry_points={"janis.extension": ["core=core"]},
    zip_safe=False,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
