{
    "$graph": [
        {
            "class": "CommandLineTool",
            "label": "compress any file with gzip",
            "hints": [
                {
                    "dockerPull": "debian:buster",
                    "class": "DockerRequirement"
                }
            ],
            "baseCommand": [
                "gzip",
                "-c"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.infile)"
                }
            ],
            "stdout": "$(inputs.infile.basename).gz",
            "inputs": [
                {
                    "type": "File",
                    "id": "#compress.cwl/infile"
                }
            ],
            "id": "#compress.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2021-00-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.infile.basename).gz"
                    },
                    "id": "#compress.cwl/output"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "label": "Genome conversion",
            "doc": "Runs Genome conversion tool from SAPP\n",
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/sapp:2.0",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "17.0.3"
                            ],
                            "specs": [
                                "https://anaconda.org/conda-forge/openjdk"
                            ],
                            "package": "sapp"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "RDF file in any format with the right extension",
                    "label": "RDF input file (hdt, ttl, nt, etc...)",
                    "inputBinding": {
                        "prefix": "-i"
                    },
                    "id": "#toHDT.cwl/input"
                },
                {
                    "type": "string",
                    "doc": "Name of the output file with the right extension",
                    "label": "RDF output file (hdt, ttl, nt, etc...)",
                    "inputBinding": {
                        "prefix": "-o"
                    },
                    "id": "#toHDT.cwl/output"
                }
            ],
            "baseCommand": [
                "java",
                "-Xmx5g",
                "-jar",
                "/SAPP-2.0.jar",
                "-convert"
            ],
            "outputs": [
                {
                    "$import": "#toHDT.cwl/output"
                }
            ],
            "id": "#toHDT.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2020-00-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "Workflow",
            "requirements": [
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                }
            ],
            "label": "Conversion and compression of RDF files",
            "doc": "Workflow to convert a RDF file to the HDT format and GZIP compress it for long term storage",
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output Destination",
                    "doc": "Optional Output destination used for cwl-prov reporting.",
                    "id": "#main/destination"
                },
                {
                    "type": "File",
                    "doc": "Genome sequence in EMBL format",
                    "label": "RDF input file",
                    "id": "#main/input"
                },
                {
                    "type": "string",
                    "doc": "Output filename",
                    "label": "HDT filename as output",
                    "id": "#main/output"
                }
            ],
            "steps": [
                {
                    "run": "#compress.cwl",
                    "in": [
                        {
                            "source": "#main/tohdt/output",
                            "id": "#main/compress/infile"
                        }
                    ],
                    "out": [
                        "#main/compress/output"
                    ],
                    "id": "#main/compress"
                },
                {
                    "run": "#toHDT.cwl",
                    "in": [
                        {
                            "source": "#main/input",
                            "id": "#main/tohdt/input"
                        },
                        {
                            "source": "#main/output",
                            "id": "#main/tohdt/output"
                        }
                    ],
                    "out": [
                        "#main/tohdt/output"
                    ],
                    "id": "#main/tohdt"
                }
            ],
            "outputs": [
                {
                    "$import": "#main/output"
                }
            ],
            "id": "#main",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "https://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "https://schema.org/name": "Bart Nijsse"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2020-00-00",
            "https://schema.org/dateModified": "2022-05-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        }
    ],
    "cwlVersion": "v1.2",
    "$namespaces": {
        "s": "https://schema.org/"
    }
}
