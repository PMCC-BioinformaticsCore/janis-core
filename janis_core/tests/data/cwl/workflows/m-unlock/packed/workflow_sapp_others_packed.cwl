{
    "$graph": [
        {
            "class": "CommandLineTool",
            "label": "Perform a grep count on a file",
            "hints": [
                {
                    "dockerPull": "debian:buster",
                    "class": "DockerRequirement"
                }
            ],
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "ShellCommandRequirement"
                }
            ],
            "baseCommand": [],
            "stdout": "count.txt",
            "inputs": [
                {
                    "type": "string",
                    "default": "grep",
                    "inputBinding": {
                        "position": 2,
                        "shellQuote": false
                    },
                    "id": "#grep_count.cwl/binary"
                },
                {
                    "type": "string",
                    "default": "-c",
                    "inputBinding": {
                        "position": 3,
                        "shellQuote": false
                    },
                    "id": "#grep_count.cwl/count"
                },
                {
                    "type": "string",
                    "inputBinding": {
                        "position": 4
                    },
                    "id": "#grep_count.cwl/grep"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "shellQuote": false,
                        "position": 5
                    },
                    "id": "#grep_count.cwl/infile"
                },
                {
                    "type": "string",
                    "default": "||",
                    "inputBinding": {
                        "shellQuote": false,
                        "position": 6
                    },
                    "id": "#grep_count.cwl/pipe"
                },
                {
                    "type": "string",
                    "default": "true",
                    "inputBinding": {
                        "shellQuote": false,
                        "position": 7
                    },
                    "id": "#grep_count.cwl/status"
                }
            ],
            "id": "#grep_count.cwl",
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
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential",
            "outputs": [
                {
                    "type": "int",
                    "outputBinding": {
                        "glob": "count.txt",
                        "loadContents": true,
                        "outputEval": "$(parseInt(self[0].contents))"
                    },
                    "id": "#grep_count.cwl/matches"
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
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Codon table used for this organism",
                    "label": "Codon table",
                    "inputBinding": {
                        "prefix": "-codon"
                    },
                    "id": "#conversion.cwl/codon"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Reference genome file used in EMBL format",
                    "label": "Reference genome",
                    "inputBinding": {
                        "prefix": "-input"
                    },
                    "id": "#conversion.cwl/embl"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Reference genome file used in fasta format",
                    "label": "Reference genome",
                    "inputBinding": {
                        "prefix": "-input"
                    },
                    "id": "#conversion.cwl/fasta"
                },
                {
                    "type": "string",
                    "doc": "Name of the sample being analysed",
                    "label": "Sample name",
                    "inputBinding": {
                        "prefix": "-id"
                    },
                    "id": "#conversion.cwl/identifier"
                }
            ],
            "baseCommand": [
                "java",
                "-Xmx5g",
                "-jar",
                "/SAPP-2.0.jar"
            ],
            "arguments": [
                {
                    "valueFrom": "${\n  if (inputs.embl) {\n    return '-embl2rdf';\n  }\n  if (inputs.fasta) {\n    return '-fasta2rdf';\n  }\n}\n"
                },
                {
                    "valueFrom": "${\n  if (inputs.fasta) {\n    return '-genome';\n  } \n  return null;\n}\n"
                },
                {
                    "prefix": "-output",
                    "valueFrom": "$(inputs.identifier).ttl"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier).ttl"
                    },
                    "id": "#conversion.cwl/output"
                }
            ],
            "id": "#conversion.cwl",
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
            "class": "CommandLineTool",
            "label": "InterProScan annotation",
            "doc": "Runs Genome annotation with InterProScan using an RDF genome file according to the GBOL ontology\n",
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
                },
                {
                    "ramMax": 10000,
                    "coresMax": 10,
                    "class": "ResourceRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "doc": "Which applications to use (e.g. TIGRFAM,SFLD,SUPERFAMILY,Gene3D,Hamap,Coils,ProSiteProfiles,SMART,PRINTS,ProSitePatterns,Pfam,ProDom,MobiDBLite,PIRSF)",
                    "label": "applications",
                    "inputBinding": {
                        "prefix": "-a"
                    },
                    "default": "Pfam",
                    "id": "#interproscan.cwl/applications"
                },
                {
                    "type": "int",
                    "doc": "Number of CPUs acccessible for InterProScan",
                    "label": "cpu limit",
                    "inputBinding": {
                        "prefix": "-cpu"
                    },
                    "default": 2,
                    "id": "#interproscan.cwl/cpu"
                },
                {
                    "type": "string",
                    "doc": "Name of the sample being analysed",
                    "label": "Sample name",
                    "id": "#interproscan.cwl/identifier"
                },
                {
                    "type": "File",
                    "doc": "Reference genome file used in RDF format",
                    "label": "Reference genome",
                    "inputBinding": {
                        "prefix": "-input"
                    },
                    "id": "#interproscan.cwl/input"
                },
                {
                    "type": "Directory",
                    "doc": "Path to the interproscan.sh folder",
                    "label": "InterProScan folder",
                    "inputBinding": {
                        "prefix": "-path"
                    },
                    "id": "#interproscan.cwl/interpro"
                }
            ],
            "baseCommand": [
                "java",
                "-Xmx5g",
                "-jar",
                "/SAPP-2.0.jar",
                "-interpro"
            ],
            "arguments": [
                {
                    "prefix": "-output",
                    "valueFrom": "$(inputs.identifier).hdt"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier).hdt"
                    },
                    "id": "#interproscan.cwl/output"
                }
            ],
            "id": "#interproscan.cwl",
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
            "label": "Genome conversion and annotation",
            "doc": "Workflow for genome annotation from EMBL format",
            "inputs": [
                {
                    "type": "int",
                    "doc": "Codon table used for gene prediction",
                    "label": "Codon table",
                    "default": 1,
                    "id": "#main/codon"
                },
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
                    "label": "EMBL input file",
                    "id": "#main/embl"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Genome sequence in FASTA format",
                    "label": "EMBL input file",
                    "id": "#main/fasta"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "FASTA to RDF conversion with genome data",
                    "id": "#main/genome"
                },
                {
                    "type": "string",
                    "doc": "Identifier of the sample being converted",
                    "label": "Sample name",
                    "id": "#main/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "number of threads to use for computational processes",
                    "label": "number of threads",
                    "default": 2,
                    "id": "#main/threads"
                }
            ],
            "steps": [
                {
                    "run": "#conversion.cwl",
                    "in": [
                        {
                            "source": "#main/codon",
                            "id": "#main/conversion/codon"
                        },
                        {
                            "source": "#main/embl",
                            "id": "#main/conversion/embl"
                        },
                        {
                            "source": "#main/fasta",
                            "id": "#main/conversion/fasta"
                        },
                        {
                            "source": "#main/genome",
                            "id": "#main/conversion/genome"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/conversion/identifier"
                        }
                    ],
                    "out": [
                        "#main/conversion/output"
                    ],
                    "id": "#main/conversion"
                },
                {
                    "run": "#grep_count.cwl",
                    "in": [
                        {
                            "default": "/Gene>",
                            "id": "#main/count_genes/grep"
                        },
                        {
                            "source": "#main/conversion/output",
                            "id": "#main/count_genes/infile"
                        }
                    ],
                    "out": [
                        "#main/count_genes/matches"
                    ],
                    "id": "#main/count_genes"
                },
                {
                    "when": "$(inputs.matches > 0)",
                    "run": "#interproscan.cwl",
                    "in": [
                        {
                            "source": "#main/threads",
                            "id": "#main/interproscan/cpu"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/interproscan/identifier"
                        },
                        {
                            "source": "#main/conversion/output",
                            "id": "#main/interproscan/input"
                        },
                        {
                            "source": "#main/count_genes/matches",
                            "id": "#main/interproscan/matches"
                        }
                    ],
                    "out": [
                        "#main/interproscan/output"
                    ],
                    "id": "#main/interproscan"
                },
                {
                    "run": "#toHDT.cwl",
                    "in": [
                        {
                            "source": "#main/interproscan/output",
                            "id": "#main/tohdt/input"
                        },
                        {
                            "valueFrom": "$(inputs.identifier).hdt",
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
                    "type": "File",
                    "outputSource": "#main/tohdt/output",
                    "id": "#main/output"
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
