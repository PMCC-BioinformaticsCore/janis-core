{
    "$graph": [
        {
            "class": "ExpressionTool",
            "doc": "Transforms the input files to a mentioned directory\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "id": "#files_to_folder.cwl/destination"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "id": "#files_to_folder.cwl/files"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "Directory"
                        }
                    ],
                    "id": "#files_to_folder.cwl/folders"
                }
            ],
            "expression": "${\n  var array = []\n  if (inputs.files != null) {\n    array = array.concat(inputs.files)\n  }\n  if (inputs.folders != null) {\n    array = array.concat(inputs.folders)\n  }\n  var r = {\n     'results':\n       { \"class\": \"Directory\",\n         \"basename\": inputs.destination,\n         \"listing\": array\n       } \n     };\n   return r; \n }\n",
            "id": "#files_to_folder.cwl",
            "http://schema.org/citation": "https://m-unlock.nl",
            "http://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "http://schema.org/dateCreated": "2020-00-00",
            "http://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "http://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential",
            "outputs": [
                {
                    "type": "Directory",
                    "id": "#files_to_folder.cwl/results"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "label": "NGTax amplicon analysis",
            "doc": "Runs NGTAX amplicon analysis\n",
            "requirements": [
                {
                    "listing": [
                        "$(inputs.reference_db)"
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/ngtax:2.2.9",
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
                            "package": "java"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "doc": "Folder containing demultiplexed files_to_folder",
                    "label": "Demultiplexed folder",
                    "inputBinding": {
                        "prefix": "-folder"
                    },
                    "id": "#ngtax.cwl/folder"
                },
                {
                    "type": "int",
                    "doc": "Read length of the reverse read",
                    "label": "Reverse read length",
                    "inputBinding": {
                        "prefix": "-for_read_len"
                    },
                    "id": "#ngtax.cwl/for_read_len"
                },
                {
                    "type": "string",
                    "doc": "Forward primer used",
                    "label": "The forward primer used",
                    "inputBinding": {
                        "prefix": "-for_p"
                    },
                    "id": "#ngtax.cwl/forward_primer"
                },
                {
                    "type": "string",
                    "doc": "Subfragment that is being analysed (e.g. V1-V3 or V5-region)",
                    "label": "Subfragment name",
                    "id": "#ngtax.cwl/fragment"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Mock3 reference selection",
                    "label": "Mock3 reference",
                    "id": "#ngtax.cwl/mock3"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Mock4 reference selection",
                    "label": "Mock4 reference",
                    "id": "#ngtax.cwl/mock4"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Wether the primers are removed or not from the input files",
                    "label": "Primers are removed",
                    "inputBinding": {
                        "prefix": "-primersRemoved"
                    },
                    "id": "#ngtax.cwl/primersRemoved"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Reference database used in FASTA format",
                    "label": "Reference database",
                    "inputBinding": {
                        "prefix": "-refdb"
                    },
                    "id": "#ngtax.cwl/reference_db"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Read length of the reverse read",
                    "label": "Reverse read length",
                    "inputBinding": {
                        "prefix": "-rev_read_len"
                    },
                    "id": "#ngtax.cwl/rev_read_len"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Reverse primer used",
                    "label": "The reverse primer used",
                    "inputBinding": {
                        "prefix": "-rev_p"
                    },
                    "id": "#ngtax.cwl/reverse_primer"
                },
                {
                    "type": "string",
                    "doc": "Name of the sample being analysed",
                    "label": "Sample name",
                    "id": "#ngtax.cwl/sample"
                }
            ],
            "baseCommand": [
                "java",
                "-jar",
                "/NGTax-2.2.9.jar",
                "-ngtax",
                "-mapFile",
                "cwl_mapping_file.txt"
            ],
            "arguments": [
                {
                    "prefix": "-t",
                    "valueFrom": "$(inputs.sample)_NG-Tax_$(inputs.for_read_len).ttl"
                },
                {
                    "prefix": "-b",
                    "valueFrom": "$(inputs.sample)_NG-Tax_$(inputs.for_read_len).biom"
                },
                {
                    "prefix": "-mock3",
                    "valueFrom": "$(inputs.mock3)"
                },
                {
                    "prefix": "-mock4",
                    "valueFrom": "$(inputs.mock4)"
                },
                {
                    "prefix": "-fragment",
                    "valueFrom": "$(inputs.fragment)"
                }
            ],
            "stdout": "ngtax2.stdout.log",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.sample)_NG-Tax_$(inputs.for_read_len).biom"
                    },
                    "id": "#ngtax.cwl/biom"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "ngtax2.stdout.log"
                    },
                    "id": "#ngtax.cwl/stdout_out"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.sample)_NG-Tax_$(inputs.for_read_len).ttl"
                    },
                    "id": "#ngtax.cwl/turtle"
                }
            ],
            "id": "#ngtax.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "https://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "https://schema.org/name": "Jasper Koehorst"
                }
            ],
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2020-00-00",
            "https://schema.org/dateModified": "2023-02-06",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "doc": "NGtax2 output conversion to prepare for biom file and ASV fasta file\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/scripts:1.0.1",
                    "class": "DockerRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "doc": "Fragment identifier that was used",
                    "label": "fragment identifier",
                    "id": "#ngtax_to_tsv-fasta.cwl/fragment"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#ngtax_to_tsv-fasta.cwl/identifier"
                },
                {
                    "type": "File",
                    "doc": "NGTax2 turtle output file",
                    "label": "NGTax2 turtle file",
                    "id": "#ngtax_to_tsv-fasta.cwl/input"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "UNLOCK assay metadata file",
                    "label": "Metadata file",
                    "id": "#ngtax_to_tsv-fasta.cwl/metadata"
                }
            ],
            "baseCommand": [
                "python3",
                "/scripts/ngtax_to_tsv-fasta.py"
            ],
            "arguments": [
                "-t",
                "$(inputs.input.path)",
                "-i",
                "$(inputs.identifier)",
                "-f",
                "$(inputs.fragment)",
                "-m",
                "$(inputs.metadata)"
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_asv.tsv"
                    },
                    "id": "#ngtax_to_tsv-fasta.cwl/physeq_asv"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_met.tsv"
                    },
                    "id": "#ngtax_to_tsv-fasta.cwl/physeq_met"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_seq.tsv"
                    },
                    "id": "#ngtax_to_tsv-fasta.cwl/physeq_seq"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*_tax.tsv"
                    },
                    "id": "#ngtax_to_tsv-fasta.cwl/physeq_tax"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.picrust.fasta"
                    },
                    "id": "#ngtax_to_tsv-fasta.cwl/picrust_fasta"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "*.picrust.tsv"
                    },
                    "id": "#ngtax_to_tsv-fasta.cwl/picrust_tsv"
                }
            ],
            "id": "#ngtax_to_tsv-fasta.cwl",
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
                    "type": "int",
                    "doc": "Read length of the reverse read",
                    "label": "Reverse read length",
                    "id": "#main/for_read_len"
                },
                {
                    "type": "string",
                    "doc": "Forward primer used",
                    "label": "Forward primer",
                    "id": "#main/forward_primer"
                },
                {
                    "type": "string",
                    "doc": "Subfragment that is being analysed (e.g. V1-V3 or V5-region)",
                    "label": "Subfragment name",
                    "id": "#main/fragment"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "UNLOCK assay metadata file",
                    "label": "Metadata file",
                    "id": "#main/metadata"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Mock3 reference selection",
                    "label": "Mock3 reference",
                    "id": "#main/mock3"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Mock4 reference selection",
                    "label": "Mock4 reference",
                    "id": "#main/mock4"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Reference database used in FASTA format",
                    "label": "Reference database",
                    "id": "#main/reference_db"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Read length of the reverse read",
                    "label": "Reverse read length",
                    "id": "#main/rev_read_len"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Reverse primer used",
                    "label": "Reverse primer",
                    "id": "#main/reverse_primer"
                },
                {
                    "type": "string",
                    "doc": "Name of the sample being analysed",
                    "label": "Sample name",
                    "id": "#main/sample"
                }
            ],
            "steps": [
                {
                    "run": "#ngtax.cwl",
                    "in": [
                        {
                            "source": "#main/for_read_len",
                            "id": "#main/ngtax/for_read_len"
                        },
                        {
                            "source": "#main/forward_primer",
                            "id": "#main/ngtax/forward_primer"
                        },
                        {
                            "source": "#main/fragment",
                            "id": "#main/ngtax/fragment"
                        },
                        {
                            "source": "#main/mock3",
                            "id": "#main/ngtax/mock3"
                        },
                        {
                            "source": "#main/mock4",
                            "id": "#main/ngtax/mock4"
                        },
                        {
                            "source": "#main/reference_db",
                            "id": "#main/ngtax/reference_db"
                        },
                        {
                            "source": "#main/rev_read_len",
                            "id": "#main/ngtax/rev_read_len"
                        },
                        {
                            "source": "#main/reverse_primer",
                            "id": "#main/ngtax/reverse_primer"
                        },
                        {
                            "source": "#main/sample",
                            "id": "#main/ngtax/sample"
                        }
                    ],
                    "out": [
                        "#main/ngtax/biom",
                        "#main/ngtax/turtle"
                    ],
                    "id": "#main/ngtax"
                },
                {
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"1_Classification\")",
                            "id": "#main/ngtax_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/ngtax/biom",
                                "#main/ngtax/turtle"
                            ],
                            "id": "#main/ngtax_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#main/ngtax_files_to_folder/results"
                    ],
                    "id": "#main/ngtax_files_to_folder"
                },
                {
                    "run": "#ngtax_to_tsv-fasta.cwl",
                    "in": [
                        {
                            "source": "#main/fragment",
                            "id": "#main/ngtax_to_tsv-fasta/fragment"
                        },
                        {
                            "source": "#main/sample",
                            "id": "#main/ngtax_to_tsv-fasta/identifier"
                        },
                        {
                            "source": "#main/ngtax/turtle",
                            "id": "#main/ngtax_to_tsv-fasta/input"
                        },
                        {
                            "source": "#main/metadata",
                            "id": "#main/ngtax_to_tsv-fasta/metadata"
                        }
                    ],
                    "out": [
                        "#main/ngtax_to_tsv-fasta/picrust_fasta",
                        "#main/ngtax_to_tsv-fasta/picrust_tsv",
                        "#main/ngtax_to_tsv-fasta/physeq_asv",
                        "#main/ngtax_to_tsv-fasta/physeq_seq",
                        "#main/ngtax_to_tsv-fasta/physeq_tax",
                        "#main/ngtax_to_tsv-fasta/physeq_met"
                    ],
                    "id": "#main/ngtax_to_tsv-fasta"
                },
                {
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"2_PHYLOSEQ\")",
                            "id": "#main/phyloseq_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/ngtax_to_tsv-fasta/physeq_asv",
                                "#main/ngtax_to_tsv-fasta/physeq_seq",
                                "#main/ngtax_to_tsv-fasta/physeq_tax",
                                "#main/ngtax_to_tsv-fasta/physeq_met"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/phyloseq_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#main/phyloseq_files_to_folder/results"
                    ],
                    "id": "#main/phyloseq_files_to_folder"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputSource": "#main/ngtax_files_to_folder/results",
                    "id": "#main/files_to_folder_ngtax"
                },
                {
                    "type": "Directory",
                    "outputSource": "#main/phyloseq_files_to_folder/results",
                    "id": "#main/files_to_folder_phyloseq"
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
