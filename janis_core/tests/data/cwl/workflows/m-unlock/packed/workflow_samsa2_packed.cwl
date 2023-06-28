{
    "$graph": [
        {
            "class": "CommandLineTool",
            "label": "compress a file multithreaded with pigz",
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/pigz:2.6",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.6"
                            ],
                            "specs": [
                                "https://anaconda.org/conda-forge/pigz"
                            ],
                            "package": "pigz"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "pigz",
                "-c"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.inputfile)"
                }
            ],
            "stdout": "$(inputs.inputfile.basename).gz",
            "inputs": [
                {
                    "type": "File",
                    "id": "#pigz.cwl/inputfile"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 1,
                    "inputBinding": {
                        "prefix": "-p"
                    },
                    "id": "#pigz.cwl/threads"
                }
            ],
            "id": "#pigz.cwl",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.inputfile.basename).gz"
                    },
                    "id": "#pigz.cwl/outfile"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "doc": "Diamond workflow implementation\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/diamond:2.0.15",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.0.15"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/diamond"
                            ],
                            "package": "diamond"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "doc": "Align amino acid (blastp) or DNA (blastx) sequences against a protein reference database",
                    "label": "Blast type",
                    "default": "blastx",
                    "id": "#diamond.cwl/align"
                },
                {
                    "type": "int",
                    "doc": "Block size in billions of sequence letters to be processed at a time. This is the main parameter for controlling the program\u2019s memory and disk space usage. Bigger numbers will increase the use of memory and temporary disk space, but also improve performance. The program can be expected to use roughly six times this number of memory (in GB).",
                    "label": "Block size",
                    "inputBinding": {
                        "prefix": "--block-size"
                    },
                    "default": 12,
                    "id": "#diamond.cwl/blocksize"
                },
                {
                    "type": "string",
                    "doc": "Path to the DIAMOND database file. Since v2.0.8, a BLAST database can also be used here.",
                    "label": "Database file path",
                    "inputBinding": {
                        "prefix": "--db"
                    },
                    "id": "#diamond.cwl/database"
                },
                {
                    "type": "File",
                    "doc": "forward sequence file locally",
                    "label": "forward reads",
                    "id": "#diamond.cwl/forward_reads"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#diamond.cwl/identifier"
                },
                {
                    "type": "int",
                    "doc": "The number of chunks for processing the seed index. This option can be additionally used to tune the performance. The default value is -c4, while setting this parameter to -c1 instead will improve the performance at the cost of increased memory use. Note that the very-sensitive and ultra-sensitive modes use -c1 by default.",
                    "label": "The number of chunks",
                    "inputBinding": {
                        "prefix": "--index-chunks"
                    },
                    "default": 1,
                    "id": "#diamond.cwl/indexchunks"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Scoring matrix. The following matrices are supported, with the default being BLOSUM62.",
                    "label": "scoring matrix",
                    "inputBinding": {
                        "prefix": "--matrix"
                    },
                    "default": "BLOSUM62",
                    "id": "#diamond.cwl/matrix"
                },
                {
                    "type": "int",
                    "doc": "maximum number of target sequences to report alignments for (default=25)",
                    "label": "Max target sequences",
                    "inputBinding": {
                        "prefix": "--max-target-seqs"
                    },
                    "default": 25,
                    "id": "#diamond.cwl/maxtargetseq"
                },
                {
                    "type": "string",
                    "doc": "Format of the output file. See the diamond manual for accepted output formats",
                    "label": "Output format",
                    "inputBinding": {
                        "prefix": "--outfmt"
                    },
                    "default": "100",
                    "id": "#diamond.cwl/outfmt"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Path to the output file. If this parameter is omitted, the results will be written to the standard output and all other program output will be suppressed.",
                    "label": "output file",
                    "inputBinding": {
                        "prefix": "--out"
                    },
                    "id": "#diamond.cwl/output"
                },
                {
                    "type": "File",
                    "doc": "reverse sequence file locally",
                    "label": "reverse reads",
                    "id": "#diamond.cwl/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 3,
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#diamond.cwl/threads"
                }
            ],
            "baseCommand": [
                "/unlock/infrastructure/binaries/diamond/diamond_v2.0.8.146/diamond"
            ],
            "arguments": [
                "$(inputs.align)",
                "--query",
                "$(inputs.forward_reads.path)",
                "$(inputs.reverse_reads.path)",
                "--salltitles"
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output).daa"
                    },
                    "id": "#diamond.cwl/output_diamond"
                }
            ],
            "id": "#diamond.cwl",
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
            "class": "CommandLineTool",
            "doc": "Diamond workflow implementation\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/diamond:2.0.15",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.0.15"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/diamond"
                            ],
                            "package": "diamond"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Diamond binary result file",
                    "label": "input file",
                    "inputBinding": {
                        "prefix": "--daa"
                    },
                    "id": "#view.cwl/inputfile"
                }
            ],
            "baseCommand": [
                "diamond"
            ],
            "arguments": [
                "view",
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
                "stitle",
                "--out",
                "$(inputs.inputfile.basename).tsv"
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.inputfile.basename).tsv"
                    },
                    "id": "#view.cwl/output_diamond_tabular"
                }
            ],
            "id": "#view.cwl",
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
            "outputs": [
                {
                    "type": "Directory",
                    "id": "#files_to_folder.cwl/results"
                }
            ],
            "id": "#files_to_folder.cwl",
            "http://schema.org/citation": "https://m-unlock.nl",
            "http://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "http://schema.org/dateCreated": "2020-00-00",
            "http://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "http://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "doc": "Samsa2 conversion workflow\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#convert.cwl/identifier"
                },
                {
                    "type": "File",
                    "doc": "diamond refseq or seed result table with salltitles",
                    "label": "diamond tabular file",
                    "id": "#convert.cwl/inputfile"
                }
            ],
            "baseCommand": [
                "python3"
            ],
            "arguments": [
                "/unlock/infrastructure/scripts/samsa2/manager.py",
                "$(inputs.inputfile.path)",
                "$(inputs.identifier)"
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_function.tsv"
                    },
                    "id": "#convert.cwl/output_refseq_function"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*_organism.tsv"
                    },
                    "id": "#convert.cwl/output_refseq_organism"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.hierarchy"
                    },
                    "id": "#convert.cwl/output_seed"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "*.reduced"
                    },
                    "id": "#convert.cwl/output_seed_reduced"
                }
            ],
            "id": "#convert.cwl",
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
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "MultipleInputFeatureRequirement"
                },
                {
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "label": "SAMSA2 pipeline",
            "doc": "SAMSA2 complete workflow for meta-omics read annotation\nSteps:\n  - Diamond read blastx\n    - Refseq\n    - SEED\n  - SAMSA2 processing\n",
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
                    "doc": "forward sequence file locally",
                    "label": "forward reads",
                    "id": "#main/forward_reads"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#main/identifier"
                },
                {
                    "type": "File",
                    "doc": "reverse sequence file locally",
                    "label": "reverse reads",
                    "id": "#main/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "CWL base step number",
                    "doc": "Step number for order of steps",
                    "default": 3,
                    "id": "#main/step"
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
            "outputs": [
                {
                    "label": "SAMSA2",
                    "doc": "functional and classification output folder by samsa2",
                    "type": "Directory",
                    "outputSource": "#main/samsa2_files_to_folder/results",
                    "id": "#main/samsa2_output"
                }
            ],
            "steps": [
                {
                    "label": "Compress large output files",
                    "doc": "Converts the diamond binary output file into a tabular output file",
                    "run": "#pigz.cwl",
                    "scatter": [
                        "#main/compress_diamond/inputfile"
                    ],
                    "scatterMethod": "dotproduct",
                    "in": [
                        {
                            "source": [
                                "#main/workflow_diamond_refseq/output_diamond",
                                "#main/workflow_diamond_seed/output_diamond"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/compress_diamond/inputfile"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/compress_diamond/threads"
                        }
                    ],
                    "out": [
                        "#main/compress_diamond/outfile"
                    ],
                    "id": "#main/compress_diamond"
                },
                {
                    "doc": "Preparation of samsa2 output files to a specific output folder",
                    "label": "SAMSA2 output files",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "${inputs.step+\"_SAMSA2\"}",
                            "id": "#main/samsa2_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/compress_diamond/outfile",
                                "#main/samsa2_postscripts/output_refseq_function",
                                "#main/samsa2_postscripts/output_refseq_organism",
                                "#main/samsa2_postscripts/output_seed",
                                "#main/samsa2_postscripts/output_seed_reduced"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/samsa2_files_to_folder/files"
                        },
                        {
                            "source": "#main/step",
                            "id": "#main/samsa2_files_to_folder/step"
                        }
                    ],
                    "out": [
                        "#main/samsa2_files_to_folder/results"
                    ],
                    "id": "#main/samsa2_files_to_folder"
                },
                {
                    "label": "Run SAMSA2 post scripts",
                    "doc": "Converts the diamond output file to samsa2 tables",
                    "run": "#convert.cwl",
                    "scatter": [
                        "#main/samsa2_postscripts/inputfile"
                    ],
                    "scatterMethod": "dotproduct",
                    "in": [
                        {
                            "source": "#main/identifier",
                            "id": "#main/samsa2_postscripts/identifier"
                        },
                        {
                            "source": [
                                "#main/workflow_diamond_view/output_diamond_tabular"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/samsa2_postscripts/inputfile"
                        }
                    ],
                    "out": [
                        "#main/samsa2_postscripts/output_refseq_function",
                        "#main/samsa2_postscripts/output_refseq_organism",
                        "#main/samsa2_postscripts/output_seed",
                        "#main/samsa2_postscripts/output_seed_reduced"
                    ],
                    "id": "#main/samsa2_postscripts"
                },
                {
                    "label": "Diamond refseq workflow",
                    "doc": "Read mapping using the refseq database",
                    "run": "#diamond.cwl",
                    "in": [
                        {
                            "default": "/unlock/references/databases/ncbi/Refseq_Bacterial/diamond/ncbi-bact-refseq_28-01-2020_proteins.dmnd",
                            "id": "#main/workflow_diamond_refseq/database"
                        },
                        {
                            "source": "#main/forward_reads",
                            "id": "#main/workflow_diamond_refseq/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_diamond_refseq/identifier"
                        },
                        {
                            "default": 1,
                            "id": "#main/workflow_diamond_refseq/maxtargetseq"
                        },
                        {
                            "valueFrom": "$(inputs.identifier)_refseq",
                            "id": "#main/workflow_diamond_refseq/output"
                        },
                        {
                            "source": "#main/reverse_reads",
                            "id": "#main/workflow_diamond_refseq/reverse_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_diamond_refseq/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_diamond_refseq/output_diamond"
                    ],
                    "id": "#main/workflow_diamond_refseq"
                },
                {
                    "label": "Diamond seed workflow",
                    "doc": "Read mapping using the seed database",
                    "run": "#diamond.cwl",
                    "in": [
                        {
                            "default": "/unlock/references/databases/SEED/diamond/seed_subsystems_db.dmnd",
                            "id": "#main/workflow_diamond_seed/database"
                        },
                        {
                            "source": "#main/forward_reads",
                            "id": "#main/workflow_diamond_seed/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_diamond_seed/identifier"
                        },
                        {
                            "default": 1,
                            "id": "#main/workflow_diamond_seed/maxtargetseq"
                        },
                        {
                            "valueFrom": "$(inputs.identifier)_seed",
                            "id": "#main/workflow_diamond_seed/output"
                        },
                        {
                            "source": "#main/reverse_reads",
                            "id": "#main/workflow_diamond_seed/reverse_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_diamond_seed/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_diamond_seed/output_diamond"
                    ],
                    "id": "#main/workflow_diamond_seed"
                },
                {
                    "label": "Change view of diamond binary to tabular",
                    "doc": "Converts the diamond binary output file into a tabular output file",
                    "run": "#view.cwl",
                    "scatter": [
                        "#main/workflow_diamond_view/inputfile"
                    ],
                    "scatterMethod": "dotproduct",
                    "in": [
                        {
                            "source": [
                                "#main/workflow_diamond_refseq/output_diamond",
                                "#main/workflow_diamond_seed/output_diamond"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/workflow_diamond_view/inputfile"
                        }
                    ],
                    "out": [
                        "#main/workflow_diamond_view/output_diamond_tabular"
                    ],
                    "id": "#main/workflow_diamond_view"
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
            "https://schema.org/dateCreated": "2021-09-00",
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
