{
    "$graph": [
        {
            "class": "CommandLineTool",
            "label": "Compress a directory (tar)",
            "hints": [
                {
                    "dockerPull": "debian:buster",
                    "class": "DockerRequirement"
                }
            ],
            "baseCommand": [
                "tar",
                "czfh"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.indir.basename).tar.gz"
                },
                {
                    "valueFrom": "-C"
                },
                {
                    "valueFrom": "$(inputs.indir.path)/.."
                },
                {
                    "valueFrom": "$(inputs.indir.basename)"
                }
            ],
            "inputs": [
                {
                    "type": "Directory",
                    "id": "#compress_directory.cwl/indir"
                }
            ],
            "id": "#compress_directory.cwl",
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
                        "glob": "$(inputs.indir.basename).tar.gz"
                    },
                    "id": "#compress_directory.cwl/outfile"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "label": "BUSCO",
            "doc": "Based on evolutionarily-informed expectations of gene content of near-universal single-copy orthologs, \nBUSCO metric is complementary to technical metrics like N50.\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "networkAccess": "$(inputs.busco_data !== undefined)",
                    "class": "NetworkAccess"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/busco:5.4.4",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "5.4.4"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/busco"
                            ],
                            "package": "busco"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "busco"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Auto-lineage detection",
                    "doc": "Run auto-lineage to find optimum lineage path",
                    "inputBinding": {
                        "prefix": "--auto-lineage"
                    },
                    "id": "#busco.cwl/auto-lineage"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Eukaryote auto-lineage detection",
                    "doc": "Run auto-placement just on eukaryote tree to find optimum lineage path.",
                    "inputBinding": {
                        "prefix": "--auto-lineage-euk"
                    },
                    "id": "#busco.cwl/auto-lineage-euk"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Prokaryote auto-lineage detection",
                    "doc": "Run auto-lineage just on non-eukaryote trees to find optimum lineage path.",
                    "inputBinding": {
                        "prefix": "--auto-lineage-prok"
                    },
                    "id": "#busco.cwl/auto-lineage-prok"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "Dataset location",
                    "doc": "This assumes --offline mode. Specify local filepath for finding BUSCO dataset downloads",
                    "inputBinding": {
                        "prefix": "--download_path"
                    },
                    "id": "#busco.cwl/busco_data"
                },
                {
                    "type": "string",
                    "label": "Name of the output file",
                    "doc": "Give your analysis run a recognisable short name. Output folders and files will be labelled with this name.",
                    "inputBinding": {
                        "prefix": "--out"
                    },
                    "id": "#busco.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Lineage",
                    "doc": "Specify the name of the BUSCO lineage to be used.",
                    "inputBinding": {
                        "prefix": "--lineage_dataset"
                    },
                    "id": "#busco.cwl/lineage"
                },
                {
                    "type": "string",
                    "label": "Input molecule type",
                    "doc": "Specify which BUSCO analysis mode to run.\nThere are three valid modes:\n- geno or genome, for genome assemblies (DNA)\n- tran or transcriptome, for transcriptome assemblies (DNA)\n- prot or proteins, for annotated gene sets (protein)\n",
                    "inputBinding": {
                        "prefix": "--mode"
                    },
                    "id": "#busco.cwl/mode"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Input fasta file",
                    "doc": "Input sequence file in FASTA format. Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. Also possible to use a path to a directory containing multiple input files.",
                    "inputBinding": {
                        "prefix": "--in"
                    },
                    "id": "#busco.cwl/sequence_file"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "Input folder",
                    "doc": "Input folder with sequence files in FASTA format. Can be an assembled genome or transcriptome (DNA), or protein sequences from an annotated gene set. Also possible to use a path to a directory containing multiple input files.",
                    "inputBinding": {
                        "prefix": "--in"
                    },
                    "id": "#busco.cwl/sequence_folder"
                },
                {
                    "type": "boolean",
                    "label": "Compress output",
                    "doc": "Compress some subdirectories with many files to save space",
                    "default": true,
                    "inputBinding": {
                        "prefix": "--tar"
                    },
                    "id": "#busco.cwl/tar_output"
                },
                {
                    "type": "int",
                    "label": "Number of threads",
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--cpu"
                    },
                    "id": "#busco.cwl/threads"
                }
            ],
            "arguments": [
                "${\n  if (inputs.busco_data){\n    return '--offline';\n  } else {\n    return null;\n  }\n}\n"
            ],
            "outputs": [
                {
                    "label": "Batch summary",
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Summary file when input is multiple files",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)/batch_summary.txt"
                    },
                    "id": "#busco.cwl/batch_summary"
                },
                {
                    "label": "BUSCO logs folder",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.identifier)/logs"
                    },
                    "id": "#busco.cwl/logs"
                },
                {
                    "label": "BUSCO short summary files",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.identifier)/short_summary.*"
                    },
                    "id": "#busco.cwl/short_summaries"
                }
            ],
            "id": "#busco.cwl",
            "https://schema.org/author": [
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0002-5516-8391",
                    "https://schema.org/email": "mailto:german.royvalgarcia@wur.nl",
                    "https://schema.org/name": "Germ\u00e1n Royval"
                },
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
            "https://schema.org/dateCreated": "2022-01-01",
            "https://schema.org/dateModified": "2022-02-28",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "CheckM",
            "doc": "CheckM provides a set of tools for assessing the quality of genomes recovered from isolates, single cells, or metagenomes\n",
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/checkm-genome:1.2.2",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "1.2.1"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/checkm-genome"
                            ],
                            "package": "checkm-genome"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "requirements": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/checkm-genome:1.2.2",
                    "class": "DockerRequirement"
                },
                {
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "/checkm_data",
                            "writable": true
                        },
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "checkm_output",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "# !/bin/bash\nexport CHECKM_DATA_PATH=/venv/checkm_data\ncheckm lineage_wf $@"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "arguments": [
                {
                    "position": 51,
                    "valueFrom": "--reduced_tree"
                },
                {
                    "position": 52,
                    "prefix": "-f",
                    "valueFrom": "$(inputs.identifier)_CheckM_report.txt"
                },
                {
                    "position": 54,
                    "valueFrom": "$(inputs.identifier)_checkm"
                }
            ],
            "baseCommand": [
                "bash",
                "-x",
                "script.sh"
            ],
            "inputs": [
                {
                    "type": "Directory",
                    "doc": "folder containing bins in fasta format from metagenomic binning",
                    "label": "bins folder",
                    "inputBinding": {
                        "position": 53
                    },
                    "id": "#checkm_lineagewf.cwl/bin_dir"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "fasta file extension",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "-x"
                    },
                    "default": "fa",
                    "id": "#checkm_lineagewf.cwl/fasta_extension"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#checkm_lineagewf.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Number of threads to use",
                    "default": 8,
                    "inputBinding": {
                        "position": 4,
                        "prefix": "-t"
                    },
                    "id": "#checkm_lineagewf.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_checkm"
                    },
                    "id": "#checkm_lineagewf.cwl/checkm_out_folder"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_CheckM_report.txt"
                    },
                    "id": "#checkm_lineagewf.cwl/checkm_out_table"
                }
            ],
            "id": "#checkm_lineagewf.cwl"
        },
        {
            "class": "CommandLineTool",
            "label": "DAS Tool",
            "doc": "Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy.",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/das_tool:1.1.5",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "1.1.5"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/das_tool"
                            ],
                            "package": "dastool"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "DAS_Tool"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Input assembly in fasta format",
                    "label": "Input assembly",
                    "inputBinding": {
                        "prefix": "--contigs"
                    },
                    "id": "#das_tool.cwl/assembly"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Comma separated list of tab separated contigs to bin tables.",
                    "label": "Bin-Contig tables",
                    "inputBinding": {
                        "itemSeparator": ",",
                        "prefix": "--bins"
                    },
                    "id": "#das_tool.cwl/bin_tables"
                },
                {
                    "type": "string",
                    "doc": "Comma separated list of the binning prediction tool names.",
                    "label": "Binner labels",
                    "inputBinding": {
                        "prefix": "--labels"
                    },
                    "id": "#das_tool.cwl/binner_labels"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#das_tool.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Number of threads to use",
                    "default": 1,
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--threads"
                    },
                    "id": "#das_tool.cwl/threads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Export bins as fasta files.",
                    "label": "Write bins",
                    "inputBinding": {
                        "prefix": "--write_bins"
                    },
                    "default": true,
                    "id": "#das_tool.cwl/write_bins"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Export unbinned contigs as fasta file",
                    "label": "Write unbinned",
                    "inputBinding": {
                        "prefix": "--write_unbinned"
                    },
                    "default": true,
                    "id": "#das_tool.cwl/write_unbinned"
                }
            ],
            "arguments": [
                {
                    "prefix": "-o",
                    "valueFrom": "$(inputs.identifier)"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "Bins",
                    "doc": "Bin fasta files.",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_DASTool_bins"
                    },
                    "id": "#das_tool.cwl/bin_dir"
                },
                {
                    "type": "File",
                    "label": "Contig to bin",
                    "doc": "Contigs to bin file table",
                    "outputBinding": {
                        "glob": "*_DASTool_contig2bin.tsv"
                    },
                    "id": "#das_tool.cwl/contig2bin"
                },
                {
                    "type": "File",
                    "label": "Log",
                    "doc": "DASTool log file",
                    "outputBinding": {
                        "glob": "*_DASTool.log"
                    },
                    "id": "#das_tool.cwl/log"
                },
                {
                    "type": "File",
                    "label": "DAS Tool run summary",
                    "doc": "Summary",
                    "outputBinding": {
                        "glob": "*_DASTool_summary.tsv"
                    },
                    "id": "#das_tool.cwl/summary"
                }
            ],
            "id": "#das_tool.cwl",
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
            "https://schema.org/dateCreated": "2022-09-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "Fasta_to_Scaffolds2Bin",
            "doc": "Converts genome bins in fasta format to scaffolds-to-bin table. (DAS Tool helper script)",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/das_tool:1.1.4",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "1.1.4"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/das_tool"
                            ],
                            "package": "dastool"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "Fasta_to_Contig2Bin.sh"
            ],
            "inputs": [
                {
                    "type": "Directory",
                    "doc": "Input assembly in fasta format",
                    "label": "Input assembly",
                    "inputBinding": {
                        "prefix": "--input_folder"
                    },
                    "id": "#fasta_to_contig2bin.cwl/bin_folder"
                },
                {
                    "type": "string",
                    "doc": "Binner name used to create the bins",
                    "label": "Binner name",
                    "id": "#fasta_to_contig2bin.cwl/binner_name"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Extension of fasta files. (default fasta)",
                    "label": "Fasta extension",
                    "inputBinding": {
                        "prefix": "--extension"
                    },
                    "id": "#fasta_to_contig2bin.cwl/extension"
                }
            ],
            "stdout": "$(inputs.binner_name)_Contig2Bin.tsv",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.binner_name)_Contig2Bin.tsv"
                    },
                    "id": "#fasta_to_contig2bin.cwl/table"
                }
            ],
            "id": "#fasta_to_contig2bin.cwl",
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
            "https://schema.org/dateCreated": "2022-09-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "EukRep",
            "doc": "EukRep, Classification of Eukaryotic and Prokaryotic sequences from metagenomic datasets",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/eukrep:0.6.7",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "0.6.7"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/eukrep"
                            ],
                            "package": "diamond"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "EukRep"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Input assembly in fasta format",
                    "label": "Input assembly",
                    "inputBinding": {
                        "prefix": "-i"
                    },
                    "id": "#eukrep.cwl/assembly"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#eukrep.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Minumum contig length",
                    "doc": "Minimum sequence length cutoff for sequences to be included in prediction. Default is 3kb",
                    "id": "#eukrep.cwl/min_contig_size"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "{strict,balanced,lenient} Default is balanced.\nHow stringent the algorithm is in identifying eukaryotic scaffolds. Strict has a lower false positive rate and true positive rate; vice verso for leneient.\n",
                    "label": "Algorithm stringency",
                    "inputBinding": {
                        "prefix": "-m"
                    },
                    "id": "#eukrep.cwl/stringency"
                }
            ],
            "arguments": [
                {
                    "prefix": "-o",
                    "valueFrom": "$(inputs.identifier)_EukRep.fasta"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_EukRep.fasta"
                    },
                    "id": "#eukrep.cwl/euk_fasta_out"
                }
            ],
            "id": "#eukrep.cwl",
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
            "https://schema.org/dateCreated": "2022-06-00",
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
            "class": "ExpressionTool",
            "doc": "Expression to filter files (by name) in a directory using a regular expression.\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "loadListing": "shallow_listing",
                    "class": "LoadListingRequirement"
                }
            ],
            "inputs": [
                {
                    "label": "Input folder",
                    "doc": "Folder with only files",
                    "type": "Directory",
                    "id": "#folder_file_regex.cwl/folder"
                },
                {
                    "type": "string",
                    "label": "Output folder name",
                    "doc": "Output folder name",
                    "id": "#folder_file_regex.cwl/output_folder_name"
                },
                {
                    "label": "Regex (JS)",
                    "doc": "JavaScript regular expression to be used on the filenames\nMetaBAT2 example: \"bin\\.[0-9]+\\.fa\"\n",
                    "type": "string",
                    "id": "#folder_file_regex.cwl/regex"
                }
            ],
            "expression": "${\n  var regex = new RegExp(inputs.regex);\n  var array = [];\n  for (var i = 0; i < inputs.folder.listing.length; i++) {\n    if (regex.test(inputs.folder.listing[i].location)){\n      array = array.concat(inputs.folder.listing[i]);\n    }\n  }\n  var r = {\n     'output_folder':\n       { \"class\": \"Directory\",\n         \"basename\": inputs.output_folder_name,\n         \"listing\": array\n       }\n     };\n   return r;\n}\n",
            "outputs": [
                {
                    "type": "Directory",
                    "id": "#folder_file_regex.cwl/output_folder"
                }
            ],
            "id": "#folder_file_regex.cwl",
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
            "https://schema.org/dateCreated": "2022-10-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "ExpressionTool",
            "doc": "Transforms the input folder to an array of files\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "loadListing": "shallow_listing",
                    "class": "LoadListingRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "Directory",
                    "id": "#folder_to_files.cwl/folder"
                }
            ],
            "expression": "${\n  var files = [];\n  for (var i = 0; i < inputs.folder.listing.length; i++) {\n    files.push(inputs.folder.listing[i]);\n  }\n  return {\"files\": files};\n}  \n",
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#folder_to_files.cwl/files"
                }
            ],
            "id": "#folder_to_files.cwl",
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
            "https://schema.org/dateCreated": "2022-10-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "GTDBTK Classify Workflow",
            "doc": "Taxonomic genome classification workflow with GTDBTK. \n",
            "baseCommand": [
                "bash",
                "script.sh"
            ],
            "requirements": [
                {
                    "listing": [
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nexport GTDBTK_DATA_PATH=$1\nshift;\ngtdbtk classify_wf $@"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/gtdbtk:2.1.1",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.1.1"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/gtdbtk"
                            ],
                            "package": "gtdbtk"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "Directory",
                    "doc": "Directory containing bins in fasta format from metagenomic binning",
                    "label": "bins with directory",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--genome_dir"
                    },
                    "id": "#gtdbtk_classify_wf.cwl/bin_dir"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "fasta file extension",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--extension"
                    },
                    "default": "fa",
                    "id": "#gtdbtk_classify_wf.cwl/fasta_extension"
                },
                {
                    "type": "Directory",
                    "doc": "Directory containing the GTDBTK repository",
                    "label": "gtdbtk data directory",
                    "loadListing": "no_listing",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#gtdbtk_classify_wf.cwl/gtdbtk_data"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#gtdbtk_classify_wf.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Number of threads to use",
                    "default": 8,
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--cpus"
                    },
                    "id": "#gtdbtk_classify_wf.cwl/threads"
                }
            ],
            "arguments": [
                {
                    "valueFrom": "--force",
                    "position": 10
                },
                {
                    "prefix": "--prefix",
                    "valueFrom": "$(inputs.identifier).gtdbtk",
                    "position": 11
                },
                {
                    "prefix": "--out_dir",
                    "position": 12,
                    "valueFrom": "$(inputs.identifier)_GTDB-Tk"
                }
            ],
            "stdout": "$(inputs.identifier)_GTDB-Tk.stdout.log",
            "outputs": [
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_GTDB-Tk"
                    },
                    "id": "#gtdbtk_classify_wf.cwl/gtdbtk_out_folder"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_GTDB-Tk/classify/$(inputs.identifier).gtdbtk.bac120.summary.tsv"
                    },
                    "id": "#gtdbtk_classify_wf.cwl/gtdbtk_summary"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_GTDB-Tk.stdout.log"
                    },
                    "id": "#gtdbtk_classify_wf.cwl/stdout_out"
                }
            ],
            "id": "#gtdbtk_classify_wf.cwl",
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
            "https://schema.org/dateModified": "2022-02-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "MaxBin2",
            "doc": "MaxBin2 is a software for binning assembled metagenomic sequences based on an Expectation-Maximization algorithm.",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/maxbin2:2.2.7",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.2.7"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/maxbin2"
                            ],
                            "package": "maxbin2"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "run_MaxBin.pl"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Abundances file",
                    "label": "Abundances",
                    "inputBinding": {
                        "prefix": "-abund"
                    },
                    "id": "#maxbin2.cwl/abundances"
                },
                {
                    "type": "File",
                    "doc": "Input assembly in fasta format",
                    "label": "Input assembly",
                    "inputBinding": {
                        "prefix": "-contig"
                    },
                    "id": "#maxbin2.cwl/contigs"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#maxbin2.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Number of threads to use",
                    "default": 1,
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-thread"
                    },
                    "id": "#maxbin2.cwl/threads"
                }
            ],
            "arguments": [
                {
                    "prefix": "-out",
                    "valueFrom": "$(inputs.identifier)_MaxBin2.bin"
                }
            ],
            "stdout": "$(inputs.identifier)_MaxBin2.log",
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "label": "Bins",
                    "doc": "Bin fasta files. The XX bin. XX are numbers, e.g. out.001.fasta",
                    "outputBinding": {
                        "glob": "*.fasta"
                    },
                    "id": "#maxbin2.cwl/bins"
                },
                {
                    "type": "File",
                    "label": "Log",
                    "doc": "Log file recording the core steps of MaxBin algorithm",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_MaxBin2.log"
                    },
                    "id": "#maxbin2.cwl/log"
                },
                {
                    "type": "File",
                    "label": "Markers",
                    "doc": "Marker gene presence numbers for each bin. This table is ready to be plotted by R or other 3rd-party software.",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_MaxBin2.bin.marker"
                    },
                    "id": "#maxbin2.cwl/markers"
                },
                {
                    "type": "File",
                    "label": "MaxBin2 Summary",
                    "doc": "Summary file describing which contigs are being classified into which bin.",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_MaxBin2.bin.summary"
                    },
                    "id": "#maxbin2.cwl/summary"
                }
            ],
            "id": "#maxbin2.cwl",
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
            "https://schema.org/dateCreated": "2022-08-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "aggregateBinDepths",
            "doc": "Aggregate bin depths using MetaBat2 using the script aggregateBinDepths.pl\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/metabat2:2.15",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.15"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/metabat2"
                            ],
                            "package": "metabat2"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Bin fasta files",
                    "label": "Bin fasta files",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#aggregateBinDepths.cwl/bins"
                },
                {
                    "type": "string",
                    "doc": "Name of the output file",
                    "label": "output file name",
                    "id": "#aggregateBinDepths.cwl/identifier"
                },
                {
                    "type": "File",
                    "doc": "Contig depths files obtained from metabat2 script jgi_summarize_bam_contig_depths",
                    "label": "contigs depths",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#aggregateBinDepths.cwl/metabatdepthsFile"
                }
            ],
            "baseCommand": [
                "aggregateBinDepths.pl"
            ],
            "stdout": "$(inputs.identifier)_binDepths.tsv",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_binDepths.tsv"
                    },
                    "id": "#aggregateBinDepths.cwl/binDepths"
                }
            ],
            "id": "#aggregateBinDepths.cwl",
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
            "label": "MetaBAT2 binning",
            "doc": "Metagenome Binning based on Abundance and Tetranucleotide frequency (MetaBat2)\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/metabat2:2.15",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.15"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/metabat2"
                            ],
                            "package": "metabat2"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "label": "The input assembly in fasta format",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--inFile"
                    },
                    "id": "#metabat2.cwl/assembly"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 5,
                        "prefix": "--abdFile"
                    },
                    "id": "#metabat2.cwl/depths"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#metabat2.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Number of threads to use",
                    "default": 1,
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--numThreads"
                    },
                    "id": "#metabat2.cwl/threads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Export unbinned contigs as fasta file",
                    "label": "Write unbinned",
                    "inputBinding": {
                        "prefix": "--unbinned"
                    },
                    "id": "#metabat2.cwl/write_unbinned"
                }
            ],
            "arguments": [
                {
                    "prefix": "--outFile",
                    "valueFrom": "MetaBAT2_bins/$(inputs.identifier)_MetaBAT2_bin"
                }
            ],
            "baseCommand": [
                "metabat2"
            ],
            "stdout": "$(inputs.identifier)_MetaBAT2.log",
            "outputs": [
                {
                    "type": "Directory",
                    "label": "Bin directory",
                    "doc": "Bin directory",
                    "outputBinding": {
                        "glob": "MetaBAT2_bins"
                    },
                    "id": "#metabat2.cwl/bin_dir"
                },
                {
                    "type": "File",
                    "label": "Log",
                    "doc": "MetaBat2 log file",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_MetaBAT2.log"
                    },
                    "id": "#metabat2.cwl/log"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Unbinned contigs",
                    "doc": "Unbinned contig fasta files",
                    "outputBinding": {
                        "glob": "MetaBAT2_bins/$(inputs.identifier)_MetaBAT2_bin.unbinned.fa"
                    },
                    "id": "#metabat2.cwl/unbinned"
                }
            ],
            "id": "#metabat2.cwl",
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
            "https://schema.org/dateCreated": "2022-10-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "jgi_summarize_bam_contig_depths",
            "doc": "Summarize contig read depth from bam file for metabat2 binning.\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/metabat2:2.15",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.15"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/metabat2"
                            ],
                            "package": "metabat2"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#metabatContigDepths.cwl/bamFile"
                },
                {
                    "type": "string",
                    "doc": "Name of the output file",
                    "label": "output file name",
                    "id": "#metabatContigDepths.cwl/identifier"
                }
            ],
            "baseCommand": [
                "jgi_summarize_bam_contig_depths"
            ],
            "arguments": [
                {
                    "position": 1,
                    "prefix": "--outputDepth",
                    "valueFrom": "$(inputs.identifier)_contigDepths.tsv"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_contigDepths.tsv"
                    },
                    "id": "#metabatContigDepths.cwl/depths"
                }
            ],
            "id": "#metabatContigDepths.cwl"
        },
        {
            "class": "CommandLineTool",
            "label": "Bin read mapping stats",
            "doc": "Table of general read mapping statistics of the bins and assembly\n\nID\nReads\nAssembly size\nContigs\nn50\nLargest contig\nMapped reads\nBins\nTotal bin size\nBinned\nReads mapped to bins\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/scripts:1.0.3",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "0.20.0"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/pysam"
                            ],
                            "package": "pysam"
                        },
                        {
                            "version": [
                                "3.10.6"
                            ],
                            "specs": [
                                "https://anaconda.org/conda-forge/python"
                            ],
                            "package": "python3"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "python3",
                "/scripts/metagenomics/assembly_bins_readstats.py"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Assembly in fasta format",
                    "label": "Assembly",
                    "inputBinding": {
                        "prefix": "--assembly"
                    },
                    "id": "#assembly_bins_readstats.cwl/assembly"
                },
                {
                    "type": "File",
                    "doc": "BAM file with reads mapped to the assembly",
                    "label": "BAM file",
                    "inputBinding": {
                        "prefix": "--bam"
                    },
                    "id": "#assembly_bins_readstats.cwl/bam_file"
                },
                {
                    "type": "File",
                    "doc": "File containing bins names and their respective (assembly) contigs. Format contig<tab>bin_name",
                    "label": "binContigs file",
                    "inputBinding": {
                        "prefix": "--binContigs"
                    },
                    "id": "#assembly_bins_readstats.cwl/binContigs"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "inputBinding": {
                        "prefix": "--identifier"
                    },
                    "id": "#assembly_bins_readstats.cwl/identifier"
                }
            ],
            "stdout": "$(inputs.identifier)_binReadStats.tsv",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_binReadStats.tsv"
                    },
                    "id": "#assembly_bins_readstats.cwl/binReadStats"
                }
            ],
            "id": "#assembly_bins_readstats.cwl",
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
            "https://schema.org/dateCreated": "2022-12-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "Bin summary",
            "doc": "Creates a summary table of the bins and their quality and taxonomy.\n\nColumns are:\nBin\nContigs\nSize\nLargest_contig\nN50\nGC\navgDepth\nGTDB-Tk_taxonomy\nBUSCO_Taxonomy\nBUSCO_score\nCheckM_Completeness\nCheckM_Contamination\nCheckM_Strain-heterogeneity    \n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/scripts:1.0.3",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "1.5.0"
                            ],
                            "specs": [
                                "https://anaconda.org/conda-forge/pandas"
                            ],
                            "package": "pandas"
                        },
                        {
                            "version": [
                                "3.10.6"
                            ],
                            "specs": [
                                "https://anaconda.org/conda-forge/python"
                            ],
                            "package": "python3"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "python3",
                "/scripts/metagenomics/bins_summary.py"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "MetaBAT2 aggregateDepths file",
                    "label": "bin depths",
                    "inputBinding": {
                        "prefix": "--bin_depths"
                    },
                    "id": "#bins_summary.cwl/bin_depths"
                },
                {
                    "type": "Directory",
                    "doc": "Directory containing bins in fasta format from metagenomic binning",
                    "label": "Bins directory",
                    "inputBinding": {
                        "prefix": "--bin_dir"
                    },
                    "id": "#bins_summary.cwl/bin_dir"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Directory containing BUSCO reports",
                    "label": "BUSCO folder",
                    "inputBinding": {
                        "prefix": "--busco_batch"
                    },
                    "id": "#bins_summary.cwl/busco_batch"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "CheckM report file",
                    "label": "CheckM report",
                    "inputBinding": {
                        "prefix": "--checkm"
                    },
                    "id": "#bins_summary.cwl/checkm"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "CheckM report file",
                    "label": "CheckM report",
                    "inputBinding": {
                        "prefix": "--gtdbtk"
                    },
                    "id": "#bins_summary.cwl/gtdbtk"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#bins_summary.cwl/identifier"
                }
            ],
            "arguments": [
                {
                    "prefix": "--output",
                    "valueFrom": "$(inputs.identifier)"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_binContigs.tsv"
                    },
                    "id": "#bins_summary.cwl/bin_contigs"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_binSummary.tsv"
                    },
                    "id": "#bins_summary.cwl/bins_summary_table"
                }
            ],
            "id": "#bins_summary.cwl",
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
            "https://schema.org/dateModified": "2022-12-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "SemiBin",
            "doc": "Metagenomic binning with semi-supervised deep learning using information from reference genomes.",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "networkAccess": true,
                    "class": "NetworkAccess"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/semibin:1.4.0",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "1.4.0"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/semibin"
                            ],
                            "package": "semibin"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "SemiBin single_easy_bin"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Input assembly in fasta format",
                    "label": "Input assembly",
                    "inputBinding": {
                        "prefix": "--input-fasta"
                    },
                    "id": "#semibin_single_easy_bin.cwl/assembly"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Mapped reads to assembly in sorted BAM format",
                    "label": "BAM file",
                    "inputBinding": {
                        "prefix": "--input-bam"
                    },
                    "id": "#semibin_single_easy_bin.cwl/bam_file"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Built-in models (human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global)",
                    "label": "Environment",
                    "inputBinding": {
                        "prefix": "--environment"
                    },
                    "id": "#semibin_single_easy_bin.cwl/environment"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#semibin_single_easy_bin.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "Contig depth file from MetaBAT2",
                    "label": "MetaBAT2 depths",
                    "inputBinding": {
                        "prefix": "--depth-metabat2"
                    },
                    "id": "#semibin_single_easy_bin.cwl/metabat2_depth_file"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "doc": "Reference Database data directory (usually, MMseqs2 GTDB)",
                    "label": "Reference Database",
                    "inputBinding": {
                        "prefix": "--reference-db"
                    },
                    "id": "#semibin_single_easy_bin.cwl/reference_database"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "An alternative binning algorithm for assemblies from long-read datasets.",
                    "label": "Long read assembly",
                    "inputBinding": {
                        "prefix": "--sequencing-type=long_read"
                    },
                    "id": "#semibin_single_easy_bin.cwl/sequencing_type_longread"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Number of threads to use",
                    "default": 1,
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--threads"
                    },
                    "id": "#semibin_single_easy_bin.cwl/threads"
                }
            ],
            "arguments": [
                {
                    "prefix": "-o",
                    "valueFrom": "$(inputs.identifier)_SemiBin"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Coverage data",
                    "doc": "Coverage data generated from depth file.",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SemiBin/\"*_cov.csv\""
                    },
                    "id": "#semibin_single_easy_bin.cwl/coverage"
                },
                {
                    "type": "File",
                    "label": "Training data",
                    "doc": "Data used in the training of deep learning model",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SemiBin/data.csv"
                    },
                    "id": "#semibin_single_easy_bin.cwl/data"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Training data",
                    "doc": "Data used in the training of deep learning model",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SemiBin/data_split.csv"
                    },
                    "id": "#semibin_single_easy_bin.cwl/data_split"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Bins info",
                    "doc": "Info on (reclustered) bins (contig,nbs,n50 etc..)",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SemiBin/recluster_bins_info.tsv"
                    },
                    "id": "#semibin_single_easy_bin.cwl/info"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "MMseqs annotation",
                    "doc": "MMseqs contig annotation",
                    "outputBinding": {
                        "glob": "mmseqs_contig_annotation"
                    },
                    "id": "#semibin_single_easy_bin.cwl/mmseqs_contig_annotation"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Deep learning model",
                    "doc": "Saved semi-supervised deep learning model.",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SemiBin/model.h5"
                    },
                    "id": "#semibin_single_easy_bin.cwl/model"
                },
                {
                    "type": "Directory",
                    "label": "Bins",
                    "doc": "Directory of all reconstructed bins before reclustering.",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SemiBin/output_prerecluster_bins"
                    },
                    "id": "#semibin_single_easy_bin.cwl/output_bins"
                },
                {
                    "type": "Directory",
                    "label": "Reclustered Bins",
                    "doc": "Directory of all reconstructed bins after reclustering",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SemiBin/output_recluster_bins"
                    },
                    "id": "#semibin_single_easy_bin.cwl/recluster_bins"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "Markers",
                    "doc": "Directory with HMM marker hits",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SemiBin/sample0"
                    },
                    "id": "#semibin_single_easy_bin.cwl/sample0"
                }
            ],
            "id": "#semibin_single_easy_bin.cwl",
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
            "https://schema.org/dateCreated": "2022-09-00",
            "https://schema.org/dateModified": "2022-12-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "Fasta statistics",
            "doc": "Fasta statistics like N50, total length, etc..",
            "hints": [
                {
                    "dockerPull": "quay.io/biocontainers/idba:1.1.3--1",
                    "class": "DockerRequirement"
                }
            ],
            "baseCommand": [
                "raw_n50"
            ],
            "stdout": "$(inputs.identifier)_stats.txt",
            "inputs": [
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#raw_n50.cwl/identifier"
                },
                {
                    "type": "File",
                    "label": "Input fasta",
                    "doc": "Input multi fasta file",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#raw_n50.cwl/input_fasta"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_stats.txt"
                    },
                    "id": "#raw_n50.cwl/output"
                }
            ],
            "id": "#raw_n50.cwl",
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
            "https://schema.org/dateCreated": "2022-00-06",
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
            "label": "Metagenomic Binning from Assembly",
            "doc": "Workflow for Metagenomics binning from assembly.<br>\n\nMinimal inputs are: Identifier, assembly (fasta) and a associated sorted BAM file\n\nSummary\n  - MetaBAT2 (binning)\n  - MaxBin2 (binning)\n  - SemiBin (binning)\n  - DAS Tool (bin merging)\n  - EukRep (eukaryotic classification)\n  - CheckM (bin completeness and contamination)\n  - BUSCO (bin completeness)\n  - GTDB-Tk (bin taxonomic classification)\n\nOther UNLOCK workflows on WorkflowHub: https://workflowhub.eu/projects/16/workflows?view=default<br><br>\n\n**All tool CWL files and other workflows can be found here:**<br>\n  Tools: https://gitlab.com/m-unlock/cwl<br>\n  Workflows: https://gitlab.com/m-unlock/cwl/workflows<br>\n\n**How to setup and use an UNLOCK workflow:**<br>\nhttps://m-unlock.gitlab.io/docs/setup/setup.html<br>\n",
            "outputs": [
                {
                    "label": "Bin files",
                    "doc": "Bins files in fasta format. To be be used in other workflows.",
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputSource": "#main/output_bin_files/bins_out",
                    "id": "#main/bins"
                },
                {
                    "label": "Assembly/Bin read stats",
                    "doc": "General assembly and bin coverage",
                    "type": "File",
                    "outputSource": "#main/bin_readstats/binReadStats",
                    "id": "#main/bins_read_stats"
                },
                {
                    "label": "Bins summary",
                    "doc": "Summary of info about the bins",
                    "type": "File",
                    "outputSource": "#main/bins_summary/bins_summary_table",
                    "id": "#main/bins_summary_table"
                },
                {
                    "label": "BUSCO",
                    "doc": "BUSCO output directory",
                    "type": "Directory",
                    "outputSource": "#main/busco_files_to_folder/results",
                    "id": "#main/busco_output"
                },
                {
                    "label": "CheckM",
                    "doc": "CheckM output directory",
                    "type": "Directory",
                    "outputSource": "#main/checkm_files_to_folder/results",
                    "id": "#main/checkm_output"
                },
                {
                    "label": "DAS Tool",
                    "doc": "DAS Tool output directory",
                    "type": "Directory",
                    "outputSource": "#main/das_tool_files_to_folder/results",
                    "id": "#main/das_tool_output"
                },
                {
                    "label": "EukRep fasta",
                    "doc": "EukRep eukaryotic classified contigs",
                    "type": "File",
                    "outputSource": "#main/eukrep/euk_fasta_out",
                    "id": "#main/eukrep_fasta"
                },
                {
                    "label": "EukRep stats",
                    "doc": "EukRep fasta statistics",
                    "type": "File",
                    "outputSource": "#main/eukrep_stats/output",
                    "id": "#main/eukrep_stats_file"
                },
                {
                    "label": "GTDB-Tk",
                    "doc": "GTDB-Tk output directory",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#main/gtdbtk_files_to_folder/results",
                    "id": "#main/gtdbtk_output"
                },
                {
                    "label": "MaxBin2",
                    "doc": "MaxBin2 output directory",
                    "type": "Directory",
                    "outputSource": "#main/maxbin2_files_to_folder/results",
                    "id": "#main/maxbin2_output"
                },
                {
                    "label": "MetaBAT2",
                    "doc": "MetaBAT2 output directory",
                    "type": "Directory",
                    "outputSource": "#main/metabat2_files_to_folder/results",
                    "id": "#main/metabat2_output"
                },
                {
                    "label": "SemiBin",
                    "doc": "MaxBin2 output directory",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#main/semibin_files_to_folder/results",
                    "id": "#main/semibin_output"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Assembly in fasta format",
                    "label": "Assembly fasta",
                    "loadListing": "no_listing",
                    "id": "#main/assembly"
                },
                {
                    "type": "File",
                    "doc": "Mapping file in sorted bam format containing reads mapped to the assembly",
                    "label": "Bam file",
                    "loadListing": "no_listing",
                    "id": "#main/bam_file"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "BUSCO dataset",
                    "doc": "Directory containing the BUSCO dataset location.",
                    "loadListing": "no_listing",
                    "id": "#main/busco_data"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output destination (not used in the workflow itself)",
                    "doc": "Optional output destination path for cwl-prov reporting.",
                    "id": "#main/destination"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "doc": "Directory containing the GTDB database. When none is given GTDB-Tk will be skipped.",
                    "label": "gtdbtk data directory",
                    "loadListing": "no_listing",
                    "id": "#main/gtdbtk_data"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "Identifier used",
                    "id": "#main/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum memory usage in megabytes",
                    "label": "memory usage (MB)",
                    "default": 4000,
                    "id": "#main/memory"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Run with SemiBin binner",
                    "label": "Run SemiBin",
                    "default": true,
                    "id": "#main/run_semibin"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Semibin Built-in models (human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global/chicken_caecum)",
                    "label": "SemiBin Environment",
                    "default": "global",
                    "id": "#main/semibin_environment"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "CWL base step number",
                    "doc": "Step number for order of steps",
                    "default": 1,
                    "id": "#main/step"
                },
                {
                    "type": "boolean",
                    "label": "Sub workflow Run",
                    "doc": "Use this when you need the output bins as File[] for subsequent analysis workflow steps in another workflow.",
                    "default": false,
                    "id": "#main/sub_workflow"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of threads to use for computational processes",
                    "label": "Threads",
                    "default": 2,
                    "id": "#main/threads"
                }
            ],
            "steps": [
                {
                    "doc": "Depths per bin",
                    "label": "Depths per bin",
                    "run": "#aggregateBinDepths.cwl",
                    "in": [
                        {
                            "source": "#main/das_tool_bins/files",
                            "id": "#main/aggregate_bin_depths/bins"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/aggregate_bin_depths/identifier"
                        },
                        {
                            "source": "#main/metabat2_contig_depths/depths",
                            "id": "#main/aggregate_bin_depths/metabatdepthsFile"
                        }
                    ],
                    "out": [
                        "#main/aggregate_bin_depths/binDepths"
                    ],
                    "id": "#main/aggregate_bin_depths"
                },
                {
                    "doc": "Table general bin and assembly read mapping stats",
                    "label": "Bin and assembly read stats",
                    "run": "#assembly_bins_readstats.cwl",
                    "in": [
                        {
                            "source": "#main/assembly",
                            "id": "#main/bin_readstats/assembly"
                        },
                        {
                            "source": "#main/bam_file",
                            "id": "#main/bin_readstats/bam_file"
                        },
                        {
                            "source": "#main/das_tool/contig2bin",
                            "id": "#main/bin_readstats/binContigs"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/bin_readstats/identifier"
                        }
                    ],
                    "out": [
                        "#main/bin_readstats/binReadStats"
                    ],
                    "id": "#main/bin_readstats"
                },
                {
                    "doc": "Table of all bins and their statistics like size, contigs, completeness etc",
                    "label": "Bins summary",
                    "run": "#bins_summary.cwl",
                    "in": [
                        {
                            "source": "#main/aggregate_bin_depths/binDepths",
                            "id": "#main/bins_summary/bin_depths"
                        },
                        {
                            "source": "#main/das_tool/bin_dir",
                            "id": "#main/bins_summary/bin_dir"
                        },
                        {
                            "source": "#main/busco/batch_summary",
                            "id": "#main/bins_summary/busco_batch"
                        },
                        {
                            "source": "#main/checkm/checkm_out_table",
                            "id": "#main/bins_summary/checkm"
                        },
                        {
                            "source": "#main/gtdbtk/gtdbtk_summary",
                            "id": "#main/bins_summary/gtdbtk"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/bins_summary/identifier"
                        }
                    ],
                    "out": [
                        "#main/bins_summary/bins_summary_table"
                    ],
                    "id": "#main/bins_summary"
                },
                {
                    "doc": "BUSCO assembly completeness workflow",
                    "label": "BUSCO",
                    "run": "#busco.cwl",
                    "when": "$(inputs.bins.length !== 0)",
                    "in": [
                        {
                            "valueFrom": "$(true)",
                            "id": "#main/busco/auto-lineage-prok"
                        },
                        {
                            "source": "#main/das_tool_bins/files",
                            "id": "#main/busco/bins"
                        },
                        {
                            "source": "#main/busco_data",
                            "id": "#main/busco/busco_data"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/busco/identifier"
                        },
                        {
                            "valueFrom": "geno",
                            "id": "#main/busco/mode"
                        },
                        {
                            "source": "#main/remove_unbinned/bins_dir",
                            "id": "#main/busco/sequence_folder"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/busco/threads"
                        }
                    ],
                    "out": [
                        "#main/busco/batch_summary"
                    ],
                    "id": "#main/busco"
                },
                {
                    "doc": "Preparation of BUSCO output files to a specific output folder",
                    "label": "BUSCO output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "BUSCO_Bin_Completeness",
                            "id": "#main/busco_files_to_folder/destination"
                        },
                        {
                            "source": "#main/busco/batch_summary",
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/busco_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#main/busco_files_to_folder/results"
                    ],
                    "id": "#main/busco_files_to_folder"
                },
                {
                    "doc": "CheckM bin quality assessment",
                    "label": "CheckM",
                    "run": "#checkm_lineagewf.cwl",
                    "in": [
                        {
                            "source": "#main/remove_unbinned/bins_dir",
                            "id": "#main/checkm/bin_dir"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/checkm/identifier"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/checkm/threads"
                        }
                    ],
                    "out": [
                        "#main/checkm/checkm_out_table",
                        "#main/checkm/checkm_out_folder"
                    ],
                    "id": "#main/checkm"
                },
                {
                    "doc": "Preparation of CheckM output files to a specific output folder",
                    "label": "CheckM output",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "CheckM_Bin_Quality",
                            "id": "#main/checkm_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/checkm/checkm_out_table"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/checkm_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/checkm/checkm_out_folder"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/checkm_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#main/checkm_files_to_folder/results"
                    ],
                    "id": "#main/checkm_files_to_folder"
                },
                {
                    "doc": "Compress GTDB-Tk output folder",
                    "label": "Compress GTDB-Tk",
                    "when": "$(inputs.gtdbtk_data !== null)",
                    "run": "#compress_directory.cwl",
                    "in": [
                        {
                            "source": "#main/gtdbtk_data",
                            "id": "#main/compress_gtdbtk/gtdbtk_data"
                        },
                        {
                            "source": "#main/gtdbtk/gtdbtk_out_folder",
                            "id": "#main/compress_gtdbtk/indir"
                        }
                    ],
                    "out": [
                        "#main/compress_gtdbtk/outfile"
                    ],
                    "id": "#main/compress_gtdbtk"
                },
                {
                    "doc": "DAS Tool",
                    "label": "DAS Tool integrate predictions from multiple binning tools",
                    "run": "#das_tool.cwl",
                    "in": [
                        {
                            "source": "#main/assembly",
                            "id": "#main/das_tool/assembly"
                        },
                        {
                            "source": [
                                "#main/metabat2_contig2bin/table",
                                "#main/maxbin2_contig2bin/table",
                                "#main/semibin_contig2bin/table"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/das_tool/bin_tables"
                        },
                        {
                            "valueFrom": "${\n  if (inputs.run_semibin) {\n    return \"MetaBAT2,MaxBin2,SemiBin\";\n  } else {\n    return \"MetaBAT2,MaxBin2\";\n  }\n}\n",
                            "id": "#main/das_tool/binner_labels"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/das_tool/identifier"
                        },
                        {
                            "source": "#main/run_semibin",
                            "id": "#main/das_tool/run_semibin"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/das_tool/threads"
                        }
                    ],
                    "out": [
                        "#main/das_tool/bin_dir",
                        "#main/das_tool/summary",
                        "#main/das_tool/contig2bin",
                        "#main/das_tool/log"
                    ],
                    "id": "#main/das_tool"
                },
                {
                    "doc": "DAS Tool bins folder to File array for further analysis",
                    "label": "Bin dir to files[]",
                    "run": "#folder_to_files.cwl",
                    "in": [
                        {
                            "source": "#main/das_tool/bin_dir",
                            "id": "#main/das_tool_bins/folder"
                        }
                    ],
                    "out": [
                        "#main/das_tool_bins/files"
                    ],
                    "id": "#main/das_tool_bins"
                },
                {
                    "doc": "Preparation of DAS Tool output files to a specific output folder.",
                    "label": "DAS Tool output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "Bin_Refinement_DAS_Tool",
                            "id": "#main/das_tool_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/das_tool/log",
                                "#main/das_tool/summary",
                                "#main/das_tool/contig2bin",
                                "#main/aggregate_bin_depths/binDepths"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/das_tool_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/das_tool/bin_dir"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/das_tool_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#main/das_tool_files_to_folder/results"
                    ],
                    "id": "#main/das_tool_files_to_folder"
                },
                {
                    "doc": "EukRep, eukaryotic sequence classification",
                    "label": "EukRep",
                    "run": "#eukrep.cwl",
                    "in": [
                        {
                            "source": "#main/assembly",
                            "id": "#main/eukrep/assembly"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/eukrep/identifier"
                        }
                    ],
                    "out": [
                        "#main/eukrep/euk_fasta_out"
                    ],
                    "id": "#main/eukrep"
                },
                {
                    "doc": "EukRep fasta statistics",
                    "label": "EukRep stats",
                    "run": "#raw_n50.cwl",
                    "in": [
                        {
                            "valueFrom": "$(inputs.tmp_id)_EukRep",
                            "id": "#main/eukrep_stats/identifier"
                        },
                        {
                            "source": "#main/eukrep/euk_fasta_out",
                            "id": "#main/eukrep_stats/input_fasta"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/eukrep_stats/tmp_id"
                        }
                    ],
                    "out": [
                        "#main/eukrep_stats/output"
                    ],
                    "id": "#main/eukrep_stats"
                },
                {
                    "doc": "Taxomic assigment of bins with GTDB-Tk",
                    "label": "GTDBTK",
                    "when": "$(inputs.gtdbtk_data !== null && inputs.bins.length !== 0)",
                    "run": "#gtdbtk_classify_wf.cwl",
                    "in": [
                        {
                            "source": "#main/remove_unbinned/bins_dir",
                            "id": "#main/gtdbtk/bin_dir"
                        },
                        {
                            "source": "#main/das_tool_bins/files",
                            "id": "#main/gtdbtk/bins"
                        },
                        {
                            "source": "#main/gtdbtk_data",
                            "id": "#main/gtdbtk/gtdbtk_data"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/gtdbtk/identifier"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/gtdbtk/threads"
                        }
                    ],
                    "out": [
                        "#main/gtdbtk/gtdbtk_summary",
                        "#main/gtdbtk/gtdbtk_out_folder"
                    ],
                    "id": "#main/gtdbtk"
                },
                {
                    "doc": "Preparation of GTDB-Tk output files to a specific output folder",
                    "label": "GTBD-Tk output folder",
                    "when": "$(inputs.gtdbtk_data !== null)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "GTDB-Tk_Bin_Taxonomy",
                            "id": "#main/gtdbtk_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/gtdbtk/gtdbtk_summary",
                                "#main/compress_gtdbtk/outfile"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/gtdbtk_files_to_folder/files"
                        },
                        {
                            "source": "#main/gtdbtk_data",
                            "id": "#main/gtdbtk_files_to_folder/gtdbtk_data"
                        }
                    ],
                    "out": [
                        "#main/gtdbtk_files_to_folder/results"
                    ],
                    "id": "#main/gtdbtk_files_to_folder"
                },
                {
                    "doc": "Binning procedure using MaxBin2",
                    "label": "MaxBin2 binning",
                    "run": "#maxbin2.cwl",
                    "in": [
                        {
                            "source": "#main/metabat2_contig_depths/depths",
                            "id": "#main/maxbin2/abundances"
                        },
                        {
                            "source": "#main/assembly",
                            "id": "#main/maxbin2/contigs"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/maxbin2/identifier"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/maxbin2/threads"
                        }
                    ],
                    "out": [
                        "#main/maxbin2/bins",
                        "#main/maxbin2/summary",
                        "#main/maxbin2/log"
                    ],
                    "id": "#main/maxbin2"
                },
                {
                    "label": "MaxBin2 to contig to bins",
                    "doc": "List the contigs and their corresponding bin.",
                    "run": "#fasta_to_contig2bin.cwl",
                    "in": [
                        {
                            "source": "#main/maxbin2_to_folder/results",
                            "id": "#main/maxbin2_contig2bin/bin_folder"
                        },
                        {
                            "valueFrom": "MaxBin2",
                            "id": "#main/maxbin2_contig2bin/binner_name"
                        }
                    ],
                    "out": [
                        "#main/maxbin2_contig2bin/table"
                    ],
                    "id": "#main/maxbin2_contig2bin"
                },
                {
                    "doc": "Preparation of maxbin2 output files to a specific output folder.",
                    "label": "MaxBin2 output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "Binner_MaxBin2",
                            "id": "#main/maxbin2_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/maxbin2/summary",
                                "#main/maxbin2/log",
                                "#main/maxbin2_contig2bin/table"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/maxbin2_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/maxbin2_to_folder/results"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/maxbin2_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#main/maxbin2_files_to_folder/results"
                    ],
                    "id": "#main/maxbin2_files_to_folder"
                },
                {
                    "doc": "Create folder with MaxBin2 bins",
                    "label": "MaxBin2 bins to folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "MaxBin2_bins",
                            "id": "#main/maxbin2_to_folder/destination"
                        },
                        {
                            "source": "#main/maxbin2/bins",
                            "id": "#main/maxbin2_to_folder/files"
                        }
                    ],
                    "out": [
                        "#main/maxbin2_to_folder/results"
                    ],
                    "id": "#main/maxbin2_to_folder"
                },
                {
                    "doc": "Binning procedure using MetaBAT2",
                    "label": "MetaBAT2 binning",
                    "run": "#metabat2.cwl",
                    "in": [
                        {
                            "source": "#main/assembly",
                            "id": "#main/metabat2/assembly"
                        },
                        {
                            "source": "#main/metabat2_contig_depths/depths",
                            "id": "#main/metabat2/depths"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/metabat2/identifier"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/metabat2/threads"
                        }
                    ],
                    "out": [
                        "#main/metabat2/bin_dir",
                        "#main/metabat2/log"
                    ],
                    "id": "#main/metabat2"
                },
                {
                    "label": "MetaBAT2 to contig to bins",
                    "doc": "List the contigs and their corresponding bin.",
                    "run": "#fasta_to_contig2bin.cwl",
                    "in": [
                        {
                            "source": "#main/metabat2_filter_bins/output_folder",
                            "id": "#main/metabat2_contig2bin/bin_folder"
                        },
                        {
                            "valueFrom": "MetaBAT2",
                            "id": "#main/metabat2_contig2bin/binner_name"
                        },
                        {
                            "valueFrom": "fa",
                            "id": "#main/metabat2_contig2bin/extension"
                        }
                    ],
                    "out": [
                        "#main/metabat2_contig2bin/table"
                    ],
                    "id": "#main/metabat2_contig2bin"
                },
                {
                    "label": "contig depths",
                    "doc": "MetabatContigDepths to obtain the depth file used in the MetaBat2 and SemiBin binning process",
                    "run": "#metabatContigDepths.cwl",
                    "in": [
                        {
                            "source": "#main/bam_file",
                            "id": "#main/metabat2_contig_depths/bamFile"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/metabat2_contig_depths/identifier"
                        }
                    ],
                    "out": [
                        "#main/metabat2_contig_depths/depths"
                    ],
                    "id": "#main/metabat2_contig_depths"
                },
                {
                    "doc": "Preparation of MetaBAT2 output files + unbinned contigs to a specific output folder",
                    "label": "MetaBAT2 output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "Binner_MetaBAT2",
                            "id": "#main/metabat2_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/metabat2/log",
                                "#main/metabat2_contig_depths/depths",
                                "#main/metabat2_contig2bin/table"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/metabat2_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/metabat2/bin_dir"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/metabat2_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#main/metabat2_files_to_folder/results"
                    ],
                    "id": "#main/metabat2_files_to_folder"
                },
                {
                    "doc": "Only keep genome bin fasta files (exlude e.g TooShort.fa)",
                    "label": "Keep MetaBAT2 genome bins",
                    "run": "#folder_file_regex.cwl",
                    "in": [
                        {
                            "source": "#main/metabat2/bin_dir",
                            "id": "#main/metabat2_filter_bins/folder"
                        },
                        {
                            "valueFrom": "MetaBAT2_bins",
                            "id": "#main/metabat2_filter_bins/output_folder_name"
                        },
                        {
                            "valueFrom": "bin\\.[0-9]+\\.fa",
                            "id": "#main/metabat2_filter_bins/regex"
                        }
                    ],
                    "out": [
                        "#main/metabat2_filter_bins/output_folder"
                    ],
                    "id": "#main/metabat2_filter_bins"
                },
                {
                    "doc": "Bin files for subsequent workflow runs when sub_worflow = true",
                    "label": "Bin files",
                    "when": "$(inputs.sub_workflow)",
                    "run": {
                        "class": "ExpressionTool",
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "inputs": [
                            {
                                "type": {
                                    "type": "array",
                                    "items": "File"
                                },
                                "id": "#main/output_bin_files/run/bins"
                            }
                        ],
                        "outputs": [
                            {
                                "type": {
                                    "type": "array",
                                    "items": "File"
                                },
                                "id": "#main/output_bin_files/run/bins_out"
                            }
                        ],
                        "expression": "${ return inputs.bins; }"
                    },
                    "in": [
                        {
                            "source": "#main/das_tool_bins/files",
                            "id": "#main/output_bin_files/bins"
                        },
                        {
                            "source": "#main/sub_workflow",
                            "id": "#main/output_bin_files/sub_workflow"
                        }
                    ],
                    "out": [
                        "#main/output_bin_files/bins_out"
                    ],
                    "id": "#main/output_bin_files"
                },
                {
                    "doc": "Remove unbinned fasta from bin directory. So analysed by subsequent tools.",
                    "label": "Remove unbinned",
                    "run": {
                        "class": "ExpressionTool",
                        "requirements": [
                            {
                                "class": "InlineJavascriptRequirement"
                            }
                        ],
                        "hints": [
                            {
                                "loadListing": "shallow_listing",
                                "class": "LoadListingRequirement"
                            }
                        ],
                        "inputs": [
                            {
                                "type": "Directory",
                                "id": "#main/remove_unbinned/run/bins"
                            }
                        ],
                        "outputs": [
                            {
                                "type": "Directory",
                                "id": "#main/remove_unbinned/run/bins_dir"
                            }
                        ],
                        "expression": "${  \n  var regex = new RegExp('.*unbinned.*');\n  var array = [];\n  for (var i = 0; i < inputs.bins.listing.length; i++) {\n    if (!regex.test(inputs.bins.listing[i].location)){\n      array = array.concat(inputs.bins.listing[i]);\n    }\n  }\n  var r = {\n    'bins_dir':\n      { \"class\": \"Directory\",\n        \"basename\": \"DAS_Tool_genome_bins\",\n        \"listing\": array\n      }\n    };\n  return r;\n}\n"
                    },
                    "in": [
                        {
                            "source": "#main/das_tool/bin_dir",
                            "id": "#main/remove_unbinned/bins"
                        }
                    ],
                    "out": [
                        "#main/remove_unbinned/bins_dir"
                    ],
                    "id": "#main/remove_unbinned"
                },
                {
                    "doc": "Binning procedure using SemiBin",
                    "label": "Semibin binning",
                    "run": "#semibin_single_easy_bin.cwl",
                    "when": "$(inputs.run_semibin)",
                    "in": [
                        {
                            "source": "#main/assembly",
                            "id": "#main/semibin/assembly"
                        },
                        {
                            "source": "#main/semibin_environment",
                            "id": "#main/semibin/environment"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/semibin/identifier"
                        },
                        {
                            "source": "#main/metabat2_contig_depths/depths",
                            "id": "#main/semibin/metabat2_depth_file"
                        },
                        {
                            "source": "#main/run_semibin",
                            "id": "#main/semibin/run_semibin"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/semibin/threads"
                        }
                    ],
                    "out": [
                        "#main/semibin/recluster_bins",
                        "#main/semibin/data",
                        "#main/semibin/data_split",
                        "#main/semibin/model",
                        "#main/semibin/coverage"
                    ],
                    "id": "#main/semibin"
                },
                {
                    "label": "SemiBin to contig to bins",
                    "doc": "List the contigs and their corresponding bin.",
                    "run": "#fasta_to_contig2bin.cwl",
                    "when": "$(inputs.run_semibin)",
                    "in": [
                        {
                            "source": "#main/semibin/recluster_bins",
                            "id": "#main/semibin_contig2bin/bin_folder"
                        },
                        {
                            "valueFrom": "SemiBin",
                            "id": "#main/semibin_contig2bin/binner_name"
                        },
                        {
                            "valueFrom": "fa",
                            "id": "#main/semibin_contig2bin/extension"
                        },
                        {
                            "source": "#main/run_semibin",
                            "id": "#main/semibin_contig2bin/run_semibin"
                        }
                    ],
                    "out": [
                        "#main/semibin_contig2bin/table"
                    ],
                    "id": "#main/semibin_contig2bin"
                },
                {
                    "doc": "Preparation of SemiBin output files to a specific output folder.",
                    "label": "SemiBin output folder",
                    "run": "#files_to_folder.cwl",
                    "when": "$(inputs.run_semibin)",
                    "in": [
                        {
                            "valueFrom": "Binner_SemiBin",
                            "id": "#main/semibin_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/semibin_contig2bin/table",
                                "#main/semibin/data",
                                "#main/semibin/data_split",
                                "#main/semibin/model",
                                "#main/semibin/coverage"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/semibin_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/semibin/recluster_bins"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/semibin_files_to_folder/folders"
                        },
                        {
                            "source": "#main/run_semibin",
                            "id": "#main/semibin_files_to_folder/run_semibin"
                        }
                    ],
                    "out": [
                        "#main/semibin_files_to_folder/results"
                    ],
                    "id": "#main/semibin_files_to_folder"
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
            "https://schema.org/dateCreated": "2022-00-00",
            "https://schema.org/dateModified": "2023-01-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        }
    ],
    "cwlVersion": "v1.2",
    "$namespaces": {
        "s": "https://schema.org/"
    }
}
