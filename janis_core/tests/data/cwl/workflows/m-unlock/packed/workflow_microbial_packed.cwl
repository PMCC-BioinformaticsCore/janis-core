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
            "label": "Prodigal",
            "doc": "Prokaryotic gene prediction using Prodigal",
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/prodigal:2.6.3",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.6.3"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/prodigal"
                            ],
                            "package": "prodigal"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "prodigal"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.input_fasta.nameroot).prodigal",
                    "prefix": "-o"
                },
                {
                    "valueFrom": "$(inputs.input_fasta.nameroot).prodigal.ffn",
                    "prefix": "-d"
                },
                {
                    "valueFrom": "$(inputs.input_fasta.nameroot).prodigal.faa",
                    "prefix": "-a"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.input_fasta.nameroot).prodigal.faa"
                    },
                    "id": "#prodigal.cwl/predicted_proteins_faa"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.input_fasta.nameroot).prodigal.ffn"
                    },
                    "id": "#prodigal.cwl/predicted_proteins_ffn"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.input_fasta.nameroot).prodigal"
                    },
                    "id": "#prodigal.cwl/predicted_proteins_out"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "-i"
                    },
                    "id": "#prodigal.cwl/input_fasta"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Input is a meta-genome",
                    "inputBinding": {
                        "prefix": "-p",
                        "valueFrom": "meta"
                    },
                    "id": "#prodigal.cwl/meta_mode"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Input is an isolate genome",
                    "inputBinding": {
                        "prefix": "-p",
                        "valueFrom": "single"
                    },
                    "id": "#prodigal.cwl/single_mode"
                }
            ],
            "id": "#prodigal.cwl",
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
                },
                {
                    "class": "https://schema.org/Person",
                    "https://schema.org/identifier": "https://orcid.org/0000-0001-6867-2039",
                    "https://schema.org/name": "Ekaterina Sakharova"
                }
            ],
            "https://schema.org/copyrightHolder'": "EMBL - European Bioinformatics Institute",
            "https://schema.org/license'": "https://www.apache.org/licenses/LICENSE-2.0",
            "https://schema.org/citation": "https://m-unlock.nl",
            "https://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "https://schema.org/dateCreated": "2022-06-00",
            "https://schema.org/dateModified": "2022-08-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential",
            "https://schema.org/copyrightNotice": " Copyright < 2022 EMBL - European Bioinformatics Institute This file has been modified by UNLOCK - Unlocking Microbial Potential "
        },
        {
            "class": "CommandLineTool",
            "label": "spades genomic assembler",
            "doc": "Runs the spades assembler using a dataset file\n",
            "requirements": [
                {
                    "listing": [
                        {
                            "entryname": "input_spades.json",
                            "entry": "[\n  {\n    orientation: \"fr\",\n    type: \"paired-end\",\n    right reads: $( inputs.forward_reads.map( function(x) {return  x.path} ) ),\n    left reads: $( inputs.reverse_reads.map( function(x) {return  x.path} ) )\n  }            \n  ${\n    var pacbio=\"\"\n      if (inputs.pacbio_reads != null) {\n       pacbio+=',{ type: \"pacbio\", single reads: [\"' + inputs.pacbio_reads.map( function(x) {return  x.path} ).join('\",\"') + '\"] }' \n    }\n    return pacbio;\n  }\n  ${\n    var nanopore=\"\"\n      if (inputs.nanopore_reads != null) {\n       nanopore+=',{ type: \"nanopore\", single reads: [\"' + inputs.nanopore_reads.map( function(x) {return  x.path} ).join('\",\"') + '\"] }'\n      //  nanopore+=',{ type: \"nanopore\", single reads: [\"' + inputs.nanopore_reads.join('\",\"') + '\"] }'\n    }\n    return nanopore;\n  }\n]"
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
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/spades:3.15.5",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "3.15.5"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/spades"
                            ],
                            "package": "spades"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "spades.py",
                "--dataset",
                "input_spades.json"
            ],
            "arguments": [
                {
                    "valueFrom": "$(runtime.outdir)/output",
                    "prefix": "-o"
                },
                {
                    "valueFrom": "$(inputs.memory / 1000)",
                    "prefix": "--memory"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is required for IonTorrent data",
                    "label": "iontorrent data",
                    "inputBinding": {
                        "prefix": "--iontorrent"
                    },
                    "id": "#spades.cwl/IonTorrent"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is required for biosyntheticSPAdes mode",
                    "label": "biosynthetic spades mode",
                    "inputBinding": {
                        "prefix": "--bio"
                    },
                    "id": "#spades.cwl/biosyntheticSPAdes"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "The file containing the forward reads",
                    "label": "Forward reads",
                    "id": "#spades.cwl/forward_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is highly recommended for high-coverage isolate and multi-cell data",
                    "label": "high-coverage mode",
                    "inputBinding": {
                        "prefix": "--isolate"
                    },
                    "id": "#spades.cwl/isolate"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory used in megabytes",
                    "id": "#spades.cwl/memory"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is required for metagenomic sample data",
                    "label": "metagenomics sample",
                    "inputBinding": {
                        "prefix": "--meta"
                    },
                    "id": "#spades.cwl/metagenome"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "Fastq file with Oxford NanoPore reads",
                    "label": "NanoPore reads",
                    "id": "#spades.cwl/nanopore_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Runs only assembling (without read error correction)",
                    "label": "Only assembler",
                    "inputBinding": {
                        "prefix": "--only-assembler"
                    },
                    "id": "#spades.cwl/only_assembler"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "Fastq file with PacBio CLR reads",
                    "label": "PacBio CLR reads",
                    "id": "#spades.cwl/pacbio_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "runs plasmidSPAdes pipeline for plasmid detection",
                    "label": "plasmid spades run",
                    "inputBinding": {
                        "prefix": "--plasmid"
                    },
                    "id": "#spades.cwl/plasmid"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "The file containing the reverse reads",
                    "label": "Reverse reads",
                    "id": "#spades.cwl/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is required for RNA-Seq data",
                    "label": "rnaseq data",
                    "inputBinding": {
                        "prefix": "--rna"
                    },
                    "id": "#spades.cwl/rna"
                },
                {
                    "type": "int",
                    "doc": "number of threads to use",
                    "label": "threads",
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#spades.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/assembly_graph.fastg"
                    },
                    "id": "#spades.cwl/assembly_graph"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/contigs.fasta"
                    },
                    "id": "#spades.cwl/contigs"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/contigs.paths"
                    },
                    "id": "#spades.cwl/contigs_assembly_paths"
                },
                {
                    "label": "contigs before repeat resolution",
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/before_rr.fasta"
                    },
                    "id": "#spades.cwl/contigs_before_rr"
                },
                {
                    "label": "internal configuration file",
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/dataset.info"
                    },
                    "id": "#spades.cwl/internal_config"
                },
                {
                    "label": "internal YAML data set file",
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/input_dataset.yaml"
                    },
                    "id": "#spades.cwl/internal_dataset"
                },
                {
                    "label": "SPAdes log",
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/spades.log"
                    },
                    "id": "#spades.cwl/log"
                },
                {
                    "label": "information about SPAdes parameters in this run",
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/params.txt"
                    },
                    "id": "#spades.cwl/params"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/scaffolds.fasta"
                    },
                    "id": "#spades.cwl/scaffolds"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "output/scaffolds.paths"
                    },
                    "id": "#spades.cwl/scaffolds_assembly_paths"
                }
            ],
            "id": "#spades.cwl",
            "http://schema.org/author": [
                {
                    "class": "http://schema.org/Person",
                    "http://schema.org/identifier": "https://orcid.org/0000-0001-8172-8981",
                    "http://schema.org/email": "mailto:jasper.koehorst@wur.nl",
                    "http://schema.org/name": "Jasper Koehorst"
                },
                {
                    "class": "http://schema.org/Person",
                    "http://schema.org/identifier": "https://orcid.org/0000-0001-9524-5964",
                    "http://schema.org/email": "mailto:bart.nijsse@wur.nl",
                    "http://schema.org/name": "Bart Nijsse"
                }
            ],
            "http://schema.org/citation": "https://m-unlock.nl",
            "http://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "http://schema.org/dateCreated": "2020-00-00",
            "http://schema.org/license": "https://spdx.org/licenses/CC0-1.0.html",
            "http://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
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
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "label": "Microbial (meta-)? genome conversion and annotation",
            "doc": "Workflow for the assembly, conversion, gene predication and annotation for microbial systems",
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is required for IonTorrent data",
                    "label": "iontorrent data",
                    "id": "#main/IonTorrent"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is required for biosyntheticSPAdes mode",
                    "label": "biosynthetic spades mode",
                    "id": "#main/biosyntheticSPAdes"
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
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "The file containing the forward reads",
                    "label": "Forward reads",
                    "id": "#main/forward_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is highly recommended for high-coverage isolate and multi-cell data",
                    "label": "high-coverage mode",
                    "id": "#main/isolate"
                },
                {
                    "type": "int",
                    "doc": "Memory used in megabytes",
                    "id": "#main/memory"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is required for metagenomic sample data",
                    "label": "metagenomics sample",
                    "id": "#main/metagenomics"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "file with PacBio reads",
                    "label": "pacbio reads",
                    "id": "#main/pacbio_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "runs plasmidSPAdes pipeline for plasmid detection",
                    "label": "plasmid spades run",
                    "id": "#main/plasmid"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "The file containing the reverse reads",
                    "label": "Reverse reads",
                    "id": "#main/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "this flag is required for RNA-Seq data",
                    "label": "rnaseq data",
                    "id": "#main/rna"
                },
                {
                    "type": "int",
                    "doc": "number of threads to use",
                    "label": "threads",
                    "id": "#main/threads"
                }
            ],
            "steps": [
                {
                    "run": "#prodigal.cwl",
                    "in": [
                        {
                            "source": "#main/spades/scaffolds",
                            "id": "#main/prodigal/input_fasta"
                        },
                        {
                            "source": "#main/metagenomics",
                            "id": "#main/prodigal/meta"
                        },
                        {
                            "source": "#main/isolate",
                            "id": "#main/prodigal/single"
                        }
                    ],
                    "out": [
                        "#main/prodigal/stdout",
                        "#main/prodigal/stderr",
                        "#main/prodigal/predicted_proteins_out",
                        "#main/prodigal/predicted_proteins_ffn",
                        "#main/prodigal/predicted_proteins_faa"
                    ],
                    "id": "#main/prodigal"
                },
                {
                    "$import": "#main/prodigal_files_to_folder"
                },
                {
                    "run": "#spades.cwl",
                    "in": [
                        {
                            "source": "#main/IonTorrent",
                            "id": "#main/spades/IonTorrent"
                        },
                        {
                            "source": "#main/biosyntheticSPAdes",
                            "id": "#main/spades/biosyntheticSPAdes"
                        },
                        {
                            "source": "#main/forward_reads",
                            "id": "#main/spades/forward_reads"
                        },
                        {
                            "source": "#main/isolate",
                            "id": "#main/spades/isolate"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/spades/memory"
                        },
                        {
                            "source": "#main/metagenomics",
                            "id": "#main/spades/metagenomics"
                        },
                        {
                            "source": "#main/pacbio_reads",
                            "id": "#main/spades/pacbio_reads"
                        },
                        {
                            "source": "#main/plasmid",
                            "id": "#main/spades/plasmid"
                        },
                        {
                            "source": "#main/reverse_reads",
                            "id": "#main/spades/reverse_reads"
                        },
                        {
                            "source": "#main/rna",
                            "id": "#main/spades/rna"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/spades/threads"
                        }
                    ],
                    "out": [
                        "#main/spades/stdout",
                        "#main/spades/stderr",
                        "#main/spades/contigs",
                        "#main/spades/scaffolds",
                        "#main/spades/assembly_graph",
                        "#main/spades/contigs_assembly_graph",
                        "#main/spades/scaffolds_assembly_graph",
                        "#main/spades/contigs_before_rr",
                        "#main/spades/params",
                        "#main/spades/log",
                        "#main/spades/internal_config"
                    ],
                    "id": "#main/spades"
                },
                {
                    "$import": "#main/spades_files_to_folder"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputSource": "#main/prodigal_files_to_folder/results",
                    "id": "#main/prodigal_files_to_folder"
                },
                {
                    "type": "Directory",
                    "outputSource": "#main/spades_files_to_folder/results",
                    "id": "#main/spades_files_to_folder"
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
