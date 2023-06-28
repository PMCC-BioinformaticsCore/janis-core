{
    "$graph": [
        {
            "class": "CommandLineTool",
            "label": "kallisto indexer",
            "doc": "Creates the index for pseudoalignment tool kallisto\nhttps://github.com/common-workflow-library/bio-cwl-tools/tree/release/Kallisto",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "kallisto",
                            "writable": true
                        }
                    ]
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/kallisto:0.48.0",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "0.48.0"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/kallisto"
                            ],
                            "package": "kallisto"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 200
                    },
                    "id": "#kallisto_index.cwl/inputFile"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--kmer-size=",
                        "separate": false
                    },
                    "id": "#kallisto_index.cwl/kmerSize"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "prefix": "--make-unique"
                    },
                    "id": "#kallisto_index.cwl/makeUnique"
                }
            ],
            "arguments": [
                {
                    "prefix": "--index=",
                    "separate": false,
                    "valueFrom": "$(\"kallisto/\" + inputs.inputFile.nameroot)_kallisto.idx"
                }
            ],
            "baseCommand": [
                "kallisto",
                "index"
            ],
            "id": "#kallisto_index.cwl",
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
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "kallisto"
                    },
                    "id": "#kallisto_index.cwl/kallisto_indexFolder"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "label": "STAR indexer",
            "doc": "Creates the Genome index for STAR spliced RNAseq aligner (single fasta)\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "STAR",
                            "writable": true
                        }
                    ]
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/star:2.7.10a",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.7.10a"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/star"
                            ],
                            "package": "star"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--genomeChrBinNbits"
                    },
                    "id": "#star_index.cwl/genomeChrBinNbits"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--genomeSAindexNbases"
                    },
                    "id": "#star_index.cwl/genomeSAindexNbases"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--genomeFastaFiles"
                    },
                    "id": "#star_index.cwl/inputFile"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "--sjdbFileChrStartEnd"
                    },
                    "id": "#star_index.cwl/sjdbFileChrStartEnd"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "--sjdbGTFfile"
                    },
                    "id": "#star_index.cwl/sjdbGTFfile"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--sjdbOverhang"
                    },
                    "id": "#star_index.cwl/sjdbOverhang"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#star_index.cwl/threads"
                }
            ],
            "baseCommand": [
                "STAR",
                "--runMode",
                "genomeGenerate"
            ],
            "arguments": [
                "--genomeDir",
                "STAR"
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "STAR"
                    },
                    "id": "#star_index.cwl/STAR"
                }
            ],
            "id": "#star_index.cwl",
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
            "label": "Bowtie2 indexer",
            "doc": "Creates the Genome index for bowtie2 (single fasta)\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "bowtie2",
                            "writable": true
                        }
                    ]
                }
            ],
            "hints": [
                {
                    "packages": [
                        {
                            "version": [
                                "2.4.5"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/bowtie2"
                            ],
                            "package": "bowtie2"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "label": "Reference file",
                    "doc": "Reference file in fasta format",
                    "type": "File",
                    "inputBinding": {
                        "prefix": "-f"
                    },
                    "id": "#bowtie2_index.cwl/reference"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#bowtie2_index.cwl/threads"
                }
            ],
            "arguments": [
                {
                    "valueFrom": "$(\"bowtie2/\" + inputs.reference.nameroot)",
                    "position": 100
                },
                {
                    "valueFrom": "$(runtime.cores)",
                    "prefix": "--threads"
                }
            ],
            "baseCommand": [
                "/unlock/infrastructure/binaries/bowtie2/bowtie2-v2.4.5/bowtie2-build"
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "bowtie2"
                    },
                    "id": "#bowtie2_index.cwl/bowtie2_index_dir"
                }
            ],
            "id": "#bowtie2_index.cwl",
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
            "https://schema.org/dateModified": "2022-02-23",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "SAPP conversion RDF2FASTA",
            "doc": "SAPP conversion tool utilizing the function RDF2FASTA\n",
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
                    "inputBinding": {
                        "prefix": "-i"
                    },
                    "id": "#conversion_rdf2fasta.cwl/inputFile"
                }
            ],
            "arguments": [
                {
                    "prefix": "-transcript",
                    "valueFrom": "$(inputs.inputFile.nameroot)_transcripts.fasta"
                },
                {
                    "prefix": "-protein",
                    "valueFrom": "$(inputs.inputFile.nameroot)_proteins.fasta"
                }
            ],
            "baseCommand": [
                "java",
                "-Xmx5G",
                "-jar",
                "/SAPP-2.0.jar",
                "-rdf2fasta"
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.inputFile.nameroot)_proteins.fasta"
                    },
                    "id": "#conversion_rdf2fasta.cwl/proteins"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.inputFile.nameroot)_transcripts.fasta"
                    },
                    "id": "#conversion_rdf2fasta.cwl/transcripts"
                }
            ],
            "id": "#conversion_rdf2fasta.cwl",
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
            "https://schema.org/dateCreated": "2020-08-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "SAPP conversion RDF2GTF",
            "doc": "SAPP conversion tool utilizing the function RDF2GTF\n",
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
                    "inputBinding": {
                        "prefix": "-input"
                    },
                    "id": "#conversion_rdf2gtf.cwl/inputFile"
                }
            ],
            "arguments": [
                {
                    "prefix": "-output",
                    "valueFrom": "$(inputs.inputFile.nameroot)"
                }
            ],
            "baseCommand": [
                "java",
                "-Xmx5G",
                "-jar",
                "/SAPP-2.0.jar",
                "-rdf2gtf"
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.inputFile.nameroot).fasta"
                    },
                    "id": "#conversion_rdf2gtf.cwl/genomefasta"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.inputFile.nameroot).gtf"
                    },
                    "id": "#conversion_rdf2gtf.cwl/gtf"
                }
            ],
            "id": "#conversion_rdf2gtf.cwl",
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
            "https://schema.org/dateCreated": "2020-08-00",
            "https://schema.org/dateModified": "2022-05-00",
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
                }
            ],
            "label": "Indices builder from GBOL RDF (TTL)",
            "doc": "Workflow to build different indices for different tools from a genome and transcriptome. \n\nThis workflow expects an (annotated) genome in GBOL ttl format.\n\nSteps:\n  - SAPP: rdf2gtf (genome fasta)\n  - SAPP: rdf2fasta (transcripts fasta)\n  - STAR index (Optional for Eukaryotic origin)\n  - bowtie2 index\n  - kallisto index\n",
            "outputs": [
                {
                    "type": "Directory",
                    "label": "STAR",
                    "doc": "STAR index folder",
                    "outputSource": "#main/STAR/STAR",
                    "id": "#main/STAR"
                },
                {
                    "type": "Directory",
                    "label": "bowtie2",
                    "doc": "bowtie2 index folder",
                    "outputSource": "#main/bowtie2/bowtie2",
                    "id": "#main/bowtie2"
                },
                {
                    "label": "Genome fasta",
                    "doc": "Genome fasta file",
                    "type": [
                        "File"
                    ],
                    "outputSource": [
                        "#main/rdf2gtf/genomefasta"
                    ],
                    "id": "#main/genomefasta"
                },
                {
                    "label": "GTF",
                    "doc": "Genes in GTF format",
                    "type": [
                        "File"
                    ],
                    "outputSource": [
                        "#main/rdf2gtf/gtf"
                    ],
                    "id": "#main/gtf"
                },
                {
                    "label": "kallisto",
                    "doc": "kallisto index folder",
                    "type": "Directory",
                    "outputSource": "#main/kallisto/kallisto_indexFolder",
                    "id": "#main/kallisto"
                },
                {
                    "label": "Proteins",
                    "doc": "Proteins fasta file",
                    "type": [
                        "File"
                    ],
                    "outputSource": [
                        "#main/rdf2fasta/proteins"
                    ],
                    "id": "#main/proteins"
                },
                {
                    "label": "Transcripts",
                    "doc": "Transcripts fasta file",
                    "type": [
                        "File"
                    ],
                    "outputSource": [
                        "#main/rdf2fasta/transcripts"
                    ],
                    "id": "#main/transcripts"
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
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "For small genomes, the parameter --genomeSAindexNbases must be scaled down.",
                    "label": "STAR parameter",
                    "id": "#main/genomeSAindexNbases"
                },
                {
                    "label": "Input File",
                    "doc": "Annotated genome in GBOL turtle file (.ttl) format",
                    "type": "File",
                    "id": "#main/inputFile"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "maximum memory usage in megabytes",
                    "label": "maximum memory usage in megabytes",
                    "default": 4000,
                    "id": "#main/memory"
                },
                {
                    "label": "Run STAR",
                    "doc": "create STAR index for genome if true. (For genomes with exons)",
                    "type": "boolean",
                    "default": false,
                    "id": "#main/run_star"
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
                    "$import": "#main/STAR"
                },
                {
                    "$import": "#main/bowtie2"
                },
                {
                    "$import": "#main/kallisto"
                },
                {
                    "label": "RDF to Fasta",
                    "doc": "Convert input RDF (turtle) file to Genome fasta file.",
                    "run": "#conversion_rdf2fasta.cwl",
                    "in": [
                        {
                            "source": "#main/inputFile",
                            "id": "#main/rdf2fasta/inputFile"
                        }
                    ],
                    "out": [
                        "#main/rdf2fasta/transcripts",
                        "#main/rdf2fasta/proteins"
                    ],
                    "id": "#main/rdf2fasta"
                },
                {
                    "label": "RDF to GTF",
                    "doc": "Convert input RDF (turtle) file to GTF file",
                    "run": "#conversion_rdf2gtf.cwl",
                    "in": [
                        {
                            "source": "#main/inputFile",
                            "id": "#main/rdf2gtf/inputFile"
                        }
                    ],
                    "out": [
                        "#main/rdf2gtf/genomefasta",
                        "#main/rdf2gtf/gtf"
                    ],
                    "id": "#main/rdf2gtf"
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
