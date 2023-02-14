{
    "$graph": [
        {
            "class": "CommandLineTool",
            "label": "Concatenate multiple files",
            "baseCommand": [
                "cat"
            ],
            "stdout": "$(inputs.outname)",
            "hints": [
                {
                    "dockerPull": "debian:buster",
                    "class": "DockerRequirement"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#concatenate.cwl/infiles"
                },
                {
                    "type": "string",
                    "id": "#concatenate.cwl/outname"
                }
            ],
            "id": "#concatenate.cwl",
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
                        "glob": "$(inputs.outname)"
                    },
                    "id": "#concatenate.cwl/output"
                }
            ]
        },
        {
            "class": "CommandLineTool",
            "label": "Prepare fasta DB",
            "doc": "Prepares fasta file for so it does not contain duplicate fasta headers.\nOnly looks at the first part of the header before any whitespace.\nAdds and incremental number in the header.\n\nExpects fasta file(s) or plaintext fasta(s). Not mixed!    \n",
            "requirements": [
                {
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "prepare_fasta_db",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\necho -e \"\\\n#/usr/bin/python3\nimport sys\\n\\\nheaders = set()\\n\\\nc = 0\\n\\\nfor line in sys.stdin:\\n\\\n  splitline = line.split()\\n\\\n  if line[0] == '>':    \\n\\\n    if splitline[0] in headers:\\n\\\n      c += 1\\n\\\n      print(splitline[0]+'.x'+str(c)+' '+' '.join(splitline[1:]))\\n\\\n    else:\\n\\\n      print(line.strip())\\n\\\n    headers.add(splitline[0])\\n\\\n  else:\\n\\\n    print(line.strip())\" > ./dup.py\nout_name=$1\nshift\n\nif file $@ | grep gzip; then\n  zcat $@ | python3 ./dup.py | gzip > $out_name\nelse\n  cat $@ | python3 ./dup.py | gzip > $out_name\nfi"
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
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/python:3.10.6",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
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
                "bash",
                "script.sh"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "fasta files",
                    "doc": "Fasta file(s) to be the prepared. Can also be gzipped (not mixe)",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#prepare_fasta_db.cwl/fasta_files"
                },
                {
                    "type": "string",
                    "label": "Output outfile",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#prepare_fasta_db.cwl/output_file_name"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.output_file_name)"
                    },
                    "id": "#prepare_fasta_db.cwl/fasta_db"
                }
            ],
            "id": "#prepare_fasta_db.cwl",
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
            "https://schema.org/dateCreated": "2022-07-00",
            "https://schema.org/dateModified": "2023-01-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "label": "Convert an array of 1 file to a file object",
            "doc": "Converts the array and returns the first file in the array. \nShould only be used when 1 file is in the array.\n",
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "id": "#array_to_file.cwl/files"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "id": "#array_to_file.cwl/file"
                }
            ],
            "expression": "${\n  var first_file = inputs.files[0];\n  return {'file': first_file}\n}",
            "id": "#array_to_file.cwl"
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
            "baseCommand": [
                "fastqc"
            ],
            "label": "FASTQC",
            "doc": "Performs quality control on FASTQ files\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "FASTQC",
                            "writable": true
                        }
                    ]
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/fastqc:0.11.9",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "0.11.9"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/fastqc"
                            ],
                            "package": "fastp"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "arguments": [
                "--outdir",
                "FASTQC"
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "FastQ file list",
                    "label": "FASTQ file list",
                    "inputBinding": {
                        "position": 100
                    },
                    "id": "#fastqc.cwl/fastq"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "doc": "FastQ files list",
                    "label": "FASTQ files list",
                    "inputBinding": {
                        "position": 101,
                        "prefix": "--nano"
                    },
                    "id": "#fastqc.cwl/nanopore_reads"
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
                    "id": "#fastqc.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "FASTQC/*.html"
                    },
                    "id": "#fastqc.cwl/html_files"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "FASTQC/*.zip"
                    },
                    "id": "#fastqc.cwl/zip_files"
                }
            ],
            "id": "#fastqc.cwl",
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
            "https://schema.org/dateCreated": "2021-11-26",
            "https://schema.org/dateModified": "2022-04-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "Filtlong",
            "doc": "Filtlong is a tool for filtering long reads by quality. It can take a set of long reads and produce a smaller, better subset. \nIt uses both read length (longer is better) and read identity (higher is better) when choosing which reads pass the filter.\n",
            "hints": [
                {
                    "dockerPull": "quay.io/biocontainers/filtlong:0.2.1--hd03093a_1",
                    "class": "DockerRequirement"
                }
            ],
            "requirements": [
                {
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "filtlong",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\noutname=$1\nlongreads=$2\nshift;shift;\nfiltlong $longreads $@ 2> >(tee -a $outname.filtlong.log>&2) | gzip > $outname.fastq.gz"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "baseCommand": [
                "bash",
                "-x",
                "script.sh"
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.output_filename).filtlong.log"
                    },
                    "id": "#filtlong.cwl/log"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.output_filename).fastq.gz"
                    },
                    "id": "#filtlong.cwl/output_reads"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Reference assembly",
                    "doc": "Reference assembly in FASTA format",
                    "inputBinding": {
                        "prefix": "--assembly",
                        "position": 13
                    },
                    "id": "#filtlong.cwl/assembly"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Forward reads",
                    "doc": "Forward reference Illumina reads in FASTQ format",
                    "inputBinding": {
                        "prefix": "-illumina_1",
                        "position": 11
                    },
                    "id": "#filtlong.cwl/forward_reads"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Keep percentage",
                    "doc": "Keep only this percentage of the best reads (measured by bases)",
                    "inputBinding": {
                        "prefix": "--keep_percent",
                        "position": 4
                    },
                    "id": "#filtlong.cwl/keep_percent"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Length weight",
                    "doc": "Weight given to the length score (default; 1)",
                    "inputBinding": {
                        "prefix": "--length_weight",
                        "position": 14
                    },
                    "id": "#filtlong.cwl/length_weight"
                },
                {
                    "type": "File",
                    "label": "Long reads",
                    "doc": "Long reads in fastq format",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#filtlong.cwl/long_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Maximum length",
                    "doc": "Maximum read length threshold",
                    "inputBinding": {
                        "prefix": "--max_length",
                        "position": 6
                    },
                    "id": "#filtlong.cwl/maximum_length"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Mean quality weight",
                    "doc": "Weight given to the mean quality score (default; 1)",
                    "inputBinding": {
                        "prefix": "--mean_q_weight",
                        "position": 15
                    },
                    "id": "#filtlong.cwl/mean_q_weight"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Minimum mean quality",
                    "doc": "Minimum mean quality threshold",
                    "inputBinding": {
                        "prefix": "--min_mean_q",
                        "position": 7
                    },
                    "id": "#filtlong.cwl/min_mean_q"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Minimum window quality",
                    "doc": "Minimum window quality threshold",
                    "inputBinding": {
                        "prefix": "--min_window_q",
                        "position": 8
                    },
                    "id": "#filtlong.cwl/min_window_q"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Minimum length",
                    "doc": "Minimum read length threshold",
                    "inputBinding": {
                        "prefix": "--min_length",
                        "position": 5
                    },
                    "id": "#filtlong.cwl/minimum_length"
                },
                {
                    "type": "string",
                    "label": "Output filename",
                    "doc": "Output filename (fastq.gz will be added by default)",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#filtlong.cwl/output_filename"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Reverse reads",
                    "doc": "Reverse reference Illumina reads in FASTQ format",
                    "inputBinding": {
                        "prefix": "-illumina_2",
                        "position": 12
                    },
                    "id": "#filtlong.cwl/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Split",
                    "doc": "Split reads at this many (or more) consecutive non-k-mer-matching bases",
                    "inputBinding": {
                        "prefix": "--trim",
                        "position": 10
                    },
                    "id": "#filtlong.cwl/split"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Target bases",
                    "doc": "Keep only the best reads up to this many total bases",
                    "inputBinding": {
                        "prefix": "--target_bases",
                        "position": 3
                    },
                    "id": "#filtlong.cwl/target_bases"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Trim",
                    "doc": "Trim non-k-mer-matching bases from start/end of reads",
                    "inputBinding": {
                        "prefix": "--trim",
                        "position": 9
                    },
                    "id": "#filtlong.cwl/trim"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Mean window weight",
                    "doc": "Weight given to the window quality score (default; 1)",
                    "inputBinding": {
                        "prefix": "--window_q_weight",
                        "position": 16
                    },
                    "id": "#filtlong.cwl/window_q_weight"
                }
            ],
            "id": "#filtlong.cwl",
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
            "https://schema.org/dateCreated": "2023-01-03",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "kraken2"
            ],
            "label": "Kraken2",
            "doc": "Kraken2 metagenomics taxomic read classification.\n\nUpdated databases available at: https://benlangmead.github.io/aws-indexes/k2 (e.g. PlusPF-8)\nOriginal db: https://ccb.jhu.edu/software/kraken2/index.shtml?t=downloads\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/kraken2:2.1.2",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.1.2"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/kraken2"
                            ],
                            "package": "kraken2"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.identifier)_$(inputs.database.path.split( '/' ).pop())_kraken2.txt",
                    "prefix": "--output"
                },
                {
                    "valueFrom": "$(inputs.identifier)_$(inputs.database.path.split( '/' ).pop())_kraken2_report.txt",
                    "prefix": "--report"
                },
                "--report-zero-counts",
                "--use-names"
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "doc": "input data is gzip compressed",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--bzip2-compressed"
                    },
                    "default": false,
                    "id": "#kraken2.cwl/bzip2"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Confidence",
                    "doc": "Confidence score threshold (default 0.0) must be in [0, 1]",
                    "inputBinding": {
                        "position": 4,
                        "prefix": "--confidence"
                    },
                    "id": "#kraken2.cwl/confidence"
                },
                {
                    "type": "Directory",
                    "label": "Database",
                    "doc": "Database location of kraken2",
                    "inputBinding": {
                        "prefix": "--db"
                    },
                    "id": "#kraken2.cwl/database"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Forward reads",
                    "doc": "Illumina forward read file",
                    "inputBinding": {
                        "position": 100
                    },
                    "id": "#kraken2.cwl/forward_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "input data is gzip compressed",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "--gzip-compressed"
                    },
                    "default": false,
                    "id": "#kraken2.cwl/gzip"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#kraken2.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Nanopore reads",
                    "doc": "Oxford Nanopore Technologies reads in FASTQ",
                    "inputBinding": {
                        "position": 102
                    },
                    "id": "#kraken2.cwl/nanopore_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Paired end",
                    "doc": "Data is paired end (separate files)",
                    "inputBinding": {
                        "position": 2,
                        "prefix": "--paired"
                    },
                    "default": false,
                    "id": "#kraken2.cwl/paired_end"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Reverse reads",
                    "doc": "Illumina reverse read file",
                    "inputBinding": {
                        "position": 101
                    },
                    "id": "#kraken2.cwl/reverse_reads"
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
                    "id": "#kraken2.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_$(inputs.database.path.split( '/' ).pop())_kraken2_report.txt"
                    },
                    "id": "#kraken2.cwl/sample_report"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_$(inputs.database.path.split( '/' ).pop())_kraken2.txt"
                    },
                    "id": "#kraken2.cwl/standard_report"
                }
            ],
            "id": "#kraken2.cwl",
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
            "https://schema.org/dateCreated": "2021-11-25",
            "https://schema.org/dateModified": "2021-11-04",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/krona:2.8.1",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.8.1"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/krona"
                            ],
                            "package": "krona"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "ktImportTaxonomy"
            ],
            "label": "Krona",
            "doc": "Visualization of Kraken2 report results.\nktImportText -o $1 $2\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "krona_output",
                            "writable": true
                        }
                    ]
                }
            ],
            "arguments": [
                {
                    "prefix": "-o",
                    "valueFrom": "krona_output/$(inputs.kraken.nameroot)_krona.html"
                }
            ],
            "inputs": [
                {
                    "type": "int",
                    "label": "Counts column",
                    "doc": "Column number for count information (default for kraken)",
                    "default": 3,
                    "inputBinding": {
                        "position": 2,
                        "prefix": "-m"
                    },
                    "id": "#krona.cwl/counts"
                },
                {
                    "type": "File",
                    "label": "Tab-delimited text file",
                    "inputBinding": {
                        "position": 10
                    },
                    "id": "#krona.cwl/kraken"
                },
                {
                    "type": "int",
                    "label": "Taxon column",
                    "doc": "Column number for taxon information (default for kraken)",
                    "default": 5,
                    "inputBinding": {
                        "position": 1,
                        "prefix": "-t"
                    },
                    "id": "#krona.cwl/taxonomy"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "krona_output/$(inputs.kraken.nameroot)_krona.html"
                    },
                    "id": "#krona.cwl/krona_html"
                }
            ],
            "id": "#krona.cwl",
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
            "https://schema.org/dateCreated": "2021-12-10",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "Minimap2 to (un)mapped long reads",
            "doc": "Get unmapped or mapped long reads reads in fastq.gz format using minimap2 and samtools. Mainly used for contamination removal.\n - requires pigz!\nminimap2 | samtools | pigz\n",
            "requirements": [
                {
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "minimap_run",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nminimap2 -a -t $5 -x $4 $2 $3 | samtools fastq -@ $5 -n $1 4 | pigz -p $5 > $6_filtered.fastq.gz"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "baseCommand": [
                "bash",
                "-x",
                "script.sh"
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/minimap2:2.24",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.24"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/minimap2"
                            ],
                            "package": "minimap2"
                        },
                        {
                            "version": [
                                "2.6"
                            ],
                            "specs": [
                                "https://anaconda.org/conda-forge/pigz"
                            ],
                            "package": "pigz"
                        },
                        {
                            "version": [
                                "1.15.1"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/samtools"
                            ],
                            "package": "samtools"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "arguments": [
                "${\n  if (inputs.output_mapped){\n    return '-F';\n  } else {\n    return '-f';\n  }\n}\n"
            ],
            "inputs": [
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier",
                    "inputBinding": {
                        "position": 5
                    },
                    "id": "#minimap2_to_fastq.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Keep reads mapped to the reference (default = output unmapped)",
                    "label": "Keep mapped",
                    "default": false,
                    "inputBinding": {
                        "position": 6
                    },
                    "id": "#minimap2_to_fastq.cwl/output_mapped"
                },
                {
                    "type": "string",
                    "doc": "- map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping\n- map-hifi - PacBio HiFi reads vs reference mapping\n- ava-pb/ava-ont - PacBio/Nanopore read overlap\n- asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence\n- splice/splice:hq - long-read/Pacbio-CCS spliced alignment\n- sr - genomic short-read mapping\n",
                    "label": "Read type",
                    "inputBinding": {
                        "position": 3
                    },
                    "id": "#minimap2_to_fastq.cwl/preset"
                },
                {
                    "type": "File",
                    "doc": "Query sequence in FASTQ/FASTA format (can be gzipped).",
                    "label": "Reads",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#minimap2_to_fastq.cwl/reads"
                },
                {
                    "type": "File",
                    "doc": "Target sequence in FASTQ/FASTA format (can be gzipped).",
                    "label": "Reference",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#minimap2_to_fastq.cwl/reference"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of CPU-threads used by minimap2.",
                    "label": "Threads",
                    "default": 4,
                    "inputBinding": {
                        "position": 4
                    },
                    "id": "#minimap2_to_fastq.cwl/threads"
                }
            ],
            "stderr": "$(inputs.identifier)_minimap2.log",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_filtered.fastq.gz"
                    },
                    "id": "#minimap2_to_fastq.cwl/fastq"
                },
                {
                    "type": "File",
                    "id": "#minimap2_to_fastq.cwl/log",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_minimap2.log"
                    }
                }
            ],
            "id": "#minimap2_to_fastq.cwl",
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
            "https://schema.org/dateCreated": "2022-03-00",
            "https://schema.org/dateModified": "2022-04-00",
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
            "label": "Long Read Quality Control and Filtering",
            "doc": "**Workflow for long read quality control and contamination filtering.**\n- FastQC before filtering (read quality control)\n- Filtlong filter on quality and length\n- Kraken2 taxonomic read classification\n- Minimap2 read filtering based on given references\n- FastQC after filtering (read quality control)<br><br>\n\nOther UNLOCK workflows on WorkflowHub: https://workflowhub.eu/projects/16/workflows?view=default<br><br>\n\n**All tool CWL files and other workflows can be found here:**<br>\n  Tools: https://gitlab.com/m-unlock/cwl<br>\n  Workflows: https://gitlab.com/m-unlock/cwl/workflows<br>\n\n**How to setup and use an UNLOCK workflow:**<br>\nhttps://m-unlock.gitlab.io/docs/setup/setup.html<br>\n",
            "outputs": [
                {
                    "type": "File",
                    "label": "Filtered long reads",
                    "doc": "Filtered long reads",
                    "outputSource": [
                        "#main/reference_filter_longreads/fastq",
                        "#main/filtlong/output_reads"
                    ],
                    "pickValue": "first_non_null",
                    "id": "#main/filtered_reads"
                },
                {
                    "type": "Directory",
                    "label": "Filtering reports folder",
                    "doc": "Folder containing all reports of filtering and quality control",
                    "outputSource": "#main/reports_files_to_folder/results",
                    "id": "#main/reports_folder"
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
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "Contamination references fasta file for contamination filtering. Gzipped or not (Not mixed)",
                    "label": "Contamination reference file",
                    "id": "#main/filter_references"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#main/identifier"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Maximum read length threshold (default 90)",
                    "label": "Maximum read length threshold",
                    "default": 90,
                    "id": "#main/keep_percent"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Keep with reads mapped to the given reference (default false)",
                    "label": "Keep mapped reads",
                    "default": false,
                    "id": "#main/keep_reference_mapped_reads"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Kraken2 confidence threshold",
                    "doc": "Confidence score threshold (default 0.0) must be between [0, 1]",
                    "id": "#main/kraken2_confidence"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "Directory"
                        }
                    ],
                    "label": "Kraken2 database",
                    "doc": "Kraken2 database location, multiple databases is possible",
                    "default": [],
                    "id": "#main/kraken2_database"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Weight given to the length score (default 10)",
                    "label": "Length weigth",
                    "default": 10,
                    "id": "#main/length_weight"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Long read sequence file locally fastq format",
                    "label": "Long reads",
                    "id": "#main/longreads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum memory usage in megabytes",
                    "label": "Maximum memory in MB",
                    "default": 4000,
                    "id": "#main/memory"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Minimum read length threshold (default 1000)",
                    "label": "Minimum read length",
                    "default": 1000,
                    "id": "#main/minimum_length"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Prepare reference with unique headers (default true)",
                    "label": "Prepare references",
                    "default": true,
                    "id": "#main/prepare_reference"
                },
                {
                    "type": "string",
                    "doc": "Type of read i.e. PacBio or Nanopore. Used for naming output files.",
                    "label": "Read type",
                    "id": "#main/readtype"
                },
                {
                    "type": "boolean",
                    "doc": "Skip FastQC analyses of raw input data (default; false)",
                    "label": "Skip FastQC before",
                    "default": false,
                    "id": "#main/skip_fastqc_before"
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
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of threads to use for computational processes",
                    "label": "Number of threads",
                    "default": 2,
                    "id": "#main/threads"
                }
            ],
            "steps": [
                {
                    "label": "FastQC after",
                    "doc": "Quality assessment and report of reads before filter",
                    "run": "#fastqc.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/reference_filter_longreads/fastq",
                                "#main/filtlong/output_reads"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/fastqc_longreads_after/nanopore_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/fastqc_longreads_after/threads"
                        }
                    ],
                    "out": [
                        "#main/fastqc_longreads_after/html_files",
                        "#main/fastqc_longreads_after/zip_files"
                    ],
                    "id": "#main/fastqc_longreads_after"
                },
                {
                    "label": "FastQC before",
                    "doc": "Quality assessment and report of reads before filter",
                    "run": "#fastqc.cwl",
                    "when": "$(inputs.skip_fastqc_before == false)",
                    "in": [
                        {
                            "source": [
                                "#main/merge_longreads_fastq/output",
                                "#main/longreads_array_to_file/file"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/fastqc_longreads_before/nanopore_reads"
                        },
                        {
                            "source": "#main/skip_fastqc_before",
                            "id": "#main/fastqc_longreads_before/skip_fastqc_before"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/fastqc_longreads_before/threads"
                        }
                    ],
                    "out": [
                        "#main/fastqc_longreads_before/html_files",
                        "#main/fastqc_longreads_before/zip_files"
                    ],
                    "id": "#main/fastqc_longreads_before"
                },
                {
                    "label": "Filtlong",
                    "doc": "Filter longreads on quality and length",
                    "run": "#filtlong.cwl",
                    "in": [
                        {
                            "source": "#main/identifier",
                            "id": "#main/filtlong/identifier"
                        },
                        {
                            "source": "#main/keep_percent",
                            "id": "#main/filtlong/keep_percent"
                        },
                        {
                            "source": "#main/length_weight",
                            "id": "#main/filtlong/length_weight"
                        },
                        {
                            "source": [
                                "#main/merge_longreads_fastq/output",
                                "#main/longreads_array_to_file/file"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/filtlong/long_reads"
                        },
                        {
                            "source": "#main/minimum_length",
                            "id": "#main/filtlong/minimum_length"
                        },
                        {
                            "valueFrom": "$(inputs.identifier)_$(inputs.readtype)_filtered",
                            "id": "#main/filtlong/output_filename"
                        },
                        {
                            "source": "#main/readtype",
                            "id": "#main/filtlong/readtype"
                        }
                    ],
                    "out": [
                        "#main/filtlong/output_reads",
                        "#main/filtlong/log"
                    ],
                    "id": "#main/filtlong"
                },
                {
                    "label": "Array to file",
                    "doc": "Pick first file of longreads when only 1 file is given",
                    "when": "$(inputs.longreads.length === 1)",
                    "run": "#array_to_file.cwl",
                    "in": [
                        {
                            "source": "#main/longreads",
                            "id": "#main/longreads_array_to_file/files"
                        },
                        {
                            "source": "#main/longreads",
                            "id": "#main/longreads_array_to_file/longreads"
                        }
                    ],
                    "out": [
                        "#main/longreads_array_to_file/file"
                    ],
                    "id": "#main/longreads_array_to_file"
                },
                {
                    "label": "Kraken2",
                    "doc": "Taxonomic classification of FASTQ reads",
                    "when": "$(inputs.kraken2_database !== null && inputs.kraken2_database.length !== 0)",
                    "run": "#kraken2.cwl",
                    "scatter": "#main/longreads_quality_kraken2/database",
                    "in": [
                        {
                            "source": "#main/kraken2_confidence",
                            "id": "#main/longreads_quality_kraken2/confidence"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/longreads_quality_kraken2/database"
                        },
                        {
                            "source": "#main/identifier",
                            "valueFrom": "$(self+\"_\"+inputs.readtype)_unfiltered",
                            "id": "#main/longreads_quality_kraken2/identifier"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/longreads_quality_kraken2/kraken2_database"
                        },
                        {
                            "source": [
                                "#main/merge_longreads_fastq/output",
                                "#main/longreads_array_to_file/file"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/longreads_quality_kraken2/nanopore_reads"
                        },
                        {
                            "source": "#main/readtype",
                            "id": "#main/longreads_quality_kraken2/readtype"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/longreads_quality_kraken2/threads"
                        }
                    ],
                    "out": [
                        "#main/longreads_quality_kraken2/sample_report"
                    ],
                    "id": "#main/longreads_quality_kraken2"
                },
                {
                    "label": "Krona",
                    "doc": "Visualization of Kraken2 classification with Krona",
                    "when": "$(inputs.kraken2_database !== null && inputs.kraken2_database.length !== 0)",
                    "run": "#krona.cwl",
                    "scatter": "#main/longreads_quality_kraken2_krona/kraken",
                    "in": [
                        {
                            "source": "#main/longreads_quality_kraken2/sample_report",
                            "id": "#main/longreads_quality_kraken2_krona/kraken"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/longreads_quality_kraken2_krona/kraken2_database"
                        }
                    ],
                    "out": [
                        "#main/longreads_quality_kraken2_krona/krona_html"
                    ],
                    "id": "#main/longreads_quality_kraken2_krona"
                },
                {
                    "label": "Merge fastq files",
                    "when": "$(inputs.longreads.length > 1)",
                    "run": "#concatenate.cwl",
                    "in": [
                        {
                            "source": "#main/identifier",
                            "id": "#main/merge_longreads_fastq/identifier"
                        },
                        {
                            "source": "#main/longreads",
                            "id": "#main/merge_longreads_fastq/infiles"
                        },
                        {
                            "source": "#main/longreads",
                            "id": "#main/merge_longreads_fastq/longreads"
                        },
                        {
                            "valueFrom": "$(inputs.identifier)_$(inputs.readtype)_merged_raw.fastq.gz",
                            "id": "#main/merge_longreads_fastq/outname"
                        },
                        {
                            "source": "#main/readtype",
                            "id": "#main/merge_longreads_fastq/readtype"
                        }
                    ],
                    "out": [
                        "#main/merge_longreads_fastq/output"
                    ],
                    "id": "#main/merge_longreads_fastq"
                },
                {
                    "label": "Prepare references",
                    "doc": "Prepare references to a single fasta file and unique headers",
                    "when": "$(inputs.fasta_input !== null && inputs.fasta_input.length !== 0)",
                    "run": "#workflow_prepare_fasta_db.cwl",
                    "in": [
                        {
                            "source": "#main/filter_references",
                            "id": "#main/prepare_fasta_db/fasta_input"
                        },
                        {
                            "source": "#main/prepare_reference",
                            "id": "#main/prepare_fasta_db/make_headers_unique"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/prepare_fasta_db/output_name"
                        }
                    ],
                    "out": [
                        "#main/prepare_fasta_db/fasta_db"
                    ],
                    "id": "#main/prepare_fasta_db"
                },
                {
                    "label": "Reference mapping",
                    "doc": "Removal of contaminated reads using minimap2 mapping",
                    "when": "$(inputs.filter_references !== null && inputs.filter_references.length !== 0)",
                    "run": "#minimap2_to_fastq.cwl",
                    "in": [
                        {
                            "source": "#main/filter_references",
                            "id": "#main/reference_filter_longreads/filter_references"
                        },
                        {
                            "source": "#main/identifier",
                            "valueFrom": "$(self+\"_\"+inputs.readtype)",
                            "id": "#main/reference_filter_longreads/identifier"
                        },
                        {
                            "source": "#main/keep_reference_mapped_reads",
                            "id": "#main/reference_filter_longreads/output_mapped"
                        },
                        {
                            "default": "map-ont",
                            "id": "#main/reference_filter_longreads/preset"
                        },
                        {
                            "source": "#main/filtlong/output_reads",
                            "id": "#main/reference_filter_longreads/reads"
                        },
                        {
                            "source": "#main/readtype",
                            "id": "#main/reference_filter_longreads/readtype"
                        },
                        {
                            "source": "#main/prepare_fasta_db/fasta_db",
                            "id": "#main/reference_filter_longreads/reference"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/reference_filter_longreads/threads"
                        }
                    ],
                    "out": [
                        "#main/reference_filter_longreads/fastq",
                        "#main/reference_filter_longreads/log"
                    ],
                    "id": "#main/reference_filter_longreads"
                },
                {
                    "label": "Reports to folder",
                    "doc": "Preparation of fastp output files to a specific output folder",
                    "in": [
                        {
                            "source": "#main/step",
                            "valueFrom": "$(self+\"_Longreads_Read_Quality\")\n",
                            "id": "#main/reports_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/fastqc_longreads_before/html_files",
                                "#main/fastqc_longreads_before/zip_files",
                                "#main/filtlong/log",
                                "#main/fastqc_longreads_after/html_files",
                                "#main/fastqc_longreads_after/zip_files",
                                "#main/longreads_quality_kraken2/sample_report",
                                "#main/longreads_quality_kraken2_krona/krona_html",
                                "#main/reference_filter_longreads/log"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/reports_files_to_folder/files"
                        }
                    ],
                    "run": "#files_to_folder.cwl",
                    "out": [
                        "#main/reports_files_to_folder/results"
                    ],
                    "id": "#main/reports_files_to_folder"
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
            "https://schema.org/dateModified": "2023-01-00",
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
                    "class": "StepInputExpressionRequirement"
                }
            ],
            "label": "Prepare (multiple) fasta files to one file.",
            "doc": "Prepare (multiple) fasta files to one file. \nWith option to make unique headers to avoid same fasta headers, which can break some tools.\n",
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "label": "Fasta input",
                    "doc": "Fasta file(s) to prepare",
                    "id": "#workflow_prepare_fasta_db.cwl/fasta_input"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Make headers unique",
                    "doc": "Make fasta headers unique avoiding same fasta headers, which can break some tools.",
                    "default": false,
                    "id": "#workflow_prepare_fasta_db.cwl/make_headers_unique"
                },
                {
                    "type": "string",
                    "doc": "Output name for this dataset used",
                    "label": "identifier used",
                    "id": "#workflow_prepare_fasta_db.cwl/output_name"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "label": "Prepared fasta file",
                    "doc": "Prepared fasta file",
                    "outputSource": [
                        "#workflow_prepare_fasta_db.cwl/fasta_array_to_file/file",
                        "#workflow_prepare_fasta_db.cwl/merge_input/output",
                        "#workflow_prepare_fasta_db.cwl/prepare_fasta_db/fasta_db"
                    ],
                    "pickValue": "first_non_null",
                    "id": "#workflow_prepare_fasta_db.cwl/fasta_db"
                }
            ],
            "steps": [
                {
                    "label": "Array to file",
                    "doc": "Pick first file of filter_reference when make_headers_unique input is false",
                    "when": "$(inputs.make_headers_unique === false && inputs.fasta_input.length === 1)",
                    "run": "#array_to_file.cwl",
                    "in": [
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/fasta_input",
                            "id": "#workflow_prepare_fasta_db.cwl/fasta_array_to_file/fasta_input"
                        },
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/fasta_input",
                            "id": "#workflow_prepare_fasta_db.cwl/fasta_array_to_file/files"
                        },
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/make_headers_unique",
                            "id": "#workflow_prepare_fasta_db.cwl/fasta_array_to_file/make_headers_unique"
                        }
                    ],
                    "out": [
                        "#workflow_prepare_fasta_db.cwl/fasta_array_to_file/file"
                    ],
                    "id": "#workflow_prepare_fasta_db.cwl/fasta_array_to_file"
                },
                {
                    "label": "Merge reference files",
                    "doc": "Only merge input when make unique is false.",
                    "when": "$(inputs.make_headers_unique === false && inputs.fasta_input.length > 1)",
                    "run": "#concatenate.cwl",
                    "in": [
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/fasta_input",
                            "id": "#workflow_prepare_fasta_db.cwl/merge_input/fasta_input"
                        },
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/fasta_input",
                            "id": "#workflow_prepare_fasta_db.cwl/merge_input/infiles"
                        },
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/make_headers_unique",
                            "id": "#workflow_prepare_fasta_db.cwl/merge_input/make_headers_unique"
                        },
                        {
                            "valueFrom": "$(inputs.output_name)_filter-reference_merged.fa",
                            "id": "#workflow_prepare_fasta_db.cwl/merge_input/outname"
                        },
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/output_name",
                            "id": "#workflow_prepare_fasta_db.cwl/merge_input/output_name"
                        }
                    ],
                    "out": [
                        "#workflow_prepare_fasta_db.cwl/merge_input/output"
                    ],
                    "id": "#workflow_prepare_fasta_db.cwl/merge_input"
                },
                {
                    "label": "Prepare references",
                    "doc": "Prepare references to a single fasta file and unique headers",
                    "when": "$(inputs.make_headers_unique)",
                    "run": "#prepare_fasta_db.cwl",
                    "in": [
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/fasta_input",
                            "id": "#workflow_prepare_fasta_db.cwl/prepare_fasta_db/fasta_files"
                        },
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/fasta_input",
                            "id": "#workflow_prepare_fasta_db.cwl/prepare_fasta_db/fasta_input"
                        },
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/make_headers_unique",
                            "id": "#workflow_prepare_fasta_db.cwl/prepare_fasta_db/make_headers_unique"
                        },
                        {
                            "valueFrom": "$(inputs.output_name)_filter-reference_uniq.fa.gz",
                            "id": "#workflow_prepare_fasta_db.cwl/prepare_fasta_db/output_file_name"
                        },
                        {
                            "source": "#workflow_prepare_fasta_db.cwl/output_name",
                            "id": "#workflow_prepare_fasta_db.cwl/prepare_fasta_db/output_name"
                        }
                    ],
                    "out": [
                        "#workflow_prepare_fasta_db.cwl/prepare_fasta_db/fasta_db"
                    ],
                    "id": "#workflow_prepare_fasta_db.cwl/prepare_fasta_db"
                }
            ],
            "id": "#workflow_prepare_fasta_db.cwl",
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
            "https://schema.org/dateCreated": "2023-01-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        }
    ],
    "cwlVersion": "v1.2",
    "$namespaces": {
        "s": "https://schema.org/"
    }
}
