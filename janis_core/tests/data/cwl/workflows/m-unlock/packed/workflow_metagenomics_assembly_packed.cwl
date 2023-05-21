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
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.outname)"
                    },
                    "id": "#concatenate.cwl/output"
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
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
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
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.inputfile.basename).gz"
                    },
                    "id": "#pigz.cwl/outfile"
                }
            ],
            "id": "#pigz.cwl"
        },
        {
            "class": "CommandLineTool",
            "label": "Filter from reads",
            "doc": "Filter reads using BBmaps bbduk tool (paired-end only)\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/bbmap:39.01",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "39.01"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/bbmap"
                            ],
                            "package": "bbmap"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "in=",
                        "separate": false
                    },
                    "id": "#bbduk_filter.cwl/forward_reads"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#bbduk_filter.cwl/identifier"
                },
                {
                    "type": "int",
                    "inputBinding": {
                        "prefix": "k=",
                        "separate": false
                    },
                    "default": 31,
                    "id": "#bbduk_filter.cwl/kmersize"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 8,
                    "id": "#bbduk_filter.cwl/memory"
                },
                {
                    "doc": "Reference contamination fasta file (can be compressed)",
                    "label": "Reference contamination file",
                    "type": [
                        "null",
                        "string"
                    ],
                    "inputBinding": {
                        "prefix": "ref=",
                        "separate": false
                    },
                    "id": "#bbduk_filter.cwl/reference"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "prefix": "in2=",
                        "separate": false
                    },
                    "id": "#bbduk_filter.cwl/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 1,
                    "inputBinding": {
                        "prefix": "threads=",
                        "separate": false
                    },
                    "id": "#bbduk_filter.cwl/threads"
                }
            ],
            "baseCommand": [
                "bbduk.sh"
            ],
            "arguments": [
                {
                    "prefix": "-Xmx",
                    "separate": false,
                    "valueFrom": "$(inputs.memory)M"
                },
                {
                    "prefix": "out=",
                    "separate": false,
                    "valueFrom": "$(inputs.identifier)_1.fq.gz"
                },
                {
                    "prefix": "out2=",
                    "separate": false,
                    "valueFrom": "$(inputs.identifier)_2.fq.gz"
                },
                {
                    "prefix": "stats=",
                    "separate": false,
                    "valueFrom": "$(inputs.identifier)_bbduk-stats.txt"
                }
            ],
            "stderr": "$(inputs.identifier)_bbduk-summary.txt",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_1.fq.gz"
                    },
                    "id": "#bbduk_filter.cwl/out_forward_reads"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_2.fq.gz"
                    },
                    "id": "#bbduk_filter.cwl/out_reverse_reads"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_bbduk-stats.txt"
                    },
                    "id": "#bbduk_filter.cwl/stats_file"
                },
                {
                    "type": "File",
                    "id": "#bbduk_filter.cwl/summary",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_bbduk-summary.txt"
                    }
                }
            ],
            "id": "#bbduk_filter.cwl",
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
            "https://schema.org/dateModified": "2023-02-07",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/bbmap:38.98",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "38.98"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/bbmap"
                            ],
                            "package": "bbmap"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "label": "BBMap",
            "doc": "Read filtering using BBMap against a (contamination) reference genome\n",
            "baseCommand": [
                "bbmap.sh"
            ],
            "arguments": [
                "-Xmx$(inputs.memory)M",
                "printunmappedcount",
                "overwrite=true",
                "statsfile=$(inputs.identifier)_BBMap_stats.txt",
                "covstats=$(inputs.identifier)_BBMap_covstats.txt",
                "out=$(inputs.identifier)_BBMap.sam"
            ],
            "inputs": [
                {
                    "type": "boolean",
                    "label": "fast mode",
                    "doc": "Sets other BBMap paramters to run faster, at reduced sensitivity. Bad for RNA-seq.",
                    "inputBinding": {
                        "prefix": "fast=t"
                    },
                    "default": true,
                    "id": "#bbmap.cwl/fast"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "in=",
                        "separate": false
                    },
                    "id": "#bbmap.cwl/forward_reads"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#bbmap.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "maximum memory usage in megabytes",
                    "label": "memory usage (mb)",
                    "default": 8000,
                    "id": "#bbmap.cwl/memory"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "ref=",
                        "separate": false
                    },
                    "id": "#bbmap.cwl/reference"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "in2=",
                        "separate": false
                    },
                    "id": "#bbmap.cwl/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "number of threads to use for computational processes",
                    "label": "number of threads",
                    "inputBinding": {
                        "prefix": "threads=",
                        "separate": false
                    },
                    "default": 2,
                    "id": "#bbmap.cwl/threads"
                }
            ],
            "stderr": "$(inputs.identifier)_BBMap_log.txt",
            "outputs": [
                {
                    "label": "Coverage per contig",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_BBMap_covstats.txt"
                    },
                    "id": "#bbmap.cwl/covstats"
                },
                {
                    "label": "BBMap log output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_BBMap_log.txt"
                    },
                    "id": "#bbmap.cwl/log"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_BBMap.sam"
                    },
                    "id": "#bbmap.cwl/sam"
                },
                {
                    "label": "Mapping statistics",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_BBMap_stats.txt"
                    },
                    "id": "#bbmap.cwl/stats"
                }
            ],
            "id": "#bbmap.cwl",
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
            "https://schema.org/dateModified": "2022-04-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/bbmap:38.98",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "38.98"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/bbmap"
                            ],
                            "package": "bbmap"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "label": "BBMap",
            "doc": "Read filtering using BBMap against a (contamination) reference genome\n",
            "inputs": [
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 1,
                        "prefix": "in=",
                        "separate": false
                    },
                    "id": "#bbmap_filter-reads.cwl/forward_reads"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#bbmap_filter-reads.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "maximum memory usage in megabytes",
                    "label": "memory usage (mb)",
                    "default": 8000,
                    "id": "#bbmap_filter-reads.cwl/memory"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "default": false,
                    "id": "#bbmap_filter-reads.cwl/output_mapped"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "position": 3,
                        "prefix": "ref=",
                        "separate": false
                    },
                    "id": "#bbmap_filter-reads.cwl/reference"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "inputBinding": {
                        "position": 2,
                        "prefix": "in2=",
                        "separate": false
                    },
                    "id": "#bbmap_filter-reads.cwl/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "number of threads to use for computational processes",
                    "label": "number of threads",
                    "inputBinding": {
                        "prefix": "threads=",
                        "separate": false
                    },
                    "default": 2,
                    "id": "#bbmap_filter-reads.cwl/threads"
                }
            ],
            "stderr": "$(inputs.identifier)_BBMap_log.txt",
            "outputs": [
                {
                    "label": "Coverage per contig",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_BBMap_covstats.txt"
                    },
                    "id": "#bbmap_filter-reads.cwl/covstats"
                },
                {
                    "label": "BBMap log output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_BBMap_log.txt"
                    },
                    "id": "#bbmap_filter-reads.cwl/log"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_filtered_1.fq.gz"
                    },
                    "id": "#bbmap_filter-reads.cwl/out_forward_reads"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_filtered_2.fq.gz"
                    },
                    "id": "#bbmap_filter-reads.cwl/out_reverse_reads"
                },
                {
                    "label": "Mapping statistics",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_BBMap_stats.txt"
                    },
                    "id": "#bbmap_filter-reads.cwl/stats"
                }
            ],
            "baseCommand": [
                "bbmap.sh"
            ],
            "arguments": [
                "-Xmx$(inputs.memory)M",
                "printunmappedcount",
                "overwrite=true",
                "bloom=t",
                "statsfile=$(inputs.identifier)_BBMap_stats.txt",
                "covstats=$(inputs.identifier)_BBMap_covstats.txt",
                "${\n  if (inputs.output_mapped){\n    return 'outm1='+inputs.identifier+'_filtered_1.fq.gz \\\n            outm2='+inputs.identifier+'_filtered_2.fq.gz';\n  } else {\n    return 'outu1='+inputs.identifier+'_filtered_1.fq.gz \\\n            outu2='+inputs.identifier+'_filtered_2.fq.gz';\n  }\n}\n"
            ],
            "id": "#bbmap_filter-reads.cwl",
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
            "https://schema.org/dateModified": "2022-04-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
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
            "baseCommand": [
                "bash",
                "-x",
                "script.sh"
            ],
            "label": "CarveMe GEMstats",
            "doc": "Small summary of a list of CarveMe genome-scale metabolic models in sbml-fbc2 format\nContains; number of mets,reactions and genes\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "gemstats_output",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nidentifier=$1\nshift;\necho \"Model Mets Reactions Genes\" > $identifier\\_CarveMe_GEMstats.tsv\nfor file in \"$@\"\ndo\n  bash /unlock/infrastructure/scripts/GEMstats.sh $file\ndone >> $identifier\\_CarveMe_GEMstats.tsv"
                        }
                    ]
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "label": "CarveMe GEMs",
                    "doc": "List of CarveMe metabolic models in sbml-fbc2 format.",
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#GEMstats.cwl/carveme_gems"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#GEMstats.cwl/identifier"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_CarveMe_GEMstats.tsv"
                    },
                    "id": "#GEMstats.cwl/carveme_GEMstats"
                }
            ],
            "id": "#GEMstats.cwl",
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
            "https://schema.org/dateCreated": "2022-01-01",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "baseCommand": [
                "carve"
            ],
            "label": "CarveMe",
            "doc": "CarveMe is a python-based tool for genome-scale metabolic model reconstruction.\n(Workflow will quit as successful even though no model can be created. Check messages.)\n    \n",
            "requirements": [
                {
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "carve_output",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nsource /root/miniconda/bin/activate\nconda init bash\nconda activate /unlock/infrastructure/conda/carveme/cplex/carveme_1.5.1\ncarve $@"
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
                    "packages": [
                        {
                            "version": [
                                "1.5.1"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/carveme"
                            ],
                            "package": "carveme"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "outputs": [
                {
                    "label": "CarveMe GEM",
                    "doc": "CarveMe metabolic model Output SBML in sbml-fbc2 format",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.protein_file.nameroot).GEM.xml"
                    },
                    "id": "#carveme.cwl/carveme_gem"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Gap fill",
                    "doc": "Gap fill model for given media",
                    "inputBinding": {
                        "prefix": "--gapfill"
                    },
                    "id": "#carveme.cwl/gapfill"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Initial media",
                    "doc": "Initialize model with given medium",
                    "inputBinding": {
                        "prefix": "--init"
                    },
                    "id": "#carveme.cwl/init"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Media database",
                    "doc": "Media database file",
                    "inputBinding": {
                        "prefix": "--mediadb"
                    },
                    "id": "#carveme.cwl/mediadb"
                },
                {
                    "type": "File",
                    "label": "Input fasta file",
                    "doc": "Proteins sequence file in FASTA format.",
                    "inputBinding": {
                        "position": 0
                    },
                    "id": "#carveme.cwl/protein_file"
                }
            ],
            "arguments": [
                "--fbc2",
                {
                    "prefix": "--output",
                    "valueFrom": "$(inputs.protein_file.nameroot).GEM.xml"
                }
            ],
            "id": "#carveme.cwl",
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
            "doc": "Modified from https://github.com/ambarishK/bio-cwl-tools/blob/release/fastp/fastp.cwl\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/fastp:0.23.2",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "0.23.2"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/fastp"
                            ],
                            "package": "fastp"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "default": true,
                    "inputBinding": {
                        "prefix": "--correction"
                    },
                    "id": "#fastp.cwl/base_correction"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "default": false,
                    "inputBinding": {
                        "prefix": "--dedup"
                    },
                    "id": "#fastp.cwl/deduplicate"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "default": true,
                    "inputBinding": {
                        "prefix": "--disable_trim_poly_g"
                    },
                    "id": "#fastp.cwl/disable_trim_poly_g"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "inputBinding": {
                        "prefix": "--trim_poly_g"
                    },
                    "id": "#fastp.cwl/force_polyg_tail_trimming"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--in1"
                    },
                    "id": "#fastp.cwl/forward_reads"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#fastp.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "default": false,
                    "inputBinding": {
                        "prefix": "--merge"
                    },
                    "id": "#fastp.cwl/merge_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 50,
                    "inputBinding": {
                        "prefix": "--length_required"
                    },
                    "id": "#fastp.cwl/min_length_required"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 20,
                    "inputBinding": {
                        "prefix": "--qualified_quality_phred"
                    },
                    "id": "#fastp.cwl/qualified_phred_quality"
                },
                {
                    "type": "File",
                    "inputBinding": {
                        "prefix": "--in2"
                    },
                    "id": "#fastp.cwl/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--thread"
                    },
                    "id": "#fastp.cwl/threads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 20,
                    "inputBinding": {
                        "prefix": "--unqualified_percent_limit"
                    },
                    "id": "#fastp.cwl/unqualified_phred_quality"
                }
            ],
            "arguments": [
                {
                    "prefix": "--out1",
                    "valueFrom": "$(inputs.identifier)_fastp_1.fq.gz"
                },
                "${\n  if (inputs.reverse_reads){\n    return '--out2';\n  } else {\n    return '';\n  }\n}\n",
                "${\n  if (inputs.reverse_reads){\n    return inputs.identifier + \"_fastp_2.fq.gz\";\n  } else {\n    return '';\n  }\n}\n",
                "${\n  if (inputs.reverse_reads_path){\n    return '--out2';\n  } else {\n    return '';\n  }\n}\n",
                "${\n  if (inputs.reverse_reads_path){\n    return inputs.identifier + \"_fastp_2.fq.gz\";\n  } else {\n    return '';\n  }\n}\n",
                "${\n  if (inputs.merge_reads){\n    return '--merged_out';\n  } else {\n    return '';\n  }\n}\n",
                "${\n  if (inputs.merge_reads){\n    return inputs.identifier + \"merged_fastp.fq.gz\";\n  } else {\n    return '';\n  }\n}\n",
                {
                    "prefix": "-h",
                    "valueFrom": "$(inputs.identifier)_fastp.html"
                },
                {
                    "prefix": "-j",
                    "valueFrom": "$(inputs.identifier)_fastp.json"
                }
            ],
            "baseCommand": [
                "fastp"
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_fastp.html"
                    },
                    "id": "#fastp.cwl/html_report"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_fastp.json"
                    },
                    "id": "#fastp.cwl/json_report"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_merged_fastp.fq.gz"
                    },
                    "id": "#fastp.cwl/merged_reads"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_fastp_1.fq.gz"
                    },
                    "id": "#fastp.cwl/out_forward_reads"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_fastp_2.fq.gz"
                    },
                    "id": "#fastp.cwl/out_reverse_reads"
                }
            ],
            "id": "#fastp.cwl",
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
            "https://schema.org/dateModified": "2022-02-22",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
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
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/flye:2.9.1",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "2.9.1"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/flye"
                            ],
                            "package": "flye"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "flye"
            ],
            "label": "Flye",
            "doc": "Flye De novo assembler for single molecule sequencing reads, with a focus in Oxford Nanopore Technologies reads",
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Debug mode",
                    "doc": "Set to true to display debug output while running",
                    "default": false,
                    "inputBinding": {
                        "prefix": "--debug"
                    },
                    "id": "#flye.cwl/debug_mode"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Genome size",
                    "doc": "Estimated genome size (for example, 5m or 2.6g)",
                    "inputBinding": {
                        "prefix": "--genome-size"
                    },
                    "id": "#flye.cwl/genome_size"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Metagenome",
                    "doc": "Set to true if assembling a metagenome",
                    "default": false,
                    "inputBinding": {
                        "prefix": "--meta"
                    },
                    "id": "#flye.cwl/metagenome"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "ONT corrected",
                    "doc": "ONT reads in FASTQ format that were corrected with other methods (<3% error)",
                    "inputBinding": {
                        "prefix": "--nano-corr"
                    },
                    "id": "#flye.cwl/nano_corrected"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "ONT high quality",
                    "doc": "ONT high-quality reads in FASTQ format, Guppy5 SUP or Q20 (<5% error)",
                    "inputBinding": {
                        "prefix": "--nano-hq"
                    },
                    "id": "#flye.cwl/nano_high_quality"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "ONT reads raw",
                    "doc": "ONT regular reads in FASTQ format, pre-Guppy5 (<20% error)",
                    "inputBinding": {
                        "prefix": "--nano-raw"
                    },
                    "id": "#flye.cwl/nano_raw"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output folder name",
                    "doc": "Output folder name",
                    "inputBinding": {
                        "prefix": "--out-dir"
                    },
                    "default": "flye_output",
                    "id": "#flye.cwl/output_folder_name"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "PacBio reads corrected",
                    "doc": "PacBio  reads in FASTQ format, that were corrected with other methods (<3% error)",
                    "inputBinding": {
                        "prefix": "--pacbio-corr"
                    },
                    "id": "#flye.cwl/pacbio_corrected"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "PacBio HiFi reads",
                    "doc": "PacBio HiFi  reads in FASTQ format, (<1% error)",
                    "inputBinding": {
                        "prefix": "--pacbio-hifi"
                    },
                    "id": "#flye.cwl/pacbio_hifi"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "PacBio reads raw",
                    "doc": "PacBio regular CLR  reads in FASTQ format, (<20% error)",
                    "inputBinding": {
                        "prefix": "--pacbio-raw"
                    },
                    "id": "#flye.cwl/pacbio_raw"
                },
                {
                    "label": "Flye will carry out polishing multiple times as determined here",
                    "type": [
                        "null",
                        "int"
                    ],
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--iterations"
                    },
                    "id": "#flye.cwl/polishing_iterations"
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
                    "id": "#flye.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/00-assembly"
                    },
                    "id": "#flye.cwl/00_assembly"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/10-consensus"
                    },
                    "id": "#flye.cwl/10_consensus"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/20-repeat"
                    },
                    "id": "#flye.cwl/20_repeat"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/30-contigger"
                    },
                    "id": "#flye.cwl/30_contigger"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/40-polishing"
                    },
                    "id": "#flye.cwl/40_polishing"
                },
                {
                    "label": "Polished assembly created by flye, main output for after polishing with next tool",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/assembly.fasta"
                    },
                    "id": "#flye.cwl/assembly"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/assembly_info.txt"
                    },
                    "id": "#flye.cwl/assembly_info"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/flye.log"
                    },
                    "id": "#flye.cwl/flye_log"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.output_folder_name)/params.json"
                    },
                    "id": "#flye.cwl/params"
                }
            ],
            "id": "#flye.cwl",
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
            "https://schema.org/dateCreated": "2021-11-29",
            "https://schema.org/dateModified": "2023-01-10",
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
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/medaka:1.7.0",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "1.7.0"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/medaka"
                            ],
                            "package": "medaka"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "baseCommand": [
                "medaka.py"
            ],
            "label": "Polishing of assembly created from ONT nanopore long reads",
            "doc": "Uses Medaka to polish an assembly constructed from ONT nanopore reads\n  \n",
            "requirements": [
                {
                    "listing": [
                        "$(inputs.draft_assembly)"
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "inputs": [
                {
                    "label": "Basecalling model that was used by guppy.",
                    "doc": "Please consult https://github.com/nanoporetech/medaka for detailed information.\nChoice of medaka model depends on how basecalling was performed. Run \"medaka tools list\\_models\".\n{pore}_{device}_{caller variant}_{caller version}\nAvailable models: r941_trans, r941_flip213, r941_flip235, r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high\nFor basecalling with guppy version >= v3.0.3, select model based on pore name and whether high or fast basecalling was used.\nFor flip flop basecalling with v3.03 > guppy version => v2.3.5 select r941_flip235.\nFor flip flop basecalling with v2.3.5 > guppy version >= 2.1.3 select r941_flip213.\nFor transducer basecaling using Albacore or non-flip-flop guppy basecalling, select r941_trans.\n\nFor test set (https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html?highlight=flowcell),\nuse \"r941_min_hac_g507\" according to the list of available models.\n",
                    "type": "string",
                    "default": "r941_min_hac_g507",
                    "inputBinding": {
                        "prefix": "--model"
                    },
                    "id": "#medaka_py.cwl/basecall_model"
                },
                {
                    "label": "Assembly that medaka will try to polish.",
                    "type": "File",
                    "inputBinding": {
                        "prefix": "-r"
                    },
                    "id": "#medaka_py.cwl/draft_assembly"
                },
                {
                    "label": "Basecalled ONT nanopore reads in fastq format.",
                    "type": "File",
                    "inputBinding": {
                        "prefix": "-i"
                    },
                    "id": "#medaka_py.cwl/reads"
                },
                {
                    "label": "Number of CPU threads used by tool.",
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "prefix": "-t"
                    },
                    "id": "#medaka_py.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "polished_*bed"
                    },
                    "id": "#medaka_py.cwl/gaps_in_draft_coords"
                },
                {
                    "label": "draft assembly polished by medaka.",
                    "type": "File",
                    "outputBinding": {
                        "glob": "polished_*fasta"
                    },
                    "id": "#medaka_py.cwl/polished_assembly"
                }
            ],
            "id": "#medaka_py.cwl",
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
            "https://schema.org/dateCreated": "2021-11-29",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "MEMOTE",
            "doc": "MEMOTE, short for metabolic model testing",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "MEMOTE",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nsource /root/miniconda/bin/activate\nconda init bash\nconda activate /unlock/infrastructure/conda/memote/cplex/memote_0.13.0\nmemote $@"
                        }
                    ]
                }
            ],
            "baseCommand": [
                "bash",
                "script.sh"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Metabolic model (xml format)",
                    "label": "Metabolic model",
                    "inputBinding": {
                        "position": 100
                    },
                    "id": "#memote.cwl/GEM"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Report",
                    "doc": "Take a snapshot of a model's state and generate a report.",
                    "inputBinding": {
                        "position": 0
                    },
                    "id": "#memote.cwl/report_snapshot"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Run",
                    "doc": "Run the test suite on a single model and collect results.",
                    "inputBinding": {
                        "position": 0
                    },
                    "id": "#memote.cwl/run"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; find incorrect thermodynamic reversibility",
                    "label": "Skip test; find incorrect thermodynamic reversibility",
                    "inputBinding": {
                        "prefix": "--skip skip_test_find_incorrect_thermodynamic_reversibility",
                        "position": 15
                    },
                    "id": "#memote.cwl/skip_test_find_incorrect_thermodynamic_reversibility"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; Find metabolites consumed with closed bounds",
                    "label": "Skip test; Find metabolites consumed with closed bounds",
                    "inputBinding": {
                        "prefix": "--skip test_find_metabolites_consumed_with_closed_bounds",
                        "position": 12
                    },
                    "id": "#memote.cwl/skip_test_find_metabolites_consumed_with_closed_bounds"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; Find metabolites not consumedwith open_bounds",
                    "label": "Skip test; Find metabolites not consumedwith open_bounds",
                    "inputBinding": {
                        "prefix": "--skip test_find_metabolites_not_consumed_with_open_bounds",
                        "position": 14
                    },
                    "id": "#memote.cwl/skip_test_find_metabolites_not_consumed_with_open_bounds"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; Find metabolites not produced with open bounds",
                    "label": "Skip test; Find metabolites not produced with open bounds",
                    "inputBinding": {
                        "prefix": "--skip test_find_metabolites_not_produced_with_open_bounds",
                        "position": 13
                    },
                    "id": "#memote.cwl/skip_test_find_metabolites_not_produced_with_open_bounds"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip test; find metabolites produced with closed bounds",
                    "label": "Skip test; find metabolites produced with closed bounds",
                    "inputBinding": {
                        "prefix": "--skip test_find_metabolites_produced_with_closed_bounds",
                        "position": 11
                    },
                    "id": "#memote.cwl/skip_test_find_metabolites_produced_with_closed_bounds"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "solver",
                    "doc": "Set the solver to be used. [cplex|glpk|gurobi|glpk_exact]. default; glpk",
                    "inputBinding": {
                        "prefix": "--solver",
                        "position": 3
                    },
                    "id": "#memote.cwl/solver"
                }
            ],
            "arguments": [
                "${\n  if (inputs.run){\n    return \"run --filename \" + inputs.GEM.basename + \"_MEMOTE.json.gz\";\n  } else {\n    return '';\n  }\n}\n",
                "${\n  if (inputs.report_snapshot){\n    return \"report snapshot --filename \" + inputs.GEM.basename + \"_MEMOTE.html\";\n  } else {\n    return '';\n  }\n}\n"
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.GEM.basename)_MEMOTE.html"
                    },
                    "id": "#memote.cwl/report_html"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "$(inputs.GEM.basename)_MEMOTE.json.gz"
                    },
                    "id": "#memote.cwl/run_json"
                }
            ],
            "id": "#memote.cwl",
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
            "class": "CommandLineTool",
            "label": "Pilon",
            "doc": "\"https://github.com/broadinstitute/pilon\n Pilon is a software tool which can be used to:\n   Automatically improve draft assemblies\n   Find variation among strains, including large event detection\"\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/pilon:1.24",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "1.24"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/pilon"
                            ],
                            "package": "pilon"
                        }
                    ],
                    "class": "SoftwareRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "label": "Assembly",
                    "doc": "Draft assembly fasta file",
                    "inputBinding": {
                        "prefix": "--genome"
                    },
                    "id": "#pilon.cwl/assembly"
                },
                {
                    "type": "File",
                    "label": "Bam file",
                    "doc": "Indexed sorted bam file with mapped reads to draft assembly",
                    "inputBinding": {
                        "prefix": "--frags"
                    },
                    "id": "#pilon.cwl/bam_file"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Fix List",
                    "doc": "A comma-separated list of categories of issues to try to fix:\n  \"snps\": try to fix individual base errors;\n  \"indels\": try to fix small indels;\n  \"gaps\": try to fill gaps;\n  \"local\": try to detect and fix local misassemblies;\n  \"all\": all of the above (default);\n  \"bases\": shorthand for \"snps\" and \"indels\" (for back compatibility);\n",
                    "inputBinding": {
                        "prefix": "--fix"
                    },
                    "id": "#pilon.cwl/fixlist"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#pilon.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Memory usage in megabytes",
                    "label": "memory usage (MB)",
                    "default": 8000,
                    "id": "#pilon.cwl/memory"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "number of threads to use for computational processes",
                    "label": "number of threads",
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "default": 2,
                    "id": "#pilon.cwl/threads"
                },
                {
                    "label": "vcf output",
                    "type": "boolean",
                    "default": true,
                    "inputBinding": {
                        "prefix": "--vcf"
                    },
                    "id": "#pilon.cwl/vcf"
                }
            ],
            "baseCommand": [
                "java"
            ],
            "arguments": [
                "-jar",
                "-Xmx$(inputs.memory)M",
                "/venv/share/pilon-1.24-0/pilon.jar",
                {
                    "valueFrom": "$(inputs.identifier)_pilon_polished",
                    "prefix": "--output"
                }
            ],
            "stdout": "$(inputs.identifier)_pilon.log",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_pilon.log"
                    },
                    "id": "#pilon.cwl/pilon_log"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_pilon_polished.fasta"
                    },
                    "id": "#pilon.cwl/pilon_polished_assembly"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_pilon_polished.vcf"
                    },
                    "id": "#pilon.cwl/pilon_vcf"
                }
            ],
            "id": "#pilon.cwl",
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
            "https://schema.org/dateCreated": "2022-02-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
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
            "baseCommand": [
                "metaquast.py"
            ],
            "label": "metaQUAST: Quality Assessment Tool for Metagenome Assemblies",
            "doc": "Runs the Quality Assessment Tool for Metagenome Assemblies application\n\nNecessary to install the pre-release to prevent issues:\nhttps://github.com/ablab/quast/releases/tag/quast_5.1.0rc1\n\nThe working installation followed the method in http://quast.sourceforge.net/docs/manual.html:\n$ wget https://github.com/ablab/quast/releases/download/quast_5.1.0rc1/quast-5.1.0rc1.tar.gz\n$ tar -xzf quast-5.1.0rc1.tar.gz\n$ cd quast-5.1.0rc1/\n$ ./setup.py install_full\n",
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/quast:5.2.0",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
                        {
                            "version": [
                                "5.2.0"
                            ],
                            "specs": [
                                "https://anaconda.org/bioconda/quast"
                            ],
                            "package": "quast"
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
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "metaQUAST_results",
                            "writable": true
                        }
                    ]
                }
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.assembly)"
                },
                {
                    "valueFrom": "metaQUAST_results",
                    "prefix": "--output-dir"
                },
                "${\n  if (inputs.blastdb){\n    return [\"--blast-db\", inputs.blastdb.path];\n  } else {\n    return [\"--blast-db\", \"/venv/silva_138.1/\"];\n  }\n}\n"
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "The input assembly in fasta format",
                    "id": "#metaquast.cwl/assembly"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "Reference BLASTdb",
                    "doc": "Reference BLASTdb. The path should point to a directory containing .nsq file or to .nsq file itself. (default is silva)",
                    "inputBinding": {
                        "prefix": "--blast-db"
                    },
                    "id": "#metaquast.cwl/blastdb"
                },
                {
                    "type": "boolean",
                    "label": "Run space efficient",
                    "doc": "Uses less diskspace. Removes aux files with some details for advanced analysis. (default true)",
                    "inputBinding": {
                        "prefix": "--space-efficient"
                    },
                    "default": true,
                    "id": "#metaquast.cwl/space_efficient"
                },
                {
                    "type": "int",
                    "default": 1,
                    "inputBinding": {
                        "prefix": "--threads"
                    },
                    "id": "#metaquast.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/basic_stats"
                    },
                    "id": "#metaquast.cwl/basicStats"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/combined_reference"
                    },
                    "id": "#metaquast.cwl/meta_combined_ref"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/quast_downloaded_references"
                    },
                    "id": "#metaquast.cwl/meta_downloaded_ref"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/icarus.html"
                    },
                    "id": "#metaquast.cwl/meta_icarus"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/icarus_viewers"
                    },
                    "id": "#metaquast.cwl/meta_icarusDir"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/summary"
                    },
                    "id": "#metaquast.cwl/meta_summary"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/krona_charts"
                    },
                    "id": "#metaquast.cwl/metaquast_krona"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/metaquast.log"
                    },
                    "id": "#metaquast.cwl/metaquast_log"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "metaQUAST_results"
                    },
                    "id": "#metaquast.cwl/metaquast_outdir"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/report.html"
                    },
                    "id": "#metaquast.cwl/metaquast_report"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/not_aligned"
                    },
                    "id": "#metaquast.cwl/not_aligned"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/quast.log"
                    },
                    "id": "#metaquast.cwl/quastLog"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/report.*"
                    },
                    "id": "#metaquast.cwl/quastReport"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/icarus_viewers"
                    },
                    "id": "#metaquast.cwl/quast_icarusDir"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/icarus.html"
                    },
                    "id": "#metaquast.cwl/quast_icarusHtml"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/runs_per_reference"
                    },
                    "id": "#metaquast.cwl/runs_per_reference"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "outputBinding": {
                        "glob": "metaQUAST_results/transposed_report.*"
                    },
                    "id": "#metaquast.cwl/transposedReport"
                }
            ],
            "id": "#metaquast.cwl",
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
            "https://schema.org/dateCreated": "2021-02-09",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "ExpressionTool",
            "requirements": [
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        "$(inputs.bam_file)",
                        "$(inputs.bai)"
                    ]
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "id": "#expression_bam_index.cwl/bam_file"
                },
                {
                    "type": "File",
                    "id": "#expression_bam_index.cwl/bam_index"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "id": "#expression_bam_index.cwl/hybrid_bamindex"
                }
            ],
            "expression": "${ var ret = inputs.bam_file; ret[\"secondaryFiles\"] = [\n    inputs.bam_index,\n]; return { \"hybrid_bamindex\": ret } ; }\n",
            "doc": "Creates an \"hybrid\" output with the index file as secondary file.\n\nAdapted from:\nhttps://github.com/vetscience/Assemblosis/blob/master/Run/expressiontoolbam.cwl\n\nLICENSE:\nBSD 3-Clause License\n\nCopyright (c) 2017, The University lf Melbourne\nAll rights reserved.\n\nRedistribution and use in source and binary forms, with or without\nmodification, are permitted provided that the following conditions are met:\n\n* Redistributions of source code must retain the above copyright notice, this\n  list of conditions and the following disclaimer.\n\n* Redistributions in binary form must reproduce the above copyright notice,\n  this list of conditions and the following disclaimer in the documentation\n  and/or other materials provided with the distribution.\n\n* Neither the name of the copyright holder nor the names of its\n  contributors may be used to endorse or promote products derived from\n  this software without specific prior written permission.\n\nTHIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"\nAND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\nIMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\nDISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE\nFOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL\nDAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR\nSERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER\nCAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,\nOR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\nOF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.",
            "id": "#expression_bam_index.cwl",
            "http://schema.org/citation": "https://m-unlock.nl",
            "http://schema.org/codeRepository": "https://gitlab.com/m-unlock/cwl",
            "http://schema.org/dateCreated": "2020-00-00",
            "http://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "http://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "sam to sorted bam",
            "doc": "samtools view -@ $2 -hu $1 | samtools sort -@ $2 -o $3.bam\n",
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
                            "entry": "#!/bin/bash\nsamtools view -@ $2 -hu $3 | samtools sort -@ $2 -o $1.sorted.bam"
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
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/samtools:1.15.1",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
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
            "inputs": [
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#sam_to_sorted-bam.cwl/identifier"
                },
                {
                    "type": "File",
                    "doc": "unsorted sam file",
                    "label": "unsorted sam file",
                    "inputBinding": {
                        "position": 3
                    },
                    "id": "#sam_to_sorted-bam.cwl/sam"
                },
                {
                    "type": "int",
                    "doc": "number of cpu threads to use",
                    "label": "cpu threads",
                    "default": 1,
                    "inputBinding": {
                        "position": 2
                    },
                    "id": "#sam_to_sorted-bam.cwl/threads"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier).sorted.bam"
                    },
                    "id": "#sam_to_sorted-bam.cwl/sortedbam"
                }
            ],
            "id": "#sam_to_sorted-bam.cwl",
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
            "https://schema.org/dateModified": "2022-02-22",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "samtools idxstats",
            "doc": "samtools idxstats - reports alignment summary statistics",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/samtools:1.15.1",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
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
            "inputs": [
                {
                    "type": "File",
                    "label": "Bam file",
                    "doc": "(sorted) Bam file",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#samtools_idxstats.cwl/bam_file"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#samtools_idxstats.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Number of threads to use",
                    "default": 2,
                    "inputBinding": {
                        "position": 0,
                        "prefix": "--threads"
                    },
                    "id": "#samtools_idxstats.cwl/threads"
                }
            ],
            "baseCommand": [
                "samtools",
                "idxstats"
            ],
            "stdout": "$(inputs.identifier)_contigReadCounts.tsv",
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_contigReadCounts.tsv"
                    },
                    "id": "#samtools_idxstats.cwl/contigReadCounts"
                }
            ],
            "id": "#samtools_idxstats.cwl",
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
            "https://schema.org/dateModified": "2022-02-22",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        },
        {
            "class": "CommandLineTool",
            "label": "samtools index",
            "doc": "samtools index - creates an index file for bam file",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/samtools:1.15.1",
                    "class": "DockerRequirement"
                },
                {
                    "packages": [
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
            "inputs": [
                {
                    "type": "File",
                    "label": "Bam file",
                    "doc": "(sorted) Bam file",
                    "inputBinding": {
                        "position": 1
                    },
                    "id": "#samtools_index.cwl/bam_file"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "Number of threads to use",
                    "default": 2,
                    "inputBinding": {
                        "position": 0,
                        "prefix": "-@"
                    },
                    "id": "#samtools_index.cwl/threads"
                }
            ],
            "baseCommand": [
                "samtools",
                "index"
            ],
            "arguments": [
                {
                    "valueFrom": "$(inputs.bam_file.basename).bai",
                    "position": 2
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.bam_file.basename).bai"
                    },
                    "id": "#samtools_index.cwl/bam_index"
                }
            ],
            "id": "#samtools_index.cwl",
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
            "https://schema.org/dateCreated": "2022-02-00",
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
            "baseCommand": [
                "bash",
                "script.sh"
            ],
            "label": "SMETANA",
            "doc": "Species METabolic interaction ANAlysis (SMETANA) is a python-based command line tool to analyse microbial communities.\nIt takes as input a microbial communtity (from a collection of genome-scale metabolic models in SBML format) and \ncomputes of several metrics that describe the potential for cross-feeding interactions between community members.\n    \n",
            "requirements": [
                {
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "smetana_output",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nsource /root/miniconda/bin/activate\nconda init bash\nconda activate /unlock/infrastructure/conda/smetana/cplex/smetana_1.2.0\nsmetana $@"
                        }
                    ],
                    "class": "InitialWorkDirRequirement"
                },
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "outputs": [
                {
                    "label": "SMETANA output",
                    "doc": "SMETANA detailed tsv output",
                    "type": "File",
                    "outputBinding": {
                        "glob": "$(inputs.identifier)_SMETANA*"
                    },
                    "id": "#smetana.cwl/detailed_output_tsv"
                }
            ],
            "inputs": [
                {
                    "label": "Metabolic model",
                    "doc": "Multiple Metabolic models (xml format)",
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "inputBinding": {
                        "position": 10
                    },
                    "id": "#smetana.cwl/GEM"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Runs MIP/MRO and is much faster, recommended when analysing multiple communities.",
                    "label": "Detailed mode",
                    "inputBinding": {
                        "prefix": "--detailed",
                        "position": 11
                    },
                    "id": "#smetana.cwl/detailed"
                },
                {
                    "type": "string",
                    "label": "Flavor",
                    "doc": "Expected SBML flavor of the input files (cobra or fbc2)",
                    "inputBinding": {
                        "prefix": "--flavor",
                        "position": 3
                    },
                    "default": "fbc2",
                    "id": "#smetana.cwl/flavor"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Calculates all inter-species interactions (much slower), check Algorithms for details.",
                    "label": "Global mode",
                    "inputBinding": {
                        "prefix": "--global",
                        "position": 11
                    },
                    "default": false,
                    "id": "#smetana.cwl/global"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "Identifier used",
                    "id": "#smetana.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Media database",
                    "doc": "Media database file",
                    "inputBinding": {
                        "prefix": "-m",
                        "position": 1
                    },
                    "id": "#smetana.cwl/media"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Media database",
                    "doc": "Media database file",
                    "inputBinding": {
                        "prefix": "--mediadb",
                        "position": 2
                    },
                    "id": "#smetana.cwl/mediadb"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "solver",
                    "doc": "Set the solver to be used. [cplex|glpk|gurobi|glpk_exact]. default; glpk",
                    "inputBinding": {
                        "prefix": "--solver",
                        "position": 4
                    },
                    "id": "#smetana.cwl/solver"
                }
            ],
            "arguments": [
                "--verbose",
                {
                    "prefix": "--output",
                    "valueFrom": "$(inputs.identifier)_SMETANA"
                }
            ],
            "id": "#smetana.cwl",
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
            "label": "Illumina read quality control, trimming and contamination filter.",
            "doc": "**Workflow for Illumina paired read quality control, trimming and filtering.**<br />\nMultiple paired datasets will be merged into single paired dataset.<br />\nSummary:\n- FastQC on raw data files<br />\n- fastp for read quality trimming<br />\n- BBduk for phiX and (optional) rRNA filtering<br />\n- Kraken2 for taxonomic classification of reads (optional)<br />\n- BBmap for (contamination) filtering using given references (optional)<br />\n- FastQC on filtered (merged) data<br />\n\nOther UNLOCK workflows on WorkflowHub: https://workflowhub.eu/projects/16/workflows?view=default<br><br>\n\n**All tool CWL files and other workflows can be found here:**<br>\n  Tools: https://gitlab.com/m-unlock/cwl<br>\n  Workflows: https://gitlab.com/m-unlock/cwl/workflows<br>\n\n**How to setup and use an UNLOCK workflow:**<br>\nhttps://m-unlock.gitlab.io/docs/setup/setup.html<br>\n",
            "outputs": [
                {
                    "type": "File",
                    "label": "Filtered forward read",
                    "doc": "Filtered forward read",
                    "outputSource": "#workflow_illumina_quality.cwl/phix_filter/out_forward_reads",
                    "id": "#workflow_illumina_quality.cwl/QC_forward_reads"
                },
                {
                    "type": "File",
                    "label": "Filtered reverse read",
                    "doc": "Filtered reverse read",
                    "outputSource": "#workflow_illumina_quality.cwl/phix_filter/out_reverse_reads",
                    "id": "#workflow_illumina_quality.cwl/QC_reverse_reads"
                },
                {
                    "type": "Directory",
                    "label": "Filtering reports folder",
                    "doc": "Folder containing all reports of filtering and quality control",
                    "outputSource": "#workflow_illumina_quality.cwl/reports_files_to_folder/results",
                    "id": "#workflow_illumina_quality.cwl/reports_folder"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Remove exact duplicate reads with fastp",
                    "label": "Deduplicate reads",
                    "default": false,
                    "id": "#workflow_illumina_quality.cwl/deduplicate"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output Destination",
                    "doc": "Optional output destination only used for cwl-prov reporting.",
                    "id": "#workflow_illumina_quality.cwl/destination"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "References fasta file(s) for filtering",
                    "label": "Filter reference file(s)",
                    "loadListing": "no_listing",
                    "id": "#workflow_illumina_quality.cwl/filter_references"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Optionally remove rRNA sequences from the reads (default false)",
                    "label": "filter rRNA",
                    "default": false,
                    "id": "#workflow_illumina_quality.cwl/filter_rrna"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Forward sequence fastq file(s) locally",
                    "label": "Forward reads",
                    "loadListing": "no_listing",
                    "id": "#workflow_illumina_quality.cwl/forward_reads"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#workflow_illumina_quality.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Keep with reads mapped to the given reference (default false)",
                    "label": "Keep mapped reads",
                    "default": false,
                    "id": "#workflow_illumina_quality.cwl/keep_reference_mapped_reads"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Kraken2 confidence threshold",
                    "doc": "Confidence score threshold (default 0.0) must be between [0, 1]",
                    "id": "#workflow_illumina_quality.cwl/kraken2_confidence"
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
                    "loadListing": "no_listing",
                    "id": "#workflow_illumina_quality.cwl/kraken2_database"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum memory usage in MegaBytes",
                    "label": "Maximum memory in MB",
                    "default": 4000,
                    "id": "#workflow_illumina_quality.cwl/memory"
                },
                {
                    "type": "boolean",
                    "doc": "Prepare references to a single fasta file and unique headers (default true).\nWhen false a single fasta file as reference is expected with unique headers\n",
                    "label": "Prepare references",
                    "default": true,
                    "id": "#workflow_illumina_quality.cwl/prepare_reference"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Reverse sequence fastq file(s) locally",
                    "label": "Reverse reads",
                    "loadListing": "no_listing",
                    "id": "#workflow_illumina_quality.cwl/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Skip FastQC analyses of raw input data (default false)",
                    "label": "Skip FastQC before",
                    "default": false,
                    "id": "#workflow_illumina_quality.cwl/skip_fastqc_before"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Step number for output folder numbering (default 1)",
                    "label": "Output Step number",
                    "default": 1,
                    "id": "#workflow_illumina_quality.cwl/step"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of threads to use for computational processes",
                    "label": "Number of threads",
                    "default": 2,
                    "id": "#workflow_illumina_quality.cwl/threads"
                }
            ],
            "steps": [
                {
                    "label": "fastp",
                    "doc": "Read quality filtering and (barcode) trimming.",
                    "run": "#fastp.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/deduplicate",
                            "id": "#workflow_illumina_quality.cwl/fastp/deduplicate"
                        },
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/fastq_merge_fwd/output",
                                "#workflow_illumina_quality.cwl/fastq_fwd_array_to_file/file"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_illumina_quality.cwl/fastp/forward_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/identifier",
                            "id": "#workflow_illumina_quality.cwl/fastp/identifier"
                        },
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/fastq_merge_rev/output",
                                "#workflow_illumina_quality.cwl/fastq_rev_array_to_file/file"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_illumina_quality.cwl/fastp/reverse_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/threads",
                            "id": "#workflow_illumina_quality.cwl/fastp/threads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/fastp/out_forward_reads",
                        "#workflow_illumina_quality.cwl/fastp/out_reverse_reads",
                        "#workflow_illumina_quality.cwl/fastp/html_report",
                        "#workflow_illumina_quality.cwl/fastp/json_report"
                    ],
                    "id": "#workflow_illumina_quality.cwl/fastp"
                },
                {
                    "label": "Fwd reads array to file",
                    "doc": "Forward file of single file array to file object",
                    "when": "$(inputs.forward_reads.length === 1)",
                    "run": "#array_to_file.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/forward_reads",
                            "id": "#workflow_illumina_quality.cwl/fastq_fwd_array_to_file/files"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/forward_reads",
                            "id": "#workflow_illumina_quality.cwl/fastq_fwd_array_to_file/forward_reads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/fastq_fwd_array_to_file/file"
                    ],
                    "id": "#workflow_illumina_quality.cwl/fastq_fwd_array_to_file"
                },
                {
                    "label": "Merge forward reads",
                    "doc": "Merge multiple forward fastq reads to a single file",
                    "when": "$(inputs.forward_reads.length > 1)",
                    "run": "#concatenate.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/forward_reads",
                            "id": "#workflow_illumina_quality.cwl/fastq_merge_fwd/forward_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/forward_reads",
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_illumina_quality.cwl/fastq_merge_fwd/infiles"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/identifier",
                            "valueFrom": "$(self)_illumina_merged_1.fq.gz",
                            "id": "#workflow_illumina_quality.cwl/fastq_merge_fwd/outname"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/fastq_merge_fwd/output"
                    ],
                    "id": "#workflow_illumina_quality.cwl/fastq_merge_fwd"
                },
                {
                    "label": "Merge reverse reads",
                    "doc": "Merge multiple reverse fastq reads to a single file",
                    "when": "$(inputs.reverse_reads.length > 1)",
                    "run": "#concatenate.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/reverse_reads",
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_illumina_quality.cwl/fastq_merge_rev/infiles"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/identifier",
                            "valueFrom": "$(self)_illumina_merged_2.fq.gz",
                            "id": "#workflow_illumina_quality.cwl/fastq_merge_rev/outname"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/reverse_reads",
                            "id": "#workflow_illumina_quality.cwl/fastq_merge_rev/reverse_reads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/fastq_merge_rev/output"
                    ],
                    "id": "#workflow_illumina_quality.cwl/fastq_merge_rev"
                },
                {
                    "label": "Rev reads array to file",
                    "doc": "Forward file of single file array to file object",
                    "when": "$(inputs.reverse_reads.length === 1)",
                    "run": "#array_to_file.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/reverse_reads",
                            "id": "#workflow_illumina_quality.cwl/fastq_rev_array_to_file/files"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/reverse_reads",
                            "id": "#workflow_illumina_quality.cwl/fastq_rev_array_to_file/reverse_reads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/fastq_rev_array_to_file/file"
                    ],
                    "id": "#workflow_illumina_quality.cwl/fastq_rev_array_to_file"
                },
                {
                    "label": "FastQC after",
                    "doc": "Quality assessment and report of reads",
                    "run": "#fastqc.cwl",
                    "in": [
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/phix_filter/out_forward_reads",
                                "#workflow_illumina_quality.cwl/phix_filter/out_reverse_reads"
                            ],
                            "id": "#workflow_illumina_quality.cwl/fastqc_illumina_after/fastq"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/threads",
                            "id": "#workflow_illumina_quality.cwl/fastqc_illumina_after/threads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/fastqc_illumina_after/html_files",
                        "#workflow_illumina_quality.cwl/fastqc_illumina_after/zip_files"
                    ],
                    "id": "#workflow_illumina_quality.cwl/fastqc_illumina_after"
                },
                {
                    "label": "FastQC before",
                    "doc": "Quality assessment and report of reads",
                    "run": "#fastqc.cwl",
                    "when": "$(inputs.skip_fastqc_before == false)",
                    "in": [
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/forward_reads",
                                "#workflow_illumina_quality.cwl/reverse_reads"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_illumina_quality.cwl/fastqc_illumina_before/fastq"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/skip_fastqc_before",
                            "id": "#workflow_illumina_quality.cwl/fastqc_illumina_before/skip_fastqc_before"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/threads",
                            "id": "#workflow_illumina_quality.cwl/fastqc_illumina_before/threads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/fastqc_illumina_before/html_files",
                        "#workflow_illumina_quality.cwl/fastqc_illumina_before/zip_files"
                    ],
                    "id": "#workflow_illumina_quality.cwl/fastqc_illumina_before"
                },
                {
                    "label": "Kraken2",
                    "doc": "Taxonomic classification of FASTQ reads",
                    "when": "$(inputs.database !== null && inputs.database.length !== 0)",
                    "run": "#kraken2.cwl",
                    "scatter": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/database",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/kraken2_confidence",
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/confidence"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/kraken2_database",
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/database"
                        },
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/rrna_filter/out_forward_reads",
                                "#workflow_illumina_quality.cwl/fastp/out_forward_reads"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/forward_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/identifier",
                            "valueFrom": "$(self+\"illumina_quality_filtered\")",
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/identifier"
                        },
                        {
                            "default": true,
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/paired_end"
                        },
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/rrna_filter/out_reverse_reads",
                                "#workflow_illumina_quality.cwl/fastp/out_reverse_reads"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/reverse_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/threads",
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/threads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/illumina_quality_kraken2/sample_report"
                    ],
                    "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2"
                },
                {
                    "label": "Krona",
                    "doc": "Visualization of Kraken2 classification with Krona",
                    "when": "$(inputs.kraken2_database !== null && inputs.kraken2_database.length !== 0)",
                    "run": "#krona.cwl",
                    "scatter": "#workflow_illumina_quality.cwl/illumina_quality_kraken2_krona/kraken",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/illumina_quality_kraken2/sample_report",
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2_krona/kraken"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/kraken2_database",
                            "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2_krona/kraken2_database"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/illumina_quality_kraken2_krona/krona_html"
                    ],
                    "id": "#workflow_illumina_quality.cwl/illumina_quality_kraken2_krona"
                },
                {
                    "label": "PhiX filter (bbduk)",
                    "doc": "Filters illumina spike-in PhiX sequences from reads using bbduk",
                    "run": "#bbduk_filter.cwl",
                    "in": [
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/reference_filter_illumina/out_forward_reads",
                                "#workflow_illumina_quality.cwl/rrna_filter/out_forward_reads",
                                "#workflow_illumina_quality.cwl/fastp/out_forward_reads"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_illumina_quality.cwl/phix_filter/forward_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/identifier",
                            "valueFrom": "$(self+\"_illumina_filtered\")",
                            "id": "#workflow_illumina_quality.cwl/phix_filter/identifier"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/memory",
                            "id": "#workflow_illumina_quality.cwl/phix_filter/memory"
                        },
                        {
                            "valueFrom": "/venv/opt/bbmap-39.01-0/resources/phix174_ill.ref.fa.gz",
                            "id": "#workflow_illumina_quality.cwl/phix_filter/reference"
                        },
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/reference_filter_illumina/out_reverse_reads",
                                "#workflow_illumina_quality.cwl/rrna_filter/out_reverse_reads",
                                "#workflow_illumina_quality.cwl/fastp/out_reverse_reads"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_illumina_quality.cwl/phix_filter/reverse_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/threads",
                            "id": "#workflow_illumina_quality.cwl/phix_filter/threads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/phix_filter/out_forward_reads",
                        "#workflow_illumina_quality.cwl/phix_filter/out_reverse_reads",
                        "#workflow_illumina_quality.cwl/phix_filter/summary",
                        "#workflow_illumina_quality.cwl/phix_filter/stats_file"
                    ],
                    "id": "#workflow_illumina_quality.cwl/phix_filter"
                },
                {
                    "label": "Prepare references",
                    "doc": "Prepare references to a single fasta file and unique headers",
                    "when": "$(inputs.fasta_input !== null && inputs.fasta_input.length !== 0)",
                    "run": "#workflow_prepare_fasta_db.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/filter_references",
                            "id": "#workflow_illumina_quality.cwl/prepare_fasta_db/fasta_input"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/prepare_reference",
                            "id": "#workflow_illumina_quality.cwl/prepare_fasta_db/make_headers_unique"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/identifier",
                            "id": "#workflow_illumina_quality.cwl/prepare_fasta_db/output_name"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/prepare_fasta_db/fasta_db"
                    ],
                    "id": "#workflow_illumina_quality.cwl/prepare_fasta_db"
                },
                {
                    "label": "Reference read mapping",
                    "doc": "Map reads against references using BBMap",
                    "when": "$(inputs.filter_references !== null && inputs.filter_references.length !== 0)",
                    "run": "#bbmap_filter-reads.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/filter_references",
                            "id": "#workflow_illumina_quality.cwl/reference_filter_illumina/filter_references"
                        },
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/rrna_filter/out_forward_reads",
                                "#workflow_illumina_quality.cwl/fastp/out_forward_reads"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_illumina_quality.cwl/reference_filter_illumina/forward_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/identifier",
                            "valueFrom": "$(self+\"_ref-filter\")",
                            "id": "#workflow_illumina_quality.cwl/reference_filter_illumina/identifier"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/memory",
                            "id": "#workflow_illumina_quality.cwl/reference_filter_illumina/memory"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/keep_reference_mapped_reads",
                            "id": "#workflow_illumina_quality.cwl/reference_filter_illumina/output_mapped"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/prepare_fasta_db/fasta_db",
                            "id": "#workflow_illumina_quality.cwl/reference_filter_illumina/reference"
                        },
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/rrna_filter/out_reverse_reads",
                                "#workflow_illumina_quality.cwl/fastp/out_reverse_reads"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_illumina_quality.cwl/reference_filter_illumina/reverse_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/threads",
                            "id": "#workflow_illumina_quality.cwl/reference_filter_illumina/threads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/reference_filter_illumina/out_forward_reads",
                        "#workflow_illumina_quality.cwl/reference_filter_illumina/out_reverse_reads",
                        "#workflow_illumina_quality.cwl/reference_filter_illumina/log",
                        "#workflow_illumina_quality.cwl/reference_filter_illumina/stats",
                        "#workflow_illumina_quality.cwl/reference_filter_illumina/covstats"
                    ],
                    "id": "#workflow_illumina_quality.cwl/reference_filter_illumina"
                },
                {
                    "label": "Reports to folder",
                    "doc": "Preparation of fastp output files to a specific output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/step",
                            "valueFrom": "$(self+\"_Illumina_Read_Quality\")",
                            "id": "#workflow_illumina_quality.cwl/reports_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_illumina_quality.cwl/fastqc_illumina_before/html_files",
                                "#workflow_illumina_quality.cwl/fastqc_illumina_before/zip_files",
                                "#workflow_illumina_quality.cwl/fastqc_illumina_after/html_files",
                                "#workflow_illumina_quality.cwl/fastqc_illumina_after/zip_files",
                                "#workflow_illumina_quality.cwl/fastp/html_report",
                                "#workflow_illumina_quality.cwl/fastp/json_report",
                                "#workflow_illumina_quality.cwl/reference_filter_illumina/stats",
                                "#workflow_illumina_quality.cwl/reference_filter_illumina/covstats",
                                "#workflow_illumina_quality.cwl/reference_filter_illumina/log",
                                "#workflow_illumina_quality.cwl/illumina_quality_kraken2/sample_report",
                                "#workflow_illumina_quality.cwl/illumina_quality_kraken2_krona/krona_html",
                                "#workflow_illumina_quality.cwl/phix_filter/summary",
                                "#workflow_illumina_quality.cwl/phix_filter/stats_file",
                                "#workflow_illumina_quality.cwl/rrna_filter/summary",
                                "#workflow_illumina_quality.cwl/rrna_filter/stats_file"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_illumina_quality.cwl/reports_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/reports_files_to_folder/results"
                    ],
                    "id": "#workflow_illumina_quality.cwl/reports_files_to_folder"
                },
                {
                    "label": "rRNA filter (bbduk)",
                    "doc": "Filters rRNA sequences from reads using bbduk",
                    "when": "$(inputs.filter_rrna)",
                    "run": "#bbduk_filter.cwl",
                    "in": [
                        {
                            "source": "#workflow_illumina_quality.cwl/filter_rrna",
                            "id": "#workflow_illumina_quality.cwl/rrna_filter/filter_rrna"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/fastp/out_forward_reads",
                            "id": "#workflow_illumina_quality.cwl/rrna_filter/forward_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/identifier",
                            "valueFrom": "$(self+\"_rRNA-filter\")",
                            "id": "#workflow_illumina_quality.cwl/rrna_filter/identifier"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/memory",
                            "id": "#workflow_illumina_quality.cwl/rrna_filter/memory"
                        },
                        {
                            "valueFrom": "/venv/opt/bbmap-39.01-0/resources/riboKmers.fa.gz",
                            "id": "#workflow_illumina_quality.cwl/rrna_filter/reference"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/fastp/out_reverse_reads",
                            "id": "#workflow_illumina_quality.cwl/rrna_filter/reverse_reads"
                        },
                        {
                            "source": "#workflow_illumina_quality.cwl/threads",
                            "id": "#workflow_illumina_quality.cwl/rrna_filter/threads"
                        }
                    ],
                    "out": [
                        "#workflow_illumina_quality.cwl/rrna_filter/out_forward_reads",
                        "#workflow_illumina_quality.cwl/rrna_filter/out_reverse_reads",
                        "#workflow_illumina_quality.cwl/rrna_filter/summary",
                        "#workflow_illumina_quality.cwl/rrna_filter/stats_file"
                    ],
                    "id": "#workflow_illumina_quality.cwl/rrna_filter"
                }
            ],
            "id": "#workflow_illumina_quality.cwl",
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
                        "#workflow_longread_quality.cwl/reference_filter_longreads/fastq",
                        "#workflow_longread_quality.cwl/filtlong/output_reads"
                    ],
                    "pickValue": "first_non_null",
                    "id": "#workflow_longread_quality.cwl/filtered_reads"
                },
                {
                    "type": "Directory",
                    "label": "Filtering reports folder",
                    "doc": "Folder containing all reports of filtering and quality control",
                    "outputSource": "#workflow_longread_quality.cwl/reports_files_to_folder/results",
                    "id": "#workflow_longread_quality.cwl/reports_folder"
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
                    "id": "#workflow_longread_quality.cwl/destination"
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
                    "id": "#workflow_longread_quality.cwl/filter_references"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#workflow_longread_quality.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Maximum read length threshold (default 90)",
                    "label": "Maximum read length threshold",
                    "default": 90,
                    "id": "#workflow_longread_quality.cwl/keep_percent"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Keep with reads mapped to the given reference (default false)",
                    "label": "Keep mapped reads",
                    "default": false,
                    "id": "#workflow_longread_quality.cwl/keep_reference_mapped_reads"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Kraken2 confidence threshold",
                    "doc": "Confidence score threshold (default 0.0) must be between [0, 1]",
                    "id": "#workflow_longread_quality.cwl/kraken2_confidence"
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
                    "id": "#workflow_longread_quality.cwl/kraken2_database"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Weight given to the length score (default 10)",
                    "label": "Length weigth",
                    "default": 10,
                    "id": "#workflow_longread_quality.cwl/length_weight"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Long read sequence file locally fastq format",
                    "label": "Long reads",
                    "id": "#workflow_longread_quality.cwl/longreads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum memory usage in megabytes",
                    "label": "Maximum memory in MB",
                    "default": 4000,
                    "id": "#workflow_longread_quality.cwl/memory"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Minimum read length threshold (default 1000)",
                    "label": "Minimum read length",
                    "default": 1000,
                    "id": "#workflow_longread_quality.cwl/minimum_length"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Prepare reference with unique headers (default true)",
                    "label": "Prepare references",
                    "default": true,
                    "id": "#workflow_longread_quality.cwl/prepare_reference"
                },
                {
                    "type": "string",
                    "doc": "Type of read i.e. PacBio or Nanopore. Used for naming output files.",
                    "label": "Read type",
                    "id": "#workflow_longread_quality.cwl/readtype"
                },
                {
                    "type": "boolean",
                    "doc": "Skip FastQC analyses of raw input data (default; false)",
                    "label": "Skip FastQC before",
                    "default": false,
                    "id": "#workflow_longread_quality.cwl/skip_fastqc_before"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "CWL base step number",
                    "doc": "Step number for order of steps",
                    "default": 1,
                    "id": "#workflow_longread_quality.cwl/step"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of threads to use for computational processes",
                    "label": "Number of threads",
                    "default": 2,
                    "id": "#workflow_longread_quality.cwl/threads"
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
                                "#workflow_longread_quality.cwl/reference_filter_longreads/fastq",
                                "#workflow_longread_quality.cwl/filtlong/output_reads"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_longread_quality.cwl/fastqc_longreads_after/nanopore_reads"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/threads",
                            "id": "#workflow_longread_quality.cwl/fastqc_longreads_after/threads"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/fastqc_longreads_after/html_files",
                        "#workflow_longread_quality.cwl/fastqc_longreads_after/zip_files"
                    ],
                    "id": "#workflow_longread_quality.cwl/fastqc_longreads_after"
                },
                {
                    "label": "FastQC before",
                    "doc": "Quality assessment and report of reads before filter",
                    "run": "#fastqc.cwl",
                    "when": "$(inputs.skip_fastqc_before == false)",
                    "in": [
                        {
                            "source": [
                                "#workflow_longread_quality.cwl/merge_longreads_fastq/output",
                                "#workflow_longread_quality.cwl/longreads_array_to_file/file"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_longread_quality.cwl/fastqc_longreads_before/nanopore_reads"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/skip_fastqc_before",
                            "id": "#workflow_longread_quality.cwl/fastqc_longreads_before/skip_fastqc_before"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/threads",
                            "id": "#workflow_longread_quality.cwl/fastqc_longreads_before/threads"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/fastqc_longreads_before/html_files",
                        "#workflow_longread_quality.cwl/fastqc_longreads_before/zip_files"
                    ],
                    "id": "#workflow_longread_quality.cwl/fastqc_longreads_before"
                },
                {
                    "label": "Filtlong",
                    "doc": "Filter longreads on quality and length",
                    "run": "#filtlong.cwl",
                    "in": [
                        {
                            "source": "#workflow_longread_quality.cwl/identifier",
                            "id": "#workflow_longread_quality.cwl/filtlong/identifier"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/keep_percent",
                            "id": "#workflow_longread_quality.cwl/filtlong/keep_percent"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/length_weight",
                            "id": "#workflow_longread_quality.cwl/filtlong/length_weight"
                        },
                        {
                            "source": [
                                "#workflow_longread_quality.cwl/merge_longreads_fastq/output",
                                "#workflow_longread_quality.cwl/longreads_array_to_file/file"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_longread_quality.cwl/filtlong/long_reads"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/minimum_length",
                            "id": "#workflow_longread_quality.cwl/filtlong/minimum_length"
                        },
                        {
                            "valueFrom": "$(inputs.identifier)_$(inputs.readtype)_filtered",
                            "id": "#workflow_longread_quality.cwl/filtlong/output_filename"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/readtype",
                            "id": "#workflow_longread_quality.cwl/filtlong/readtype"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/filtlong/output_reads",
                        "#workflow_longread_quality.cwl/filtlong/log"
                    ],
                    "id": "#workflow_longread_quality.cwl/filtlong"
                },
                {
                    "label": "Array to file",
                    "doc": "Pick first file of longreads when only 1 file is given",
                    "when": "$(inputs.longreads.length === 1)",
                    "run": "#array_to_file.cwl",
                    "in": [
                        {
                            "source": "#workflow_longread_quality.cwl/longreads",
                            "id": "#workflow_longread_quality.cwl/longreads_array_to_file/files"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/longreads",
                            "id": "#workflow_longread_quality.cwl/longreads_array_to_file/longreads"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/longreads_array_to_file/file"
                    ],
                    "id": "#workflow_longread_quality.cwl/longreads_array_to_file"
                },
                {
                    "label": "Kraken2",
                    "doc": "Taxonomic classification of FASTQ reads",
                    "when": "$(inputs.kraken2_database !== null && inputs.kraken2_database.length !== 0)",
                    "run": "#kraken2.cwl",
                    "scatter": "#workflow_longread_quality.cwl/longreads_quality_kraken2/database",
                    "in": [
                        {
                            "source": "#workflow_longread_quality.cwl/kraken2_confidence",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2/confidence"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/kraken2_database",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2/database"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/identifier",
                            "valueFrom": "$(self+\"_\"+inputs.readtype)_unfiltered",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2/identifier"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/kraken2_database",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2/kraken2_database"
                        },
                        {
                            "source": [
                                "#workflow_longread_quality.cwl/merge_longreads_fastq/output",
                                "#workflow_longread_quality.cwl/longreads_array_to_file/file"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2/nanopore_reads"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/readtype",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2/readtype"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/threads",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2/threads"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/longreads_quality_kraken2/sample_report"
                    ],
                    "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2"
                },
                {
                    "label": "Krona",
                    "doc": "Visualization of Kraken2 classification with Krona",
                    "when": "$(inputs.kraken2_database !== null && inputs.kraken2_database.length !== 0)",
                    "run": "#krona.cwl",
                    "scatter": "#workflow_longread_quality.cwl/longreads_quality_kraken2_krona/kraken",
                    "in": [
                        {
                            "source": "#workflow_longread_quality.cwl/longreads_quality_kraken2/sample_report",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2_krona/kraken"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/kraken2_database",
                            "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2_krona/kraken2_database"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/longreads_quality_kraken2_krona/krona_html"
                    ],
                    "id": "#workflow_longread_quality.cwl/longreads_quality_kraken2_krona"
                },
                {
                    "label": "Merge fastq files",
                    "when": "$(inputs.longreads.length > 1)",
                    "run": "#concatenate.cwl",
                    "in": [
                        {
                            "source": "#workflow_longread_quality.cwl/identifier",
                            "id": "#workflow_longread_quality.cwl/merge_longreads_fastq/identifier"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/longreads",
                            "id": "#workflow_longread_quality.cwl/merge_longreads_fastq/infiles"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/longreads",
                            "id": "#workflow_longread_quality.cwl/merge_longreads_fastq/longreads"
                        },
                        {
                            "valueFrom": "$(inputs.identifier)_$(inputs.readtype)_merged_raw.fastq.gz",
                            "id": "#workflow_longread_quality.cwl/merge_longreads_fastq/outname"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/readtype",
                            "id": "#workflow_longread_quality.cwl/merge_longreads_fastq/readtype"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/merge_longreads_fastq/output"
                    ],
                    "id": "#workflow_longread_quality.cwl/merge_longreads_fastq"
                },
                {
                    "label": "Prepare references",
                    "doc": "Prepare references to a single fasta file and unique headers",
                    "when": "$(inputs.fasta_input !== null && inputs.fasta_input.length !== 0)",
                    "run": "#workflow_prepare_fasta_db.cwl",
                    "in": [
                        {
                            "source": "#workflow_longread_quality.cwl/filter_references",
                            "id": "#workflow_longread_quality.cwl/prepare_fasta_db/fasta_input"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/prepare_reference",
                            "id": "#workflow_longread_quality.cwl/prepare_fasta_db/make_headers_unique"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/identifier",
                            "id": "#workflow_longread_quality.cwl/prepare_fasta_db/output_name"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/prepare_fasta_db/fasta_db"
                    ],
                    "id": "#workflow_longread_quality.cwl/prepare_fasta_db"
                },
                {
                    "label": "Reference mapping",
                    "doc": "Removal of contaminated reads using minimap2 mapping",
                    "when": "$(inputs.filter_references !== null && inputs.filter_references.length !== 0)",
                    "run": "#minimap2_to_fastq.cwl",
                    "in": [
                        {
                            "source": "#workflow_longread_quality.cwl/filter_references",
                            "id": "#workflow_longread_quality.cwl/reference_filter_longreads/filter_references"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/identifier",
                            "valueFrom": "$(self+\"_\"+inputs.readtype)",
                            "id": "#workflow_longread_quality.cwl/reference_filter_longreads/identifier"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/keep_reference_mapped_reads",
                            "id": "#workflow_longread_quality.cwl/reference_filter_longreads/output_mapped"
                        },
                        {
                            "default": "map-ont",
                            "id": "#workflow_longread_quality.cwl/reference_filter_longreads/preset"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/filtlong/output_reads",
                            "id": "#workflow_longread_quality.cwl/reference_filter_longreads/reads"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/readtype",
                            "id": "#workflow_longread_quality.cwl/reference_filter_longreads/readtype"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/prepare_fasta_db/fasta_db",
                            "id": "#workflow_longread_quality.cwl/reference_filter_longreads/reference"
                        },
                        {
                            "source": "#workflow_longread_quality.cwl/threads",
                            "id": "#workflow_longread_quality.cwl/reference_filter_longreads/threads"
                        }
                    ],
                    "out": [
                        "#workflow_longread_quality.cwl/reference_filter_longreads/fastq",
                        "#workflow_longread_quality.cwl/reference_filter_longreads/log"
                    ],
                    "id": "#workflow_longread_quality.cwl/reference_filter_longreads"
                },
                {
                    "label": "Reports to folder",
                    "doc": "Preparation of fastp output files to a specific output folder",
                    "in": [
                        {
                            "source": "#workflow_longread_quality.cwl/step",
                            "valueFrom": "$(self+\"_Longreads_Read_Quality\")\n",
                            "id": "#workflow_longread_quality.cwl/reports_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_longread_quality.cwl/fastqc_longreads_before/html_files",
                                "#workflow_longread_quality.cwl/fastqc_longreads_before/zip_files",
                                "#workflow_longread_quality.cwl/filtlong/log",
                                "#workflow_longread_quality.cwl/fastqc_longreads_after/html_files",
                                "#workflow_longread_quality.cwl/fastqc_longreads_after/zip_files",
                                "#workflow_longread_quality.cwl/longreads_quality_kraken2/sample_report",
                                "#workflow_longread_quality.cwl/longreads_quality_kraken2_krona/krona_html",
                                "#workflow_longread_quality.cwl/reference_filter_longreads/log"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_longread_quality.cwl/reports_files_to_folder/files"
                        }
                    ],
                    "run": "#files_to_folder.cwl",
                    "out": [
                        "#workflow_longread_quality.cwl/reports_files_to_folder/results"
                    ],
                    "id": "#workflow_longread_quality.cwl/reports_files_to_folder"
                }
            ],
            "id": "#workflow_longread_quality.cwl",
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
                    "class": "ScatterFeatureRequirement"
                },
                {
                    "class": "StepInputExpressionRequirement"
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "label": "Metagenomic GEM construction from assembly",
            "doc": "Workflow for Metagenomics from bins to metabolic model.<br>\nSummary\n  - Prodigal gene prediction\n  - CarveMe genome scale metabolic model reconstruction\n  - MEMOTE for metabolic model testing\n  - SMETANA Species METabolic interaction ANAlysis\n\nOther UNLOCK workflows on WorkflowHub: https://workflowhub.eu/projects/16/workflows?view=default<br><br>\n\n**All tool CWL files and other workflows can be found here:**<br>\n  Tools: https://gitlab.com/m-unlock/cwl<br>\n  Workflows: https://gitlab.com/m-unlock/cwl/workflows<br>\n\n**How to setup and use an UNLOCK workflow:**<br>\nhttps://m-unlock.gitlab.io/docs/setup/setup.html<br>\n",
            "outputs": [
                {
                    "label": "CarveMe GEMs folder",
                    "doc": "CarveMe metabolic models folder",
                    "type": "Directory",
                    "outputSource": "#workflow_metagenomics_GEM.cwl/carveme_files_to_folder/results",
                    "id": "#workflow_metagenomics_GEM.cwl/carveme_gems_folder"
                },
                {
                    "label": "GEMstats",
                    "doc": "CarveMe GEM statistics",
                    "type": "File",
                    "outputSource": "#workflow_metagenomics_GEM.cwl/gemstats/carveme_GEMstats",
                    "id": "#workflow_metagenomics_GEM.cwl/gemstats_out"
                },
                {
                    "label": "MEMOTE outputs folder",
                    "doc": "MEMOTE outputs folder",
                    "type": "Directory",
                    "outputSource": "#workflow_metagenomics_GEM.cwl/memote_files_to_folder/results",
                    "id": "#workflow_metagenomics_GEM.cwl/memote_folder"
                },
                {
                    "label": "Protein files folder",
                    "doc": "Prodigal predicted proteins (compressed) fasta files",
                    "type": "Directory",
                    "outputSource": "#workflow_metagenomics_GEM.cwl/prodigal_files_to_folder/results",
                    "id": "#workflow_metagenomics_GEM.cwl/protein_fasta_folder"
                },
                {
                    "label": "SMETANA output",
                    "doc": "SMETANA detailed output table",
                    "type": [
                        "null",
                        "File"
                    ],
                    "outputSource": "#workflow_metagenomics_GEM.cwl/smetana/detailed_output_tsv",
                    "id": "#workflow_metagenomics_GEM.cwl/smetana_output"
                }
            ],
            "inputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Bin/genome fasta files",
                    "label": "Genome/bin",
                    "id": "#workflow_metagenomics_GEM.cwl/bins"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output Destination (prov only)",
                    "doc": "Not used in this workflow. Output destination used for cwl-prov reporting only.",
                    "id": "#workflow_metagenomics_GEM.cwl/destination"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Gap fill",
                    "doc": "Gap fill model for given media",
                    "id": "#workflow_metagenomics_GEM.cwl/gapfill"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "Identifier used",
                    "id": "#workflow_metagenomics_GEM.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "File"
                    ],
                    "label": "Media database",
                    "doc": "Media database file",
                    "id": "#workflow_metagenomics_GEM.cwl/mediadb"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Run SMETANA",
                    "doc": "Run SMETANA (Species METabolic interaction ANAlysis)",
                    "default": false,
                    "id": "#workflow_metagenomics_GEM.cwl/run_smetana"
                },
                {
                    "type": "string",
                    "doc": "Solver to be used in MEMOTE and SMETANA (defaul; cplex)",
                    "default": "cplex",
                    "id": "#workflow_metagenomics_GEM.cwl/solver"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of threads to use for computational processes",
                    "label": "number of threads",
                    "default": 2,
                    "id": "#workflow_metagenomics_GEM.cwl/threads"
                }
            ],
            "steps": [
                {
                    "label": "CarveMe",
                    "doc": "Genome-scale metabolic models reconstruction with CarveMe",
                    "run": "#carveme.cwl",
                    "scatter": [
                        "#workflow_metagenomics_GEM.cwl/carveme/protein_file"
                    ],
                    "in": [
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/gapfill",
                            "id": "#workflow_metagenomics_GEM.cwl/carveme/gapfill"
                        },
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/mediadb",
                            "id": "#workflow_metagenomics_GEM.cwl/carveme/mediadb"
                        },
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/prodigal/predicted_proteins_faa",
                            "id": "#workflow_metagenomics_GEM.cwl/carveme/protein_file"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/carveme/carveme_gem"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/carveme"
                },
                {
                    "doc": "Preparation of workflow output files to a specific output folder",
                    "label": "CarveMe GEMs to folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "CarveMe_GEMs",
                            "id": "#workflow_metagenomics_GEM.cwl/carveme_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_GEM.cwl/compress_carveme/outfile"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_GEM.cwl/carveme_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/carveme_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/carveme_files_to_folder"
                },
                {
                    "label": "Compress GEM",
                    "doc": "Compress CarveMe GEM",
                    "run": "#pigz.cwl",
                    "scatter": "#workflow_metagenomics_GEM.cwl/compress_carveme/inputfile",
                    "in": [
                        {
                            "source": [
                                "#workflow_metagenomics_GEM.cwl/carveme/carveme_gem"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_GEM.cwl/compress_carveme/inputfile"
                        },
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/threads",
                            "id": "#workflow_metagenomics_GEM.cwl/compress_carveme/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/compress_carveme/outfile"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/compress_carveme"
                },
                {
                    "label": "Compress proteins",
                    "doc": "Compress prodigal protein files",
                    "run": "#pigz.cwl",
                    "scatter": "#workflow_metagenomics_GEM.cwl/compress_prodigal/inputfile",
                    "in": [
                        {
                            "source": [
                                "#workflow_metagenomics_GEM.cwl/prodigal/predicted_proteins_faa"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_GEM.cwl/compress_prodigal/inputfile"
                        },
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/threads",
                            "id": "#workflow_metagenomics_GEM.cwl/compress_prodigal/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/compress_prodigal/outfile"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/compress_prodigal"
                },
                {
                    "label": "GEM stats",
                    "doc": "CarveMe GEM statistics",
                    "run": "#GEMstats.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/carveme/carveme_gem",
                            "id": "#workflow_metagenomics_GEM.cwl/gemstats/carveme_gems"
                        },
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/identifier",
                            "id": "#workflow_metagenomics_GEM.cwl/gemstats/identifier"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/gemstats/carveme_GEMstats"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/gemstats"
                },
                {
                    "doc": "Preparation of workflow output files to a specific output folder",
                    "label": "MEMOTE output",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "MEMOTE",
                            "id": "#workflow_metagenomics_GEM.cwl/memote_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/report_html",
                                "#workflow_metagenomics_GEM.cwl/memote_run/run_json"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_GEM.cwl/memote_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/memote_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/memote_files_to_folder"
                },
                {
                    "label": "MEMOTE report snapshot",
                    "doc": "Take a snapshot of a model's state and generate a report.",
                    "run": "#memote.cwl",
                    "scatter": [
                        "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/GEM"
                    ],
                    "in": [
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/carveme/carveme_gem",
                            "id": "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/GEM"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/report_snapshot"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/skip_test_find_metabolites_consumed_with_closed_bounds"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/skip_test_find_metabolites_not_consumed_with_open_bounds"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/skip_test_find_metabolites_not_produced_with_open_bounds"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/skip_test_find_metabolites_produced_with_closed_bounds"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/memote_report_snapshot/report_html"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/memote_report_snapshot"
                },
                {
                    "label": "MEMOTE report snapshot",
                    "doc": "MEMOTE run analsis",
                    "run": "#memote.cwl",
                    "scatter": [
                        "#workflow_metagenomics_GEM.cwl/memote_run/GEM"
                    ],
                    "in": [
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/carveme/carveme_gem",
                            "id": "#workflow_metagenomics_GEM.cwl/memote_run/GEM"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_run/run"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_run/skip_test_find_metabolites_consumed_with_closed_bounds"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_run/skip_test_find_metabolites_not_consumed_with_open_bounds"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_run/skip_test_find_metabolites_not_produced_with_open_bounds"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/memote_run/skip_test_find_metabolites_produced_with_closed_bounds"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/memote_run/run_json"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/memote_run"
                },
                {
                    "label": "prodigal",
                    "doc": "prodigal gene/protein prediction",
                    "run": "#prodigal.cwl",
                    "scatter": [
                        "#workflow_metagenomics_GEM.cwl/prodigal/input_fasta"
                    ],
                    "in": [
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/bins",
                            "id": "#workflow_metagenomics_GEM.cwl/prodigal/input_fasta"
                        },
                        {
                            "default": true,
                            "id": "#workflow_metagenomics_GEM.cwl/prodigal/single_mode"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/prodigal/predicted_proteins_faa"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/prodigal"
                },
                {
                    "doc": "Preparation of workflow output files to a specific output folder",
                    "label": "Prodigal proteins to folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "Prodigal_proteins",
                            "id": "#workflow_metagenomics_GEM.cwl/prodigal_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_GEM.cwl/compress_prodigal/outfile"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_GEM.cwl/prodigal_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/prodigal_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/prodigal_files_to_folder"
                },
                {
                    "label": "SMETANA",
                    "doc": "Species METabolic interaction ANAlysis",
                    "when": "$(inputs.run_smetana)",
                    "run": "#smetana.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/carveme/carveme_gem",
                            "id": "#workflow_metagenomics_GEM.cwl/smetana/GEM"
                        },
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/identifier",
                            "id": "#workflow_metagenomics_GEM.cwl/smetana/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/run_smetana",
                            "id": "#workflow_metagenomics_GEM.cwl/smetana/run_smetana"
                        },
                        {
                            "source": "#workflow_metagenomics_GEM.cwl/solver",
                            "id": "#workflow_metagenomics_GEM.cwl/smetana/solver"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_GEM.cwl/smetana/detailed_output_tsv"
                    ],
                    "id": "#workflow_metagenomics_GEM.cwl/smetana"
                }
            ],
            "id": "#workflow_metagenomics_GEM.cwl",
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
            "label": "(Hybrid) Metagenomics workflow",
            "doc": "**Workflow (hybrid) metagenomic assembly and binning  **<br>\n  - Workflow Illumina Quality: https://workflowhub.eu/workflows/336?version=1\t\n    - FastQC (control)\n    - fastp (quality trimming)\n    - kraken2 (taxonomy)\n    - bbmap contamination filter\n  - Workflow Longread Quality:\t\n    - FastQC (control)\n    - filtlong (quality trimming)\n    - kraken2 (taxonomy)\n    - minimap2 contamination filter\n  - Kraken2 taxonomic classification of FASTQ reads\n  - SPAdes/Flye (Assembly)\n  - QUAST (Assembly quality report)\n\n  (optional)\n  - Workflow binnning https://workflowhub.eu/workflows/64?version=11\n    - Metabat2/MaxBin2/SemiBin\n    - DAS Tool\n    - CheckM\n    - BUSCO\n    - GTDB-Tk\n\n  (optional)\n  - Workflow Genome-scale metabolic models https://workflowhub.eu/workflows/372\n    - CarveMe (GEM generation)\n    - MEMOTE (GEM test suite)\n    - SMETANA (Species METabolic interaction ANAlysis)\n\nOther UNLOCK workflows on WorkflowHub: https://workflowhub.eu/projects/16/workflows?view=default<br><br>\n\n**All tool CWL files and other workflows can be found here:**<br>\n  Tools: https://gitlab.com/m-unlock/cwl<br>\n  Workflows: https://gitlab.com/m-unlock/cwl/workflows<br>\n\n**How to setup and use an UNLOCK workflow:**<br>\nhttps://m-unlock.gitlab.io/docs/setup/setup.html<br>\n",
            "outputs": [
                {
                    "label": "Assembly output",
                    "doc": "Output from different assembly steps",
                    "type": "Directory",
                    "outputSource": "#main/assembly_files_to_folder/results",
                    "id": "#main/assembly_output"
                },
                {
                    "label": "Binning output",
                    "doc": "Binning outputfolders",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#main/binning_files_to_folder/results",
                    "id": "#main/binning_output"
                },
                {
                    "label": "Community GEM output",
                    "doc": "Community GEM output folder",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#main/GEM_files_to_folder/results",
                    "id": "#main/gem_output"
                },
                {
                    "label": "Kraken2 reports",
                    "doc": "Kraken2 taxonomic classification reports",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#main/kraken2_files_to_folder/results",
                    "id": "#main/kraken2_output"
                },
                {
                    "label": "Read filtering output",
                    "doc": "Read filtering stats + filtered reads",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#main/readfilter_files_to_folder/results",
                    "id": "#main/read_filtering_output"
                },
                {
                    "label": "Read filtering output",
                    "doc": "Read filtering stats + filtered reads",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#main/keep_readfilter_files_to_folder/results",
                    "id": "#main/read_filtering_output_keep"
                }
            ],
            "inputs": [
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Run binning workflow",
                    "doc": "Run with contig binning workflow",
                    "default": false,
                    "id": "#main/binning"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "BUSCO dataset",
                    "doc": "Path to the BUSCO dataset download location",
                    "loadListing": "no_listing",
                    "id": "#main/busco_data"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Remove exact duplicate reads Illumina reads with fastp",
                    "label": "Deduplicate reads",
                    "default": false,
                    "id": "#main/deduplicate"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output Destination (prov only)",
                    "doc": "Not used in this workflow. Output destination used for cwl-prov reporting only.",
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
                    "doc": "Reference fasta file(s) used for pre-filtering. Can be gzipped (not mixed)",
                    "label": "Reference file(s)",
                    "loadListing": "no_listing",
                    "id": "#main/filter_references"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Genome Size",
                    "doc": "Estimated genome size (for example, 5m or 2.6g)",
                    "id": "#main/genome_size"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "doc": "Directory containing the GTDBTK repository",
                    "label": "gtdbtk data directory",
                    "loadListing": "no_listing",
                    "id": "#main/gtdbtk_data"
                },
                {
                    "type": "string",
                    "label": "Identifier",
                    "doc": "Identifier for this dataset used in this workflow",
                    "id": "#main/identifier"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "doc": "Illumina Forward sequence file(s)",
                    "label": "Forward reads",
                    "loadListing": "no_listing",
                    "id": "#main/illumina_forward_reads"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "Reverse reads",
                    "doc": "Illumina Reverse sequence file(s)",
                    "loadListing": "no_listing",
                    "id": "#main/illumina_reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Keep filtered reads in the final output",
                    "label": "Keep filtered reads",
                    "default": false,
                    "id": "#main/keep_filtered_reads"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "label": "Kraken2 confidence threshold",
                    "doc": "Confidence score threshold (default 0.0) must be in [0, 1]",
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
                    "doc": "Database location of kraken2",
                    "label": "Kraken2 database",
                    "default": [],
                    "loadListing": "no_listing",
                    "id": "#main/kraken2_database"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Keep only this percentage of the best reads (measured by bases) (default 90)",
                    "label": "Keep percentage",
                    "default": 90,
                    "id": "#main/longread_keep_percent"
                },
                {
                    "type": [
                        "null",
                        "float"
                    ],
                    "doc": "Weight given to the length score (default 10)",
                    "label": "Length weigth",
                    "default": 10,
                    "id": "#main/longread_length_weight"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Minimum read length threshold (default 1000)",
                    "label": "Minimum read length",
                    "default": 1000,
                    "id": "#main/longread_minimum_length"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum memory usage in megabytes",
                    "label": "Memory usage (MB)",
                    "default": 4000,
                    "id": "#main/memory"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "default": true,
                    "doc": "Metagenome option for assemblers",
                    "label": "When working with metagenomes",
                    "id": "#main/metagenome"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "Oxford Nanopore reads",
                    "doc": "File(s) with Oxford Nanopore reads in FASTQ format",
                    "loadListing": "no_listing",
                    "id": "#main/nanopore_reads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Only spades assembler",
                    "doc": "Run spades in only assembler mode (without read error correction)",
                    "default": false,
                    "id": "#main/only_assembler_mode_spades"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "ONT Basecalling model for MEDAKA",
                    "doc": "Used in MEDAKA\nBasecalling model used with guppy default r941_min_high. \nAvailable: r941_trans, r941_flip213, r941_flip235, r941_min_fast, r941_min_high, r941_prom_fast, r941_prom_high. (required with ONT data)\n",
                    "id": "#main/ont_basecall_model"
                },
                {
                    "type": [
                        "null",
                        {
                            "type": "array",
                            "items": "File"
                        }
                    ],
                    "label": "PacBio reads",
                    "doc": "File(s) with PacBio reads in FASTQ format",
                    "loadListing": "no_listing",
                    "id": "#main/pacbio_reads"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Pilon fix list",
                    "doc": "A comma-separated list of categories of issues to try to fix:\n  \"snps\": try to fix individual base errors;\n  \"indels\": try to fix small indels;\n  \"gaps\": try to fill gaps;\n  \"local\": try to detect and fix local misassemblies;\n  \"all\": all of the above (default);\n  \"bases\": shorthand for \"snps\" and \"indels\" (for back compatibility);\n  default; snps,gaps,local (conservative)\n",
                    "default": "snps,gaps,local",
                    "id": "#main/pilon_fixlist"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Run GEM workflow",
                    "doc": "Run the community genomescale metabolic models workflow on bins",
                    "default": false,
                    "id": "#main/run_GEM"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Use Flye",
                    "doc": "Run with Flye assembler",
                    "default": false,
                    "id": "#main/run_flye"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Use Medaka",
                    "doc": "Run with Mekada assembly polishing with nanopore reads",
                    "default": false,
                    "id": "#main/run_medaka"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Use Pilon",
                    "doc": "Run with Pilon illumina assembly polishing",
                    "default": false,
                    "id": "#main/run_pilon"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Run SMETANA",
                    "doc": "Run SMETANA (Species METabolic interaction ANAlysis)",
                    "default": false,
                    "id": "#main/run_smetana"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Use SPAdes",
                    "doc": "Run with SPAdes assembler",
                    "default": true,
                    "id": "#main/run_spades"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "SemiBin built-in models (human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/chicken_caecum/global)",
                    "label": "SemiBin Environment",
                    "default": "global",
                    "id": "#main/semibin_environment"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
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
                    "doc": "Number of threads to use for computational processes",
                    "label": "Number of threads",
                    "default": 2,
                    "id": "#main/threads"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Continue with reads mapped to the given reference",
                    "label": "Keep mapped reads",
                    "default": false,
                    "id": "#main/use_reference_mapped_reads"
                }
            ],
            "steps": [
                {
                    "doc": "Preparation of GEM workflow output files and folders to a specific output folder",
                    "label": "GEM workflow output to folder",
                    "when": "$(inputs.binning && inputs.run_GEM)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "source": "#main/binning",
                            "id": "#main/GEM_files_to_folder/binning"
                        },
                        {
                            "valueFrom": "$(\"5_metaGEM\")",
                            "id": "#main/GEM_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/workflow_GEM/smetana_output",
                                "#main/workflow_GEM/gemstats_out"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/GEM_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/workflow_GEM/carveme_gems_folder",
                                "#main/workflow_GEM/protein_fasta_folder",
                                "#main/workflow_GEM/memote_folder"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/GEM_files_to_folder/folders"
                        },
                        {
                            "source": "#main/run_GEM",
                            "id": "#main/GEM_files_to_folder/run_GEM"
                        }
                    ],
                    "out": [
                        "#main/GEM_files_to_folder/results"
                    ],
                    "id": "#main/GEM_files_to_folder"
                },
                {
                    "doc": "Preparation of Flye output files to a specific output folder",
                    "label": "Flye output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"3_Assembly\")",
                            "id": "#main/assembly_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/spades_files_to_folder/results",
                                "#main/flye_files_to_folder/results",
                                "#main/medaka_files_to_folder/results",
                                "#main/pilon_files_to_folder/results"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/assembly_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#main/assembly_files_to_folder/results"
                    ],
                    "id": "#main/assembly_files_to_folder"
                },
                {
                    "label": "BBmap read mapping",
                    "doc": "Illumina read mapping using BBmap on assembled contigs",
                    "when": "$(inputs.binning && inputs.forward_reads !== null && inputs.forward_reads.length !== 0)",
                    "run": "#bbmap.cwl",
                    "in": [
                        {
                            "source": "#main/binning",
                            "id": "#main/bbmap/binning"
                        },
                        {
                            "source": "#main/workflow_quality_illumina/QC_forward_reads",
                            "id": "#main/bbmap/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/bbmap/identifier"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/bbmap/memory"
                        },
                        {
                            "source": [
                                "#main/workflow_pilon/pilon_polished_assembly",
                                "#main/medaka/polished_assembly",
                                "#main/flye/assembly",
                                "#main/spades/scaffolds"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/bbmap/reference"
                        },
                        {
                            "source": "#main/workflow_quality_illumina/QC_reverse_reads",
                            "id": "#main/bbmap/reverse_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/bbmap/threads"
                        }
                    ],
                    "out": [
                        "#main/bbmap/sam",
                        "#main/bbmap/stats",
                        "#main/bbmap/covstats",
                        "#main/bbmap/log"
                    ],
                    "id": "#main/bbmap"
                },
                {
                    "doc": "Preparation of binning output files and folders to a specific output folder",
                    "label": "Binning output to folder",
                    "when": "$(inputs.binning)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "source": "#main/binning",
                            "id": "#main/binning_files_to_folder/binning"
                        },
                        {
                            "valueFrom": "$(\"4_Binning\")",
                            "id": "#main/binning_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/contig_read_counts/contigReadCounts",
                                "#main/workflow_binning/bins_read_stats",
                                "#main/workflow_binning/bins_summary_table",
                                "#main/workflow_binning/eukrep_fasta",
                                "#main/workflow_binning/eukrep_stats_file"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/binning_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/workflow_binning/das_tool_output",
                                "#main/workflow_binning/metabat2_output",
                                "#main/workflow_binning/maxbin2_output",
                                "#main/workflow_binning/semibin_output",
                                "#main/workflow_binning/checkm_output",
                                "#main/workflow_binning/gtdbtk_output",
                                "#main/workflow_binning/busco_output"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/binning_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#main/binning_files_to_folder/results"
                    ],
                    "id": "#main/binning_files_to_folder"
                },
                {
                    "label": "SPAdes compressed",
                    "doc": "Compress the large Spades assembly output files",
                    "when": "$(inputs.run_spades && inputs.forward_reads !== null && inputs.forward_reads.length !== 0)",
                    "run": "#pigz.cwl",
                    "scatter": [
                        "#main/compress_spades/inputfile"
                    ],
                    "scatterMethod": "dotproduct",
                    "in": [
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/compress_spades/forward_reads"
                        },
                        {
                            "source": [
                                "#main/spades/contigs",
                                "#main/spades/scaffolds",
                                "#main/spades/assembly_graph",
                                "#main/spades/contigs_before_rr",
                                "#main/spades/contigs_assembly_paths",
                                "#main/spades/scaffolds_assembly_paths"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/compress_spades/inputfile"
                        },
                        {
                            "source": "#main/run_spades",
                            "id": "#main/compress_spades/run_spades"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/compress_spades/threads"
                        }
                    ],
                    "out": [
                        "#main/compress_spades/outfile"
                    ],
                    "id": "#main/compress_spades"
                },
                {
                    "label": "Samtools idxstats",
                    "doc": "Reports alignment summary statistics",
                    "when": "$(inputs.binning && inputs.forward_reads !== null && inputs.forward_reads.length !== 0)",
                    "run": "#samtools_idxstats.cwl",
                    "in": [
                        {
                            "source": "#main/sam_to_sorted_bam/sortedbam",
                            "id": "#main/contig_read_counts/bam_file"
                        },
                        {
                            "source": "#main/binning",
                            "id": "#main/contig_read_counts/binning"
                        },
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/contig_read_counts/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/contig_read_counts/identifier"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/contig_read_counts/threads"
                        }
                    ],
                    "out": [
                        "#main/contig_read_counts/contigReadCounts"
                    ],
                    "id": "#main/contig_read_counts"
                },
                {
                    "label": "Flye assembly",
                    "doc": "De novo assembly of single-molecule reads with Flye",
                    "when": "$(inputs.run_flye)",
                    "run": "#flye.cwl",
                    "in": [
                        {
                            "source": "#main/genome_size",
                            "id": "#main/flye/genome_size"
                        },
                        {
                            "source": "#main/metagenome",
                            "id": "#main/flye/metagenome"
                        },
                        {
                            "source": "#main/workflow_quality_nanopore/filtered_reads",
                            "id": "#main/flye/nano_raw"
                        },
                        {
                            "source": "#main/workflow_quality_pacbio/filtered_reads",
                            "id": "#main/flye/pacbio_raw"
                        },
                        {
                            "source": "#main/run_flye",
                            "id": "#main/flye/run_flye"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/flye/threads"
                        }
                    ],
                    "out": [
                        "#main/flye/00_assembly",
                        "#main/flye/10_consensus",
                        "#main/flye/20_repeat",
                        "#main/flye/30_contigger",
                        "#main/flye/40_polishing",
                        "#main/flye/assembly",
                        "#main/flye/assembly_info",
                        "#main/flye/flye_log",
                        "#main/flye/params"
                    ],
                    "id": "#main/flye"
                },
                {
                    "doc": "Preparation of Flye output files to a specific output folder",
                    "label": "Flye output folder",
                    "when": "$(inputs.run_flye)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"Flye_Assembly\")",
                            "id": "#main/flye_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/flye/assembly",
                                "#main/flye/assembly_info",
                                "#main/flye/flye_log",
                                "#main/flye/params"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/flye_files_to_folder/files"
                        },
                        {
                            "source": "#main/run_flye",
                            "id": "#main/flye_files_to_folder/run_flye"
                        }
                    ],
                    "out": [
                        "#main/flye_files_to_folder/results"
                    ],
                    "id": "#main/flye_files_to_folder"
                },
                {
                    "label": "Kraken2 Illumina",
                    "doc": "Taxonomic classification of illumina FASTQ reads",
                    "when": "$(inputs.database !== null && inputs.database.length !== 0 && inputs.illumina_raw_forward_reads !== null && inputs.illumina_raw_forward_reads.length !== 0)",
                    "run": "#kraken2.cwl",
                    "scatter": "#main/illumina_kraken2/database",
                    "in": [
                        {
                            "source": "#main/kraken2_confidence",
                            "id": "#main/illumina_kraken2/confidence"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/illumina_kraken2/database"
                        },
                        {
                            "source": "#main/workflow_quality_illumina/QC_forward_reads",
                            "id": "#main/illumina_kraken2/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "valueFrom": "$(self)_Illumina_filtered",
                            "id": "#main/illumina_kraken2/identifier"
                        },
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/illumina_kraken2/illumina_raw_forward_reads"
                        },
                        {
                            "default": true,
                            "id": "#main/illumina_kraken2/paired_end"
                        },
                        {
                            "source": "#main/workflow_quality_illumina/QC_reverse_reads",
                            "id": "#main/illumina_kraken2/reverse_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/illumina_kraken2/threads"
                        }
                    ],
                    "out": [
                        "#main/illumina_kraken2/standard_report",
                        "#main/illumina_kraken2/sample_report"
                    ],
                    "id": "#main/illumina_kraken2"
                },
                {
                    "doc": "Preparation of read filtering output files to a specific output folder",
                    "label": "Read filtering output folder",
                    "when": "$(inputs.keep_filtered_reads)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"1_Read_Filtering\")",
                            "id": "#main/keep_readfilter_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/workflow_quality_nanopore/filtered_reads",
                                "#main/workflow_quality_pacbio/filtered_reads",
                                "#main/workflow_quality_illumina/QC_forward_reads",
                                "#main/workflow_quality_illumina/QC_reverse_reads"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/keep_readfilter_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/workflow_quality_nanopore/reports_folder",
                                "#main/workflow_quality_pacbio/reports_folder",
                                "#main/workflow_quality_illumina/reports_folder"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/keep_readfilter_files_to_folder/folders"
                        },
                        {
                            "source": "#main/keep_filtered_reads",
                            "id": "#main/keep_readfilter_files_to_folder/keep_filtered_reads"
                        }
                    ],
                    "out": [
                        "#main/keep_readfilter_files_to_folder/results"
                    ],
                    "id": "#main/keep_readfilter_files_to_folder"
                },
                {
                    "label": "Compress kraken2",
                    "doc": "Compress large kraken2 report file",
                    "when": "$(inputs.kraken2_database !== null && inputs.kraken2_database.length !== 0)",
                    "run": "#pigz.cwl",
                    "scatter": "#main/kraken2_compress/inputfile",
                    "in": [
                        {
                            "source": [
                                "#main/nanopore_kraken2/standard_report",
                                "#main/pacbio_kraken2/standard_report",
                                "#main/illumina_kraken2/standard_report"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/kraken2_compress/inputfile"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/kraken2_compress/kraken2_database"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/kraken2_compress/threads"
                        }
                    ],
                    "out": [
                        "#main/kraken2_compress/outfile"
                    ],
                    "id": "#main/kraken2_compress"
                },
                {
                    "doc": "Preparation of Kraken2 output files to a specific output folder",
                    "label": "Kraken2 output folder",
                    "when": "$(inputs.kraken2_database !== null && inputs.kraken2_database.length !== 0)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"2_Kraken2_classification\")",
                            "id": "#main/kraken2_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/kraken2_compress/outfile",
                                "#main/kraken2_krona/krona_html",
                                "#main/nanopore_kraken2/sample_report",
                                "#main/illumina_kraken2/sample_report"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/kraken2_files_to_folder/files"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/kraken2_files_to_folder/kraken2_database"
                        }
                    ],
                    "out": [
                        "#main/kraken2_files_to_folder/results"
                    ],
                    "id": "#main/kraken2_files_to_folder"
                },
                {
                    "label": "Krona Kraken2",
                    "doc": "Visualization of kraken2 with Krona",
                    "when": "$(inputs.kraken2_database !== null && inputs.kraken2_database.length !== 0)",
                    "run": "#krona.cwl",
                    "scatter": "#main/kraken2_krona/kraken",
                    "in": [
                        {
                            "source": [
                                "#main/nanopore_kraken2/sample_report",
                                "#main/pacbio_kraken2/sample_report",
                                "#main/illumina_kraken2/sample_report"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/kraken2_krona/kraken"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/kraken2_krona/kraken2_database"
                        }
                    ],
                    "out": [
                        "#main/kraken2_krona/krona_html"
                    ],
                    "id": "#main/kraken2_krona"
                },
                {
                    "label": "Medaka polishing of assembly",
                    "doc": "Medaka for (ont reads) polishing of a assembled genome",
                    "when": "$(inputs.run_medaka && inputs.nanopore_reads !== null && inputs.nanopore_reads.length !== 0)",
                    "run": "#medaka_py.cwl",
                    "in": [
                        {
                            "source": "#main/ont_basecall_model",
                            "id": "#main/medaka/basecall_model"
                        },
                        {
                            "source": [
                                "#main/flye/assembly",
                                "#main/spades/scaffolds"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/medaka/draft_assembly"
                        },
                        {
                            "source": "#main/nanopore_reads",
                            "id": "#main/medaka/nanopore_reads"
                        },
                        {
                            "source": "#main/workflow_quality_nanopore/filtered_reads",
                            "id": "#main/medaka/reads"
                        },
                        {
                            "source": "#main/run_medaka",
                            "id": "#main/medaka/run_medaka"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/medaka/threads"
                        }
                    ],
                    "out": [
                        "#main/medaka/polished_assembly",
                        "#main/medaka/gaps_in_draft_coords"
                    ],
                    "id": "#main/medaka"
                },
                {
                    "doc": "Preparation of Medaka output files to a specific output folder",
                    "label": "Medaka output folder",
                    "when": "$(inputs.run_medaka && inputs.nanopore_reads !== null && inputs.nanopore_reads.length !== 0)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"Medaka_assembly_polishing\")",
                            "id": "#main/medaka_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/medaka/polished_assembly",
                                "#main/medaka/gaps_in_draft_coords"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/medaka_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/metaquast_medaka_files_to_folder/results"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/medaka_files_to_folder/folders"
                        },
                        {
                            "source": "#main/nanopore_reads",
                            "id": "#main/medaka_files_to_folder/nanopore_reads"
                        },
                        {
                            "source": "#main/run_medaka",
                            "id": "#main/medaka_files_to_folder/run_medaka"
                        }
                    ],
                    "out": [
                        "#main/medaka_files_to_folder/results"
                    ],
                    "id": "#main/medaka_files_to_folder"
                },
                {
                    "label": "assembly evaluation",
                    "doc": "evaluation of polished assembly with metaQUAST",
                    "when": "$(inputs.run_medaka && inputs.nanopore_reads !== null && inputs.nanopore_reads !== 0)",
                    "run": "#metaquast.cwl",
                    "in": [
                        {
                            "source": "#main/medaka/polished_assembly",
                            "id": "#main/metaquast_medaka/assembly"
                        },
                        {
                            "source": "#main/nanopore_reads",
                            "id": "#main/metaquast_medaka/nanopore_reads"
                        },
                        {
                            "source": "#main/run_medaka",
                            "id": "#main/metaquast_medaka/run_medaka"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/metaquast_medaka/threads"
                        }
                    ],
                    "out": [
                        "#main/metaquast_medaka/metaquast_outdir",
                        "#main/metaquast_medaka/meta_combined_ref",
                        "#main/metaquast_medaka/meta_icarusDir",
                        "#main/metaquast_medaka/metaquast_krona",
                        "#main/metaquast_medaka/not_aligned",
                        "#main/metaquast_medaka/meta_downloaded_ref",
                        "#main/metaquast_medaka/runs_per_reference",
                        "#main/metaquast_medaka/meta_summary",
                        "#main/metaquast_medaka/meta_icarus",
                        "#main/metaquast_medaka/metaquast_log",
                        "#main/metaquast_medaka/metaquast_report",
                        "#main/metaquast_medaka/basicStats",
                        "#main/metaquast_medaka/quast_icarusDir",
                        "#main/metaquast_medaka/quast_icarusHtml",
                        "#main/metaquast_medaka/quastReport",
                        "#main/metaquast_medaka/quastLog",
                        "#main/metaquast_medaka/transposedReport"
                    ],
                    "id": "#main/metaquast_medaka"
                },
                {
                    "doc": "Preparation of metaQUAST output files to a specific output folder",
                    "label": "Medaka metaQUAST output folder",
                    "when": "$(inputs.run_medaka && inputs.nanopore_reads !== null && inputs.nanopore_reads.length !== 0)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"QUAST_Medaka_assembly_quality\")",
                            "id": "#main/metaquast_medaka_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/metaquast_medaka/metaquast_report",
                                "#main/metaquast_medaka/quastReport"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/metaquast_medaka_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/metaquast_medaka/metaquast_krona",
                                "#main/metaquast_medaka/not_aligned",
                                "#main/metaquast_medaka/runs_per_reference"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/metaquast_medaka_files_to_folder/folders"
                        },
                        {
                            "source": "#main/nanopore_reads",
                            "id": "#main/metaquast_medaka_files_to_folder/nanopore_reads"
                        },
                        {
                            "source": "#main/run_medaka",
                            "id": "#main/metaquast_medaka_files_to_folder/run_medaka"
                        }
                    ],
                    "out": [
                        "#main/metaquast_medaka_files_to_folder/results"
                    ],
                    "id": "#main/metaquast_medaka_files_to_folder"
                },
                {
                    "label": "Illumina assembly evaluation",
                    "doc": "Illumina evaluation of pilon polished assembly with metaQUAST",
                    "when": "$(inputs.run_pilon && inputs.forward_reads !== null && inputs.forward_reads.length !== 0)",
                    "run": "#metaquast.cwl",
                    "in": [
                        {
                            "source": "#main/workflow_pilon/pilon_polished_assembly",
                            "id": "#main/metaquast_pilon/assembly"
                        },
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/metaquast_pilon/forward_reads"
                        },
                        {
                            "source": "#main/run_pilon",
                            "id": "#main/metaquast_pilon/run_pilon"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/metaquast_pilon/threads"
                        }
                    ],
                    "out": [
                        "#main/metaquast_pilon/metaquast_outdir",
                        "#main/metaquast_pilon/meta_combined_ref",
                        "#main/metaquast_pilon/meta_icarusDir",
                        "#main/metaquast_pilon/metaquast_krona",
                        "#main/metaquast_pilon/not_aligned",
                        "#main/metaquast_pilon/meta_downloaded_ref",
                        "#main/metaquast_pilon/runs_per_reference",
                        "#main/metaquast_pilon/meta_summary",
                        "#main/metaquast_pilon/meta_icarus",
                        "#main/metaquast_pilon/metaquast_log",
                        "#main/metaquast_pilon/metaquast_report",
                        "#main/metaquast_pilon/basicStats",
                        "#main/metaquast_pilon/quast_icarusDir",
                        "#main/metaquast_pilon/quast_icarusHtml",
                        "#main/metaquast_pilon/quastReport",
                        "#main/metaquast_pilon/quastLog",
                        "#main/metaquast_pilon/transposedReport"
                    ],
                    "id": "#main/metaquast_pilon"
                },
                {
                    "doc": "Preparation of QUAST output files to a specific output folder",
                    "label": "Illumina metaQUAST output folder",
                    "when": "$(inputs.run_pilon && inputs.illumina_forward_reads !== null && inputs.illumina_forward_reads.length !== 0)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"QUAST_Illumina_polished_assembly_quality\")",
                            "id": "#main/metaquast_pilon_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/metaquast_pilon/metaquast_report",
                                "#main/metaquast_pilon/quastReport"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/metaquast_pilon_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/metaquast_pilon/metaquast_krona",
                                "#main/metaquast_pilon/not_aligned"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/metaquast_pilon_files_to_folder/folders"
                        },
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/metaquast_pilon_files_to_folder/illumina_forward_reads"
                        },
                        {
                            "source": "#main/run_pilon",
                            "id": "#main/metaquast_pilon_files_to_folder/run_pilon"
                        }
                    ],
                    "out": [
                        "#main/metaquast_pilon_files_to_folder/results"
                    ],
                    "id": "#main/metaquast_pilon_files_to_folder"
                },
                {
                    "label": "Kraken2 Nanopore",
                    "doc": "Taxonomic classification of nanopore FASTQ reads",
                    "when": "$(inputs.database !== null && inputs.database.length !== 0 && inputs.nanopore_reads !== null && inputs.nanopore_reads !== 0)",
                    "run": "#kraken2.cwl",
                    "scatter": "#main/nanopore_kraken2/database",
                    "in": [
                        {
                            "source": "#main/kraken2_confidence",
                            "id": "#main/nanopore_kraken2/confidence"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/nanopore_kraken2/database"
                        },
                        {
                            "source": "#main/identifier",
                            "valueFrom": "$(self)_Nanopore_filtered",
                            "id": "#main/nanopore_kraken2/identifier"
                        },
                        {
                            "source": "#main/workflow_quality_nanopore/filtered_reads",
                            "id": "#main/nanopore_kraken2/nanopore_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/nanopore_kraken2/threads"
                        }
                    ],
                    "out": [
                        "#main/nanopore_kraken2/standard_report",
                        "#main/nanopore_kraken2/sample_report"
                    ],
                    "id": "#main/nanopore_kraken2"
                },
                {
                    "label": "Kraken2 PacBio",
                    "doc": "Taxonomic classification of PacBio FASTQ reads",
                    "when": "$(inputs.database !== null && inputs.database.length !== 0 && inputs.nanopore_reads !== null && inputs.nanopore_reads !== 0)",
                    "run": "#kraken2.cwl",
                    "scatter": "#main/pacbio_kraken2/database",
                    "in": [
                        {
                            "source": "#main/kraken2_confidence",
                            "id": "#main/pacbio_kraken2/confidence"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/pacbio_kraken2/database"
                        },
                        {
                            "source": "#main/identifier",
                            "valueFrom": "$(self)_PacBio_filtered",
                            "id": "#main/pacbio_kraken2/identifier"
                        },
                        {
                            "source": "#main/workflow_quality_nanopore/filtered_reads",
                            "id": "#main/pacbio_kraken2/nanopore_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/pacbio_kraken2/threads"
                        }
                    ],
                    "out": [
                        "#main/pacbio_kraken2/standard_report",
                        "#main/pacbio_kraken2/sample_report"
                    ],
                    "id": "#main/pacbio_kraken2"
                },
                {
                    "doc": "Preparation of pilon output files to a specific output folder",
                    "label": "Pilon output folder",
                    "when": "$(inputs.run_pilon && inputs.illumina_forward_reads !== null && inputs.illumina_forward_reads.length !== 0)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"Illumina_polished_assembly\")",
                            "id": "#main/pilon_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/workflow_pilon/vcf",
                                "#main/workflow_pilon/pilon_polished_assembly",
                                "#main/workflow_pilon/log"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/pilon_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/metaquast_pilon_files_to_folder/results"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/pilon_files_to_folder/folders"
                        },
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/pilon_files_to_folder/illumina_forward_reads"
                        },
                        {
                            "source": "#main/run_pilon",
                            "id": "#main/pilon_files_to_folder/run_pilon"
                        }
                    ],
                    "out": [
                        "#main/pilon_files_to_folder/results"
                    ],
                    "id": "#main/pilon_files_to_folder"
                },
                {
                    "label": "Prepare references",
                    "doc": "Prepare references to a single fasta file and unique headers",
                    "when": "$(inputs.fasta_files !== null && inputs.fasta_files.length !== 0)",
                    "run": "#prepare_fasta_db.cwl",
                    "in": [
                        {
                            "source": "#main/filter_references",
                            "id": "#main/prepare_fasta_db/fasta_files"
                        },
                        {
                            "valueFrom": "filter-reference_prepared.fa.gz",
                            "id": "#main/prepare_fasta_db/output_file_name"
                        }
                    ],
                    "out": [
                        "#main/prepare_fasta_db/fasta_db"
                    ],
                    "id": "#main/prepare_fasta_db"
                },
                {
                    "doc": "Preparation of read filtering output files to a specific output folder",
                    "label": "Read filtering output folder",
                    "when": "$(inputs.keep_filtered_reads === false)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"1_Read_Filtering\")",
                            "id": "#main/readfilter_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/workflow_quality_nanopore/reports_folder",
                                "#main/workflow_quality_pacbio/reports_folder",
                                "#main/workflow_quality_illumina/reports_folder"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/readfilter_files_to_folder/folders"
                        },
                        {
                            "source": "#main/keep_filtered_reads",
                            "id": "#main/readfilter_files_to_folder/keep_filtered_reads"
                        }
                    ],
                    "out": [
                        "#main/readfilter_files_to_folder/results"
                    ],
                    "id": "#main/readfilter_files_to_folder"
                },
                {
                    "label": "sam conversion to sorted bam",
                    "doc": "Sam file conversion to a sorted indexed bam file",
                    "when": "$(inputs.binning && inputs.forward_reads !== null && inputs.forward_reads.length !== 0)",
                    "run": "#sam_to_sorted-bam.cwl",
                    "in": [
                        {
                            "source": "#main/binning",
                            "id": "#main/sam_to_sorted_bam/binning"
                        },
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/sam_to_sorted_bam/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/sam_to_sorted_bam/identifier"
                        },
                        {
                            "source": "#main/bbmap/sam",
                            "id": "#main/sam_to_sorted_bam/sam"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/sam_to_sorted_bam/threads"
                        }
                    ],
                    "out": [
                        "#main/sam_to_sorted_bam/sortedbam"
                    ],
                    "id": "#main/sam_to_sorted_bam"
                },
                {
                    "doc": "Genome assembly using SPAdes with illumina and or long reads",
                    "label": "SPAdes assembly",
                    "when": "$(inputs.run_spades && inputs.forward_reads !== null && inputs.forward_reads.length !== 0)",
                    "run": "#spades.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/workflow_quality_illumina/QC_forward_reads"
                            ],
                            "linkMerge": "merge_nested",
                            "id": "#main/spades/forward_reads"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/spades/memory"
                        },
                        {
                            "source": "#main/metagenome",
                            "id": "#main/spades/metagenome"
                        },
                        {
                            "source": "#main/workflow_quality_nanopore/filtered_reads",
                            "valueFrom": "${ var reads = null; if (self !== null) { reads = [self]; } return reads; }",
                            "id": "#main/spades/nanopore_reads"
                        },
                        {
                            "source": "#main/only_assembler_mode_spades",
                            "id": "#main/spades/only_assembler"
                        },
                        {
                            "source": "#main/workflow_quality_pacbio/filtered_reads",
                            "valueFrom": "${ var reads = null; if (self !== null) { reads = [self]; } return reads; }",
                            "id": "#main/spades/pacbio_reads"
                        },
                        {
                            "source": [
                                "#main/workflow_quality_illumina/QC_reverse_reads"
                            ],
                            "linkMerge": "merge_nested",
                            "id": "#main/spades/reverse_reads"
                        },
                        {
                            "source": "#main/run_spades",
                            "id": "#main/spades/run_spades"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/spades/threads"
                        }
                    ],
                    "out": [
                        "#main/spades/contigs",
                        "#main/spades/scaffolds",
                        "#main/spades/assembly_graph",
                        "#main/spades/contigs_assembly_paths",
                        "#main/spades/scaffolds_assembly_paths",
                        "#main/spades/contigs_before_rr",
                        "#main/spades/params",
                        "#main/spades/log",
                        "#main/spades/internal_config",
                        "#main/spades/internal_dataset"
                    ],
                    "id": "#main/spades"
                },
                {
                    "doc": "Preparation of SPAdes output files to a specific output folder",
                    "label": "SPADES output to folder",
                    "when": "$(inputs.run_spades)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"SPAdes_Assembly\")",
                            "id": "#main/spades_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/compress_spades/outfile",
                                "#main/spades/params",
                                "#main/spades/log",
                                "#main/spades/internal_config",
                                "#main/spades/internal_dataset"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/spades_files_to_folder/files"
                        },
                        {
                            "source": "#main/run_spades",
                            "id": "#main/spades_files_to_folder/run_spades"
                        }
                    ],
                    "out": [
                        "#main/spades_files_to_folder/results"
                    ],
                    "id": "#main/spades_files_to_folder"
                },
                {
                    "label": "GEM workflow",
                    "doc": "CarveMe community genomescale metabolic models workflow from bins",
                    "when": "$(inputs.binning && inputs.run_GEM)",
                    "run": "#workflow_metagenomics_GEM.cwl",
                    "in": [
                        {
                            "source": "#main/binning",
                            "id": "#main/workflow_GEM/binning"
                        },
                        {
                            "source": "#main/workflow_binning/bins",
                            "id": "#main/workflow_GEM/bins"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_GEM/identifier"
                        },
                        {
                            "source": "#main/run_GEM",
                            "id": "#main/workflow_GEM/run_GEM"
                        },
                        {
                            "source": "#main/run_smetana",
                            "id": "#main/workflow_GEM/run_smetana"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_GEM/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_GEM/carveme_gems_folder",
                        "#main/workflow_GEM/protein_fasta_folder",
                        "#main/workflow_GEM/memote_folder",
                        "#main/workflow_GEM/smetana_output",
                        "#main/workflow_GEM/gemstats_out"
                    ],
                    "id": "#main/workflow_GEM"
                },
                {
                    "label": "Binning workflow",
                    "doc": "Binning workflow to create bins",
                    "when": "$(inputs.binning && inputs.forward_reads !== null && inputs.forward_reads.length !== 0)",
                    "run": "#workflow_metagenomics_binning.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/workflow_pilon/pilon_polished_assembly",
                                "#main/medaka/polished_assembly",
                                "#main/flye/assembly",
                                "#main/spades/scaffolds"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/workflow_binning/assembly"
                        },
                        {
                            "source": "#main/sam_to_sorted_bam/sortedbam",
                            "id": "#main/workflow_binning/bam_file"
                        },
                        {
                            "source": "#main/binning",
                            "id": "#main/workflow_binning/binning"
                        },
                        {
                            "source": "#main/busco_data",
                            "id": "#main/workflow_binning/busco_data"
                        },
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/workflow_binning/forward_reads"
                        },
                        {
                            "source": "#main/gtdbtk_data",
                            "id": "#main/workflow_binning/gtdbtk_data"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_binning/identifier"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/workflow_binning/memory"
                        },
                        {
                            "source": "#main/semibin_environment",
                            "id": "#main/workflow_binning/semibin_environment"
                        },
                        {
                            "default": 1,
                            "id": "#main/workflow_binning/step"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_binning/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_binning/bins",
                        "#main/workflow_binning/das_tool_output",
                        "#main/workflow_binning/maxbin2_output",
                        "#main/workflow_binning/semibin_output",
                        "#main/workflow_binning/metabat2_output",
                        "#main/workflow_binning/checkm_output",
                        "#main/workflow_binning/gtdbtk_output",
                        "#main/workflow_binning/busco_output",
                        "#main/workflow_binning/bins_summary_table",
                        "#main/workflow_binning/bins_read_stats",
                        "#main/workflow_binning/eukrep_fasta",
                        "#main/workflow_binning/eukrep_stats_file"
                    ],
                    "id": "#main/workflow_binning"
                },
                {
                    "label": "Pilon worklow",
                    "doc": "Illumina reads assembly polishing with Pilon",
                    "when": "$(inputs.run_pilon && inputs.illumina_forward_reads !== null && inputs.illumina_forward_reads.length !== 0)",
                    "run": "#workflow_pilon_mapping.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/medaka/polished_assembly",
                                "#main/flye/assembly",
                                "#main/spades/scaffolds"
                            ],
                            "pickValue": "first_non_null",
                            "id": "#main/workflow_pilon/assembly"
                        },
                        {
                            "source": "#main/pilon_fixlist",
                            "id": "#main/workflow_pilon/fixlist"
                        },
                        {
                            "source": "#main/identifier",
                            "valueFrom": "$(self)_scaffolds",
                            "id": "#main/workflow_pilon/identifier"
                        },
                        {
                            "source": "#main/workflow_quality_illumina/QC_forward_reads",
                            "id": "#main/workflow_pilon/illumina_forward_reads"
                        },
                        {
                            "source": "#main/workflow_quality_illumina/QC_reverse_reads",
                            "id": "#main/workflow_pilon/illumina_reverse_reads"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/workflow_pilon/memory"
                        },
                        {
                            "source": "#main/run_pilon",
                            "id": "#main/workflow_pilon/run_pilon"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_pilon/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_pilon/pilon_polished_assembly",
                        "#main/workflow_pilon/vcf",
                        "#main/workflow_pilon/log"
                    ],
                    "id": "#main/workflow_pilon"
                },
                {
                    "label": "Quality and filtering workflow",
                    "doc": "Quality assessment of illumina reads with rRNA filtering option",
                    "when": "$(inputs.forward_reads !== null && inputs.forward_reads.length !== 0)",
                    "run": "#workflow_illumina_quality.cwl",
                    "in": [
                        {
                            "source": "#main/deduplicate",
                            "id": "#main/workflow_quality_illumina/deduplicate"
                        },
                        {
                            "source": [
                                "#main/prepare_fasta_db/fasta_db"
                            ],
                            "id": "#main/workflow_quality_illumina/filter_references"
                        },
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/workflow_quality_illumina/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_quality_illumina/identifier"
                        },
                        {
                            "source": "#main/use_reference_mapped_reads",
                            "id": "#main/workflow_quality_illumina/keep_reference_mapped_reads"
                        },
                        {
                            "source": "#main/kraken2_confidence",
                            "id": "#main/workflow_quality_illumina/kraken2_confidence"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/workflow_quality_illumina/kraken2_database"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/workflow_quality_illumina/memory"
                        },
                        {
                            "default": false,
                            "id": "#main/workflow_quality_illumina/prepare_reference"
                        },
                        {
                            "source": "#main/illumina_reverse_reads",
                            "id": "#main/workflow_quality_illumina/reverse_reads"
                        },
                        {
                            "source": "#main/skip_fastqc_before",
                            "id": "#main/workflow_quality_illumina/skip_fastqc_before"
                        },
                        {
                            "default": 1,
                            "id": "#main/workflow_quality_illumina/step"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_quality_illumina/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_quality_illumina/QC_reverse_reads",
                        "#main/workflow_quality_illumina/QC_forward_reads",
                        "#main/workflow_quality_illumina/reports_folder"
                    ],
                    "id": "#main/workflow_quality_illumina"
                },
                {
                    "label": "Oxford Nanopore quality and filtering workflow",
                    "doc": "Quality and filtering workflow for Oxford Nanopore reads",
                    "when": "$(inputs.longreads !== null && inputs.longreads !== 0)",
                    "run": "#workflow_longread_quality.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/prepare_fasta_db/fasta_db"
                            ],
                            "id": "#main/workflow_quality_nanopore/filter_references"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_quality_nanopore/identifier"
                        },
                        {
                            "source": "#main/longread_keep_percent",
                            "id": "#main/workflow_quality_nanopore/keep_percent"
                        },
                        {
                            "source": "#main/use_reference_mapped_reads",
                            "id": "#main/workflow_quality_nanopore/keep_reference_mapped_reads"
                        },
                        {
                            "source": "#main/kraken2_confidence",
                            "id": "#main/workflow_quality_nanopore/kraken2_confidence"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/workflow_quality_nanopore/kraken2_database"
                        },
                        {
                            "source": "#main/longread_length_weight",
                            "id": "#main/workflow_quality_nanopore/length_weight"
                        },
                        {
                            "source": "#main/nanopore_reads",
                            "id": "#main/workflow_quality_nanopore/longreads"
                        },
                        {
                            "source": "#main/longread_minimum_length",
                            "id": "#main/workflow_quality_nanopore/minimum_length"
                        },
                        {
                            "default": false,
                            "id": "#main/workflow_quality_nanopore/prepare_reference"
                        },
                        {
                            "default": "Nanopore",
                            "id": "#main/workflow_quality_nanopore/readtype"
                        },
                        {
                            "source": "#main/skip_fastqc_before",
                            "id": "#main/workflow_quality_nanopore/skip_fastqc_before"
                        },
                        {
                            "default": 1,
                            "id": "#main/workflow_quality_nanopore/step"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_quality_nanopore/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_quality_nanopore/filtered_reads",
                        "#main/workflow_quality_nanopore/reports_folder"
                    ],
                    "id": "#main/workflow_quality_nanopore"
                },
                {
                    "label": "PacBio quality and filtering workflow",
                    "doc": "Quality and filtering workflow for PacBio reads",
                    "when": "$(inputs.longreads !== null && inputs.longreads !== 0)",
                    "run": "#workflow_longread_quality.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/prepare_fasta_db/fasta_db"
                            ],
                            "id": "#main/workflow_quality_pacbio/filter_references"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_quality_pacbio/identifier"
                        },
                        {
                            "source": "#main/longread_keep_percent",
                            "id": "#main/workflow_quality_pacbio/keep_percent"
                        },
                        {
                            "source": "#main/use_reference_mapped_reads",
                            "id": "#main/workflow_quality_pacbio/keep_reference_mapped_reads"
                        },
                        {
                            "source": "#main/kraken2_confidence",
                            "id": "#main/workflow_quality_pacbio/kraken2_confidence"
                        },
                        {
                            "source": "#main/kraken2_database",
                            "id": "#main/workflow_quality_pacbio/kraken2_database"
                        },
                        {
                            "source": "#main/longread_length_weight",
                            "id": "#main/workflow_quality_pacbio/length_weight"
                        },
                        {
                            "source": "#main/pacbio_reads",
                            "id": "#main/workflow_quality_pacbio/longreads"
                        },
                        {
                            "source": "#main/longread_minimum_length",
                            "id": "#main/workflow_quality_pacbio/minimum_length"
                        },
                        {
                            "default": false,
                            "id": "#main/workflow_quality_pacbio/prepare_reference"
                        },
                        {
                            "default": "PacBio",
                            "id": "#main/workflow_quality_pacbio/readtype"
                        },
                        {
                            "source": "#main/skip_fastqc_before",
                            "id": "#main/workflow_quality_pacbio/skip_fastqc_before"
                        },
                        {
                            "default": 1,
                            "id": "#main/workflow_quality_pacbio/step"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_quality_pacbio/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_quality_pacbio/filtered_reads",
                        "#main/workflow_quality_pacbio/reports_folder"
                    ],
                    "id": "#main/workflow_quality_pacbio"
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
                    "outputSource": "#workflow_metagenomics_binning.cwl/output_bin_files/bins_out",
                    "id": "#workflow_metagenomics_binning.cwl/bins"
                },
                {
                    "label": "Assembly/Bin read stats",
                    "doc": "General assembly and bin coverage",
                    "type": "File",
                    "outputSource": "#workflow_metagenomics_binning.cwl/bin_readstats/binReadStats",
                    "id": "#workflow_metagenomics_binning.cwl/bins_read_stats"
                },
                {
                    "label": "Bins summary",
                    "doc": "Summary of info about the bins",
                    "type": "File",
                    "outputSource": "#workflow_metagenomics_binning.cwl/bins_summary/bins_summary_table",
                    "id": "#workflow_metagenomics_binning.cwl/bins_summary_table"
                },
                {
                    "label": "BUSCO",
                    "doc": "BUSCO output directory",
                    "type": "Directory",
                    "outputSource": "#workflow_metagenomics_binning.cwl/busco_files_to_folder/results",
                    "id": "#workflow_metagenomics_binning.cwl/busco_output"
                },
                {
                    "label": "CheckM",
                    "doc": "CheckM output directory",
                    "type": "Directory",
                    "outputSource": "#workflow_metagenomics_binning.cwl/checkm_files_to_folder/results",
                    "id": "#workflow_metagenomics_binning.cwl/checkm_output"
                },
                {
                    "label": "DAS Tool",
                    "doc": "DAS Tool output directory",
                    "type": "Directory",
                    "outputSource": "#workflow_metagenomics_binning.cwl/das_tool_files_to_folder/results",
                    "id": "#workflow_metagenomics_binning.cwl/das_tool_output"
                },
                {
                    "label": "EukRep fasta",
                    "doc": "EukRep eukaryotic classified contigs",
                    "type": "File",
                    "outputSource": "#workflow_metagenomics_binning.cwl/eukrep/euk_fasta_out",
                    "id": "#workflow_metagenomics_binning.cwl/eukrep_fasta"
                },
                {
                    "label": "EukRep stats",
                    "doc": "EukRep fasta statistics",
                    "type": "File",
                    "outputSource": "#workflow_metagenomics_binning.cwl/eukrep_stats/output",
                    "id": "#workflow_metagenomics_binning.cwl/eukrep_stats_file"
                },
                {
                    "label": "GTDB-Tk",
                    "doc": "GTDB-Tk output directory",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#workflow_metagenomics_binning.cwl/gtdbtk_files_to_folder/results",
                    "id": "#workflow_metagenomics_binning.cwl/gtdbtk_output"
                },
                {
                    "label": "MaxBin2",
                    "doc": "MaxBin2 output directory",
                    "type": "Directory",
                    "outputSource": "#workflow_metagenomics_binning.cwl/maxbin2_files_to_folder/results",
                    "id": "#workflow_metagenomics_binning.cwl/maxbin2_output"
                },
                {
                    "label": "MetaBAT2",
                    "doc": "MetaBAT2 output directory",
                    "type": "Directory",
                    "outputSource": "#workflow_metagenomics_binning.cwl/metabat2_files_to_folder/results",
                    "id": "#workflow_metagenomics_binning.cwl/metabat2_output"
                },
                {
                    "label": "SemiBin",
                    "doc": "MaxBin2 output directory",
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "outputSource": "#workflow_metagenomics_binning.cwl/semibin_files_to_folder/results",
                    "id": "#workflow_metagenomics_binning.cwl/semibin_output"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Assembly in fasta format",
                    "label": "Assembly fasta",
                    "loadListing": "no_listing",
                    "id": "#workflow_metagenomics_binning.cwl/assembly"
                },
                {
                    "type": "File",
                    "doc": "Mapping file in sorted bam format containing reads mapped to the assembly",
                    "label": "Bam file",
                    "loadListing": "no_listing",
                    "id": "#workflow_metagenomics_binning.cwl/bam_file"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "label": "BUSCO dataset",
                    "doc": "Directory containing the BUSCO dataset location.",
                    "loadListing": "no_listing",
                    "id": "#workflow_metagenomics_binning.cwl/busco_data"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output destination (not used in the workflow itself)",
                    "doc": "Optional output destination path for cwl-prov reporting.",
                    "id": "#workflow_metagenomics_binning.cwl/destination"
                },
                {
                    "type": [
                        "null",
                        "Directory"
                    ],
                    "doc": "Directory containing the GTDB database. When none is given GTDB-Tk will be skipped.",
                    "label": "gtdbtk data directory",
                    "loadListing": "no_listing",
                    "id": "#workflow_metagenomics_binning.cwl/gtdbtk_data"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "Identifier used",
                    "id": "#workflow_metagenomics_binning.cwl/identifier"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum memory usage in megabytes",
                    "label": "memory usage (MB)",
                    "default": 4000,
                    "id": "#workflow_metagenomics_binning.cwl/memory"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Run with SemiBin binner",
                    "label": "Run SemiBin",
                    "default": true,
                    "id": "#workflow_metagenomics_binning.cwl/run_semibin"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "doc": "Semibin Built-in models (human_gut/dog_gut/ocean/soil/cat_gut/human_oral/mouse_gut/pig_gut/built_environment/wastewater/global/chicken_caecum)",
                    "label": "SemiBin Environment",
                    "default": "global",
                    "id": "#workflow_metagenomics_binning.cwl/semibin_environment"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "label": "CWL base step number",
                    "doc": "Step number for order of steps",
                    "default": 1,
                    "id": "#workflow_metagenomics_binning.cwl/step"
                },
                {
                    "type": "boolean",
                    "label": "Sub workflow Run",
                    "doc": "Use this when you need the output bins as File[] for subsequent analysis workflow steps in another workflow.",
                    "default": false,
                    "id": "#workflow_metagenomics_binning.cwl/sub_workflow"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of threads to use for computational processes",
                    "label": "Threads",
                    "default": 2,
                    "id": "#workflow_metagenomics_binning.cwl/threads"
                }
            ],
            "steps": [
                {
                    "doc": "Depths per bin",
                    "label": "Depths per bin",
                    "run": "#aggregateBinDepths.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/das_tool_bins/files",
                            "id": "#workflow_metagenomics_binning.cwl/aggregate_bin_depths/bins"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/aggregate_bin_depths/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/metabat2_contig_depths/depths",
                            "id": "#workflow_metagenomics_binning.cwl/aggregate_bin_depths/metabatdepthsFile"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/aggregate_bin_depths/binDepths"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/aggregate_bin_depths"
                },
                {
                    "doc": "Table general bin and assembly read mapping stats",
                    "label": "Bin and assembly read stats",
                    "run": "#assembly_bins_readstats.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/assembly",
                            "id": "#workflow_metagenomics_binning.cwl/bin_readstats/assembly"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/bam_file",
                            "id": "#workflow_metagenomics_binning.cwl/bin_readstats/bam_file"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/das_tool/contig2bin",
                            "id": "#workflow_metagenomics_binning.cwl/bin_readstats/binContigs"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/bin_readstats/identifier"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/bin_readstats/binReadStats"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/bin_readstats"
                },
                {
                    "doc": "Table of all bins and their statistics like size, contigs, completeness etc",
                    "label": "Bins summary",
                    "run": "#bins_summary.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/aggregate_bin_depths/binDepths",
                            "id": "#workflow_metagenomics_binning.cwl/bins_summary/bin_depths"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/das_tool/bin_dir",
                            "id": "#workflow_metagenomics_binning.cwl/bins_summary/bin_dir"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/busco/batch_summary",
                            "id": "#workflow_metagenomics_binning.cwl/bins_summary/busco_batch"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/checkm/checkm_out_table",
                            "id": "#workflow_metagenomics_binning.cwl/bins_summary/checkm"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/gtdbtk/gtdbtk_summary",
                            "id": "#workflow_metagenomics_binning.cwl/bins_summary/gtdbtk"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/bins_summary/identifier"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/bins_summary/bins_summary_table"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/bins_summary"
                },
                {
                    "doc": "BUSCO assembly completeness workflow",
                    "label": "BUSCO",
                    "run": "#busco.cwl",
                    "when": "$(inputs.bins.length !== 0)",
                    "in": [
                        {
                            "valueFrom": "$(true)",
                            "id": "#workflow_metagenomics_binning.cwl/busco/auto-lineage-prok"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/das_tool_bins/files",
                            "id": "#workflow_metagenomics_binning.cwl/busco/bins"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/busco_data",
                            "id": "#workflow_metagenomics_binning.cwl/busco/busco_data"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/busco/identifier"
                        },
                        {
                            "valueFrom": "geno",
                            "id": "#workflow_metagenomics_binning.cwl/busco/mode"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/remove_unbinned/bins_dir",
                            "id": "#workflow_metagenomics_binning.cwl/busco/sequence_folder"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/threads",
                            "id": "#workflow_metagenomics_binning.cwl/busco/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/busco/batch_summary"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/busco"
                },
                {
                    "doc": "Preparation of BUSCO output files to a specific output folder",
                    "label": "BUSCO output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "BUSCO_Bin_Completeness",
                            "id": "#workflow_metagenomics_binning.cwl/busco_files_to_folder/destination"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/busco/batch_summary",
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_metagenomics_binning.cwl/busco_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/busco_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/busco_files_to_folder"
                },
                {
                    "doc": "CheckM bin quality assessment",
                    "label": "CheckM",
                    "run": "#checkm_lineagewf.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/remove_unbinned/bins_dir",
                            "id": "#workflow_metagenomics_binning.cwl/checkm/bin_dir"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/checkm/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/threads",
                            "id": "#workflow_metagenomics_binning.cwl/checkm/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/checkm/checkm_out_table",
                        "#workflow_metagenomics_binning.cwl/checkm/checkm_out_folder"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/checkm"
                },
                {
                    "doc": "Preparation of CheckM output files to a specific output folder",
                    "label": "CheckM output",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "CheckM_Bin_Quality",
                            "id": "#workflow_metagenomics_binning.cwl/checkm_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/checkm/checkm_out_table"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_binning.cwl/checkm_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/checkm/checkm_out_folder"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_binning.cwl/checkm_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/checkm_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/checkm_files_to_folder"
                },
                {
                    "doc": "Compress GTDB-Tk output folder",
                    "label": "Compress GTDB-Tk",
                    "when": "$(inputs.gtdbtk_data !== null)",
                    "run": "#compress_directory.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/gtdbtk_data",
                            "id": "#workflow_metagenomics_binning.cwl/compress_gtdbtk/gtdbtk_data"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/gtdbtk/gtdbtk_out_folder",
                            "id": "#workflow_metagenomics_binning.cwl/compress_gtdbtk/indir"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/compress_gtdbtk/outfile"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/compress_gtdbtk"
                },
                {
                    "doc": "DAS Tool",
                    "label": "DAS Tool integrate predictions from multiple binning tools",
                    "run": "#das_tool.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/assembly",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool/assembly"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/metabat2_contig2bin/table",
                                "#workflow_metagenomics_binning.cwl/maxbin2_contig2bin/table",
                                "#workflow_metagenomics_binning.cwl/semibin_contig2bin/table"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool/bin_tables"
                        },
                        {
                            "valueFrom": "${\n  if (inputs.run_semibin) {\n    return \"MetaBAT2,MaxBin2,SemiBin\";\n  } else {\n    return \"MetaBAT2,MaxBin2\";\n  }\n}\n",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool/binner_labels"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/run_semibin",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool/run_semibin"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/threads",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/das_tool/bin_dir",
                        "#workflow_metagenomics_binning.cwl/das_tool/summary",
                        "#workflow_metagenomics_binning.cwl/das_tool/contig2bin",
                        "#workflow_metagenomics_binning.cwl/das_tool/log"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/das_tool"
                },
                {
                    "doc": "DAS Tool bins folder to File array for further analysis",
                    "label": "Bin dir to files[]",
                    "run": "#folder_to_files.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/das_tool/bin_dir",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool_bins/folder"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/das_tool_bins/files"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/das_tool_bins"
                },
                {
                    "doc": "Preparation of DAS Tool output files to a specific output folder.",
                    "label": "DAS Tool output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "Bin_Refinement_DAS_Tool",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/das_tool/log",
                                "#workflow_metagenomics_binning.cwl/das_tool/summary",
                                "#workflow_metagenomics_binning.cwl/das_tool/contig2bin",
                                "#workflow_metagenomics_binning.cwl/aggregate_bin_depths/binDepths"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/das_tool/bin_dir"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_binning.cwl/das_tool_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/das_tool_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/das_tool_files_to_folder"
                },
                {
                    "doc": "EukRep, eukaryotic sequence classification",
                    "label": "EukRep",
                    "run": "#eukrep.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/assembly",
                            "id": "#workflow_metagenomics_binning.cwl/eukrep/assembly"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/eukrep/identifier"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/eukrep/euk_fasta_out"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/eukrep"
                },
                {
                    "doc": "EukRep fasta statistics",
                    "label": "EukRep stats",
                    "run": "#raw_n50.cwl",
                    "in": [
                        {
                            "valueFrom": "$(inputs.tmp_id)_EukRep",
                            "id": "#workflow_metagenomics_binning.cwl/eukrep_stats/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/eukrep/euk_fasta_out",
                            "id": "#workflow_metagenomics_binning.cwl/eukrep_stats/input_fasta"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/eukrep_stats/tmp_id"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/eukrep_stats/output"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/eukrep_stats"
                },
                {
                    "doc": "Taxomic assigment of bins with GTDB-Tk",
                    "label": "GTDBTK",
                    "when": "$(inputs.gtdbtk_data !== null && inputs.bins.length !== 0)",
                    "run": "#gtdbtk_classify_wf.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/remove_unbinned/bins_dir",
                            "id": "#workflow_metagenomics_binning.cwl/gtdbtk/bin_dir"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/das_tool_bins/files",
                            "id": "#workflow_metagenomics_binning.cwl/gtdbtk/bins"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/gtdbtk_data",
                            "id": "#workflow_metagenomics_binning.cwl/gtdbtk/gtdbtk_data"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/gtdbtk/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/threads",
                            "id": "#workflow_metagenomics_binning.cwl/gtdbtk/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/gtdbtk/gtdbtk_summary",
                        "#workflow_metagenomics_binning.cwl/gtdbtk/gtdbtk_out_folder"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/gtdbtk"
                },
                {
                    "doc": "Preparation of GTDB-Tk output files to a specific output folder",
                    "label": "GTBD-Tk output folder",
                    "when": "$(inputs.gtdbtk_data !== null)",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "GTDB-Tk_Bin_Taxonomy",
                            "id": "#workflow_metagenomics_binning.cwl/gtdbtk_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/gtdbtk/gtdbtk_summary",
                                "#workflow_metagenomics_binning.cwl/compress_gtdbtk/outfile"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_metagenomics_binning.cwl/gtdbtk_files_to_folder/files"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/gtdbtk_data",
                            "id": "#workflow_metagenomics_binning.cwl/gtdbtk_files_to_folder/gtdbtk_data"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/gtdbtk_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/gtdbtk_files_to_folder"
                },
                {
                    "doc": "Binning procedure using MaxBin2",
                    "label": "MaxBin2 binning",
                    "run": "#maxbin2.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/metabat2_contig_depths/depths",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2/abundances"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/assembly",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2/contigs"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/threads",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/maxbin2/bins",
                        "#workflow_metagenomics_binning.cwl/maxbin2/summary",
                        "#workflow_metagenomics_binning.cwl/maxbin2/log"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/maxbin2"
                },
                {
                    "label": "MaxBin2 to contig to bins",
                    "doc": "List the contigs and their corresponding bin.",
                    "run": "#fasta_to_contig2bin.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/maxbin2_to_folder/results",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2_contig2bin/bin_folder"
                        },
                        {
                            "valueFrom": "MaxBin2",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2_contig2bin/binner_name"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/maxbin2_contig2bin/table"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/maxbin2_contig2bin"
                },
                {
                    "doc": "Preparation of maxbin2 output files to a specific output folder.",
                    "label": "MaxBin2 output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "Binner_MaxBin2",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/maxbin2/summary",
                                "#workflow_metagenomics_binning.cwl/maxbin2/log",
                                "#workflow_metagenomics_binning.cwl/maxbin2_contig2bin/table"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/maxbin2_to_folder/results"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/maxbin2_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/maxbin2_files_to_folder"
                },
                {
                    "doc": "Create folder with MaxBin2 bins",
                    "label": "MaxBin2 bins to folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "MaxBin2_bins",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2_to_folder/destination"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/maxbin2/bins",
                            "id": "#workflow_metagenomics_binning.cwl/maxbin2_to_folder/files"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/maxbin2_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/maxbin2_to_folder"
                },
                {
                    "doc": "Binning procedure using MetaBAT2",
                    "label": "MetaBAT2 binning",
                    "run": "#metabat2.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/assembly",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2/assembly"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/metabat2_contig_depths/depths",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2/depths"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/threads",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/metabat2/bin_dir",
                        "#workflow_metagenomics_binning.cwl/metabat2/log"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/metabat2"
                },
                {
                    "label": "MetaBAT2 to contig to bins",
                    "doc": "List the contigs and their corresponding bin.",
                    "run": "#fasta_to_contig2bin.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/metabat2_filter_bins/output_folder",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_contig2bin/bin_folder"
                        },
                        {
                            "valueFrom": "MetaBAT2",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_contig2bin/binner_name"
                        },
                        {
                            "valueFrom": "fa",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_contig2bin/extension"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/metabat2_contig2bin/table"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/metabat2_contig2bin"
                },
                {
                    "label": "contig depths",
                    "doc": "MetabatContigDepths to obtain the depth file used in the MetaBat2 and SemiBin binning process",
                    "run": "#metabatContigDepths.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/bam_file",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_contig_depths/bamFile"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_contig_depths/identifier"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/metabat2_contig_depths/depths"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/metabat2_contig_depths"
                },
                {
                    "doc": "Preparation of MetaBAT2 output files + unbinned contigs to a specific output folder",
                    "label": "MetaBAT2 output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "Binner_MetaBAT2",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/metabat2/log",
                                "#workflow_metagenomics_binning.cwl/metabat2_contig_depths/depths",
                                "#workflow_metagenomics_binning.cwl/metabat2_contig2bin/table"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/metabat2/bin_dir"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_files_to_folder/folders"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/metabat2_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/metabat2_files_to_folder"
                },
                {
                    "doc": "Only keep genome bin fasta files (exlude e.g TooShort.fa)",
                    "label": "Keep MetaBAT2 genome bins",
                    "run": "#folder_file_regex.cwl",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/metabat2/bin_dir",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_filter_bins/folder"
                        },
                        {
                            "valueFrom": "MetaBAT2_bins",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_filter_bins/output_folder_name"
                        },
                        {
                            "valueFrom": "bin\\.[0-9]+\\.fa",
                            "id": "#workflow_metagenomics_binning.cwl/metabat2_filter_bins/regex"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/metabat2_filter_bins/output_folder"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/metabat2_filter_bins"
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
                                "id": "#workflow_metagenomics_binning.cwl/output_bin_files/run/bins"
                            }
                        ],
                        "outputs": [
                            {
                                "type": {
                                    "type": "array",
                                    "items": "File"
                                },
                                "id": "#workflow_metagenomics_binning.cwl/output_bin_files/run/bins_out"
                            }
                        ],
                        "expression": "${ return inputs.bins; }"
                    },
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/das_tool_bins/files",
                            "id": "#workflow_metagenomics_binning.cwl/output_bin_files/bins"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/sub_workflow",
                            "id": "#workflow_metagenomics_binning.cwl/output_bin_files/sub_workflow"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/output_bin_files/bins_out"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/output_bin_files"
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
                                "id": "#workflow_metagenomics_binning.cwl/remove_unbinned/run/bins"
                            }
                        ],
                        "outputs": [
                            {
                                "type": "Directory",
                                "id": "#workflow_metagenomics_binning.cwl/remove_unbinned/run/bins_dir"
                            }
                        ],
                        "expression": "${  \n  var regex = new RegExp('.*unbinned.*');\n  var array = [];\n  for (var i = 0; i < inputs.bins.listing.length; i++) {\n    if (!regex.test(inputs.bins.listing[i].location)){\n      array = array.concat(inputs.bins.listing[i]);\n    }\n  }\n  var r = {\n    'bins_dir':\n      { \"class\": \"Directory\",\n        \"basename\": \"DAS_Tool_genome_bins\",\n        \"listing\": array\n      }\n    };\n  return r;\n}\n"
                    },
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/das_tool/bin_dir",
                            "id": "#workflow_metagenomics_binning.cwl/remove_unbinned/bins"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/remove_unbinned/bins_dir"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/remove_unbinned"
                },
                {
                    "doc": "Binning procedure using SemiBin",
                    "label": "Semibin binning",
                    "run": "#semibin_single_easy_bin.cwl",
                    "when": "$(inputs.run_semibin)",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/assembly",
                            "id": "#workflow_metagenomics_binning.cwl/semibin/assembly"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/semibin_environment",
                            "id": "#workflow_metagenomics_binning.cwl/semibin/environment"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/identifier",
                            "id": "#workflow_metagenomics_binning.cwl/semibin/identifier"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/metabat2_contig_depths/depths",
                            "id": "#workflow_metagenomics_binning.cwl/semibin/metabat2_depth_file"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/run_semibin",
                            "id": "#workflow_metagenomics_binning.cwl/semibin/run_semibin"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/threads",
                            "id": "#workflow_metagenomics_binning.cwl/semibin/threads"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/semibin/recluster_bins",
                        "#workflow_metagenomics_binning.cwl/semibin/data",
                        "#workflow_metagenomics_binning.cwl/semibin/data_split",
                        "#workflow_metagenomics_binning.cwl/semibin/model",
                        "#workflow_metagenomics_binning.cwl/semibin/coverage"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/semibin"
                },
                {
                    "label": "SemiBin to contig to bins",
                    "doc": "List the contigs and their corresponding bin.",
                    "run": "#fasta_to_contig2bin.cwl",
                    "when": "$(inputs.run_semibin)",
                    "in": [
                        {
                            "source": "#workflow_metagenomics_binning.cwl/semibin/recluster_bins",
                            "id": "#workflow_metagenomics_binning.cwl/semibin_contig2bin/bin_folder"
                        },
                        {
                            "valueFrom": "SemiBin",
                            "id": "#workflow_metagenomics_binning.cwl/semibin_contig2bin/binner_name"
                        },
                        {
                            "valueFrom": "fa",
                            "id": "#workflow_metagenomics_binning.cwl/semibin_contig2bin/extension"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/run_semibin",
                            "id": "#workflow_metagenomics_binning.cwl/semibin_contig2bin/run_semibin"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/semibin_contig2bin/table"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/semibin_contig2bin"
                },
                {
                    "doc": "Preparation of SemiBin output files to a specific output folder.",
                    "label": "SemiBin output folder",
                    "run": "#files_to_folder.cwl",
                    "when": "$(inputs.run_semibin)",
                    "in": [
                        {
                            "valueFrom": "Binner_SemiBin",
                            "id": "#workflow_metagenomics_binning.cwl/semibin_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/semibin_contig2bin/table",
                                "#workflow_metagenomics_binning.cwl/semibin/data",
                                "#workflow_metagenomics_binning.cwl/semibin/data_split",
                                "#workflow_metagenomics_binning.cwl/semibin/model",
                                "#workflow_metagenomics_binning.cwl/semibin/coverage"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_metagenomics_binning.cwl/semibin_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#workflow_metagenomics_binning.cwl/semibin/recluster_bins"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#workflow_metagenomics_binning.cwl/semibin_files_to_folder/folders"
                        },
                        {
                            "source": "#workflow_metagenomics_binning.cwl/run_semibin",
                            "id": "#workflow_metagenomics_binning.cwl/semibin_files_to_folder/run_semibin"
                        }
                    ],
                    "out": [
                        "#workflow_metagenomics_binning.cwl/semibin_files_to_folder/results"
                    ],
                    "id": "#workflow_metagenomics_binning.cwl/semibin_files_to_folder"
                }
            ],
            "id": "#workflow_metagenomics_binning.cwl",
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
                },
                {
                    "class": "SubworkflowFeatureRequirement"
                }
            ],
            "label": "Metagenomics workflow",
            "doc": "Workflow pilon assembly polishing\nSteps:\n  - BBmap (Read mapping to assembly)\n  - Pilon\n",
            "outputs": [
                {
                    "label": "Pilon log",
                    "doc": "Pilon log",
                    "type": "File",
                    "outputSource": "#workflow_pilon_mapping.cwl/pilon/pilon_log",
                    "id": "#workflow_pilon_mapping.cwl/log"
                },
                {
                    "label": "Polished genome",
                    "type": "File",
                    "outputSource": "#workflow_pilon_mapping.cwl/pilon/pilon_polished_assembly",
                    "id": "#workflow_pilon_mapping.cwl/pilon_polished_assembly"
                },
                {
                    "label": "VCF file",
                    "doc": "Compressed VCF file containing the changes.",
                    "type": "File",
                    "outputSource": "#workflow_pilon_mapping.cwl/vcf_compress/outfile",
                    "id": "#workflow_pilon_mapping.cwl/vcf"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Assembly in fasta format",
                    "label": "Assembly",
                    "id": "#workflow_pilon_mapping.cwl/assembly"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Output Destination",
                    "doc": "Optional Output destination used for cwl-prov reporting.",
                    "id": "#workflow_pilon_mapping.cwl/destination"
                },
                {
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Pilon fix list",
                    "doc": "A comma-separated list of categories of issues to try to fix",
                    "default": "all",
                    "id": "#workflow_pilon_mapping.cwl/fixlist"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#workflow_pilon_mapping.cwl/identifier"
                },
                {
                    "type": "File",
                    "doc": "forward sequence file locally",
                    "label": "forward reads",
                    "id": "#workflow_pilon_mapping.cwl/illumina_forward_reads"
                },
                {
                    "type": "File",
                    "doc": "reverse sequence file locally",
                    "label": "reverse reads",
                    "id": "#workflow_pilon_mapping.cwl/illumina_reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum memory usage in megabytes",
                    "label": "Maximum memory in MB",
                    "default": 40000,
                    "id": "#workflow_pilon_mapping.cwl/memory"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "number of threads to use for computational processes",
                    "label": "number of threads",
                    "default": 2,
                    "id": "#workflow_pilon_mapping.cwl/threads"
                }
            ],
            "steps": [
                {
                    "label": "samtools index",
                    "doc": "Index file generation for sorted bam file",
                    "run": "#samtools_index.cwl",
                    "in": [
                        {
                            "source": "#workflow_pilon_mapping.cwl/sam_to_sorted_bam/sortedbam",
                            "id": "#workflow_pilon_mapping.cwl/bam_index/bam_file"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/threads",
                            "id": "#workflow_pilon_mapping.cwl/bam_index/threads"
                        }
                    ],
                    "out": [
                        "#workflow_pilon_mapping.cwl/bam_index/bam_index"
                    ],
                    "id": "#workflow_pilon_mapping.cwl/bam_index"
                },
                {
                    "label": "CWL hybrid bam/bai file",
                    "run": "#expression_bam_index.cwl",
                    "in": [
                        {
                            "source": "#workflow_pilon_mapping.cwl/sam_to_sorted_bam/sortedbam",
                            "id": "#workflow_pilon_mapping.cwl/expressiontool_bam_index/bam_file"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/bam_index/bam_index",
                            "id": "#workflow_pilon_mapping.cwl/expressiontool_bam_index/bam_index"
                        }
                    ],
                    "out": [
                        "#workflow_pilon_mapping.cwl/expressiontool_bam_index/hybrid_bamindex"
                    ],
                    "id": "#workflow_pilon_mapping.cwl/expressiontool_bam_index"
                },
                {
                    "label": "pilon",
                    "doc": "Pilon draft assembly polishing with the mapped reads",
                    "run": "#pilon.cwl",
                    "in": [
                        {
                            "source": "#workflow_pilon_mapping.cwl/assembly",
                            "id": "#workflow_pilon_mapping.cwl/pilon/assembly"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/expressiontool_bam_index/hybrid_bamindex",
                            "id": "#workflow_pilon_mapping.cwl/pilon/bam_file"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/fixlist",
                            "id": "#workflow_pilon_mapping.cwl/pilon/fixlist"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/identifier",
                            "id": "#workflow_pilon_mapping.cwl/pilon/identifier"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/memory",
                            "id": "#workflow_pilon_mapping.cwl/pilon/memory"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/threads",
                            "id": "#workflow_pilon_mapping.cwl/pilon/threads"
                        }
                    ],
                    "out": [
                        "#workflow_pilon_mapping.cwl/pilon/pilon_polished_assembly",
                        "#workflow_pilon_mapping.cwl/pilon/pilon_vcf",
                        "#workflow_pilon_mapping.cwl/pilon/pilon_log"
                    ],
                    "id": "#workflow_pilon_mapping.cwl/pilon"
                },
                {
                    "label": "BBMap read mapping",
                    "doc": "Illumina read mapping using BBmap on assembled contigs",
                    "run": "#bbmap.cwl",
                    "in": [
                        {
                            "source": "#workflow_pilon_mapping.cwl/illumina_forward_reads",
                            "id": "#workflow_pilon_mapping.cwl/readmapping_pilon/forward_reads"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/identifier",
                            "id": "#workflow_pilon_mapping.cwl/readmapping_pilon/identifier"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/memory",
                            "id": "#workflow_pilon_mapping.cwl/readmapping_pilon/memory"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/assembly",
                            "id": "#workflow_pilon_mapping.cwl/readmapping_pilon/reference"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/illumina_reverse_reads",
                            "id": "#workflow_pilon_mapping.cwl/readmapping_pilon/reverse_reads"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/threads",
                            "id": "#workflow_pilon_mapping.cwl/readmapping_pilon/threads"
                        }
                    ],
                    "out": [
                        "#workflow_pilon_mapping.cwl/readmapping_pilon/sam",
                        "#workflow_pilon_mapping.cwl/readmapping_pilon/stats",
                        "#workflow_pilon_mapping.cwl/readmapping_pilon/covstats",
                        "#workflow_pilon_mapping.cwl/readmapping_pilon/log"
                    ],
                    "id": "#workflow_pilon_mapping.cwl/readmapping_pilon"
                },
                {
                    "label": "sam conversion to sorted bam",
                    "doc": "Sam file conversion to a sorted bam file",
                    "run": "#sam_to_sorted-bam.cwl",
                    "in": [
                        {
                            "source": "#workflow_pilon_mapping.cwl/identifier",
                            "id": "#workflow_pilon_mapping.cwl/sam_to_sorted_bam/identifier"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/readmapping_pilon/sam",
                            "id": "#workflow_pilon_mapping.cwl/sam_to_sorted_bam/sam"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/threads",
                            "id": "#workflow_pilon_mapping.cwl/sam_to_sorted_bam/threads"
                        }
                    ],
                    "out": [
                        "#workflow_pilon_mapping.cwl/sam_to_sorted_bam/sortedbam"
                    ],
                    "id": "#workflow_pilon_mapping.cwl/sam_to_sorted_bam"
                },
                {
                    "run": "#pigz.cwl",
                    "in": [
                        {
                            "source": "#workflow_pilon_mapping.cwl/pilon/pilon_vcf",
                            "id": "#workflow_pilon_mapping.cwl/vcf_compress/inputfile"
                        },
                        {
                            "source": "#workflow_pilon_mapping.cwl/threads",
                            "id": "#workflow_pilon_mapping.cwl/vcf_compress/threads"
                        }
                    ],
                    "out": [
                        "#workflow_pilon_mapping.cwl/vcf_compress/outfile"
                    ],
                    "id": "#workflow_pilon_mapping.cwl/vcf_compress"
                }
            ],
            "id": "#workflow_pilon_mapping.cwl",
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
            "https://schema.org/dateCreated": "2022-04-00",
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
