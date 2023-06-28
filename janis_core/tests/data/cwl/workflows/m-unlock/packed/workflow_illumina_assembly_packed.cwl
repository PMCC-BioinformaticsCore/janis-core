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
            "baseCommand": [
                "quast.py"
            ],
            "label": "QUAST: Quality Assessment Tool for Genome Assemblies",
            "doc": "Runs the Quality Assessment Tool for Genome Assemblies application\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                },
                {
                    "class": "InitialWorkDirRequirement",
                    "listing": [
                        {
                            "entry": "$({class: 'Directory', listing: []})",
                            "entryname": "QUAST_results",
                            "writable": true
                        },
                        {
                            "entryname": "script.sh",
                            "entry": "#!/bin/bash\nsource /root/miniconda/bin/activate\nconda init bash\nconda activate /unlock/infrastructure/conda/quast/quast_v5.2.0\nquast.py $@"
                        }
                    ]
                }
            ],
            "arguments": [
                {
                    "valueFrom": "QUAST_results",
                    "prefix": "--output-dir"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "The input assembly in fasta format",
                    "label": "assembly fasta file",
                    "inputBinding": {
                        "position": 999
                    },
                    "id": "#quast.cwl/assembly"
                }
            ],
            "outputs": [
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "QUAST_results/basic_stats"
                    },
                    "id": "#quast.cwl/basicStats"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "QUAST_results/icarus_viewers"
                    },
                    "id": "#quast.cwl/icarusDir"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "QUAST_results/icarus.html"
                    },
                    "id": "#quast.cwl/icarusHtml"
                },
                {
                    "type": "File",
                    "outputBinding": {
                        "glob": "QUAST_results/quast.log"
                    },
                    "id": "#quast.cwl/quastLog"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "QUAST_results/report.*"
                    },
                    "id": "#quast.cwl/quastReport"
                },
                {
                    "type": "Directory",
                    "outputBinding": {
                        "glob": "QUAST_results"
                    },
                    "id": "#quast.cwl/quast_outdir"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "QUAST_results/transposed_report.*"
                    },
                    "id": "#quast.cwl/transposedReport"
                }
            ],
            "id": "#quast.cwl",
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
            "label": "Metagenomics workflow",
            "doc": "Workflow for Metagenomics from raw reads to annotated bins.\nSteps:\n  - workflow_illumina_quality.cwl:\n      - FastQC (control)\n      - fastp (quality trimming)\n      - bbmap contamination filter\n  - SPAdes isolate (Assembly)\n  - Pilon (Assembly polishing)\n  - QUAST (Assembly quality report)\n  - BUSCO (Assembly completeness)\n",
            "outputs": [
                {
                    "label": "BUSCO",
                    "doc": "BUSCO analysis output folder",
                    "type": "Directory",
                    "outputSource": "#main/busco_files_to_folder/results",
                    "id": "#main/busco_output"
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
                    "label": "Filtered statistics",
                    "doc": "Statistics on quality and preprocessing of the reads",
                    "type": "Directory",
                    "outputSource": "#main/workflow_quality/reports_folder",
                    "id": "#main/filtered_stats"
                },
                {
                    "type": "Directory",
                    "outputSource": "#main/pilon_files_to_folder/results",
                    "id": "#main/pilon_output"
                },
                {
                    "label": "QUAST",
                    "doc": "Quast analysis output folder",
                    "type": "Directory",
                    "outputSource": "#main/quast_files_to_folder/results",
                    "id": "#main/quast_output"
                },
                {
                    "label": "SPAdes",
                    "doc": "Metagenome assembly output by SPADES",
                    "type": "Directory",
                    "outputSource": "#main/spades_files_to_folder/results",
                    "id": "#main/spades_output"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "default": "/unlock/references/databases/BUSCO/BUSCO_odb10",
                    "label": "BUSCO dataset",
                    "doc": "Path to the BUSCO dataset download location",
                    "id": "#main/busco_dataset"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "default": false,
                    "label": "BUSCO eukaryote",
                    "doc": "run BUSCO with auto eukaryote lineage",
                    "id": "#main/busco_euk"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "default": false,
                    "label": "BUSCO prokaryote",
                    "doc": "Run BUSCO with auto prokaryote lineage",
                    "id": "#main/busco_prok"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "doc": "Remove exact duplicate reads with fastp",
                    "label": "Deduplicate reads",
                    "default": false,
                    "id": "#main/deduplicate"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "string"
                    },
                    "doc": "Reference fasta file for contamination filtering",
                    "label": "Reference files (filter)",
                    "id": "#main/filter_references"
                },
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Forward sequence file locally",
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
                    "type": "boolean",
                    "doc": "Continue with mapped reads the reference",
                    "label": "Keep mapped reads",
                    "default": false,
                    "id": "#main/mapped_reads"
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
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "doc": "Reverse sequence file locally",
                    "label": "reverse reads",
                    "id": "#main/reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Number of threads to use for computational processes",
                    "label": "number of threads",
                    "default": 2,
                    "id": "#main/threads"
                }
            ],
            "steps": [
                {
                    "doc": "Preparation of BUSCO output files to a specific output folder",
                    "label": "BUSCO output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"5_BUSCO_AsssemblyCompleteness\")",
                            "id": "#main/busco_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/workflow_busco/short_summaries"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/busco_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/workflow_busco/run_folders",
                                "#main/workflow_busco/logs"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/busco_files_to_folder/folders"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/busco_files_to_folder/identifier"
                        }
                    ],
                    "out": [
                        "#main/busco_files_to_folder/results"
                    ],
                    "id": "#main/busco_files_to_folder"
                },
                {
                    "doc": "Preparation of quast output files to a specific output folder",
                    "label": "QUAST output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"3_Polished_assembly\")",
                            "id": "#main/pilon_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/workflow_pilon/pilon_vcf",
                                "#main/workflow_pilon/pilon_polished_assembly"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/pilon_files_to_folder/files"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/pilon_files_to_folder/identifier"
                        }
                    ],
                    "out": [
                        "#main/pilon_files_to_folder/results"
                    ],
                    "id": "#main/pilon_files_to_folder"
                },
                {
                    "doc": "Preparation of quast output files to a specific output folder",
                    "label": "QUAST output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"4_Quast_AsssemblyQuality_\")$(inputs.identifier)",
                            "id": "#main/quast_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/workflow_quast/quastLog",
                                "#main/workflow_quast/icarusHtml",
                                "#main/workflow_quast/quastReport",
                                "#main/workflow_quast/transposedReport"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/quast_files_to_folder/files"
                        },
                        {
                            "source": [
                                "#main/workflow_quast/basicStats",
                                "#main/workflow_quast/icarusDir"
                            ],
                            "linkMerge": "merge_flattened",
                            "id": "#main/quast_files_to_folder/folders"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/quast_files_to_folder/identifier"
                        }
                    ],
                    "out": [
                        "#main/quast_files_to_folder/results"
                    ],
                    "id": "#main/quast_files_to_folder"
                },
                {
                    "doc": "Preparation of SPAdes output files to a specific output folder",
                    "label": "SPADES output folder",
                    "run": "#files_to_folder.cwl",
                    "in": [
                        {
                            "valueFrom": "$(\"2_SPAdes_Assembly\")",
                            "id": "#main/spades_files_to_folder/destination"
                        },
                        {
                            "source": [
                                "#main/workflow_spades/scaffolds",
                                "#main/workflow_spades/contigs",
                                "#main/workflow_spades/params",
                                "#main/workflow_spades/log",
                                "#main/workflow_spades/internal_config",
                                "#main/workflow_spades/internal_dataset",
                                "#main/workflow_bbmap/stats",
                                "#main/workflow_bbmap/covstats",
                                "#main/workflow_bbmap/log"
                            ],
                            "linkMerge": "merge_flattened",
                            "pickValue": "all_non_null",
                            "id": "#main/spades_files_to_folder/files"
                        }
                    ],
                    "out": [
                        "#main/spades_files_to_folder/results"
                    ],
                    "id": "#main/spades_files_to_folder"
                },
                {
                    "label": "bbmap read mapping",
                    "doc": "Illumina read mapping using BBmap on assembled contigs",
                    "run": "#bbmap.cwl",
                    "in": [
                        {
                            "source": "#main/workflow_quality/QC_forward_reads",
                            "id": "#main/workflow_bbmap/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_bbmap/identifier"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/workflow_bbmap/memory"
                        },
                        {
                            "source": "#main/workflow_spades/scaffolds",
                            "id": "#main/workflow_bbmap/reference"
                        },
                        {
                            "source": "#main/workflow_quality/QC_reverse_reads",
                            "id": "#main/workflow_bbmap/reverse_reads"
                        }
                    ],
                    "out": [
                        "#main/workflow_bbmap/sam",
                        "#main/workflow_bbmap/stats",
                        "#main/workflow_bbmap/covstats",
                        "#main/workflow_bbmap/log"
                    ],
                    "id": "#main/workflow_bbmap"
                },
                {
                    "doc": "BUSCO assembly completeness workflow",
                    "label": "BUSCO workflow",
                    "run": "#busco.cwl",
                    "in": [
                        {
                            "source": "#main/busco_euk",
                            "id": "#main/workflow_busco/auto-lineage-euk"
                        },
                        {
                            "source": "#main/busco_prok",
                            "id": "#main/workflow_busco/auto-lineage-prok"
                        },
                        {
                            "source": "#main/busco_dataset",
                            "id": "#main/workflow_busco/dataset"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_busco/identifier"
                        },
                        {
                            "valueFrom": "geno",
                            "id": "#main/workflow_busco/mode"
                        },
                        {
                            "source": "#main/workflow_pilon/pilon_polished_assembly",
                            "id": "#main/workflow_busco/sequence_file"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_busco/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_busco/logs",
                        "#main/workflow_busco/run_folders",
                        "#main/workflow_busco/short_summaries"
                    ],
                    "id": "#main/workflow_busco"
                },
                {
                    "run": "#expression_bam_index.cwl",
                    "in": [
                        {
                            "source": "#main/workflow_sam_to_sorted_bam/sortedbam",
                            "id": "#main/workflow_expressiontool_bam_index/bam_file"
                        },
                        {
                            "source": "#main/workflow_sam_index/bam_index",
                            "id": "#main/workflow_expressiontool_bam_index/bam_index"
                        }
                    ],
                    "out": [
                        "#main/workflow_expressiontool_bam_index/hybrid_bamindex"
                    ],
                    "id": "#main/workflow_expressiontool_bam_index"
                },
                {
                    "label": "pilon",
                    "doc": "Pilon draft assembly polishing with the mapped reads",
                    "run": "#pilon.cwl",
                    "in": [
                        {
                            "source": "#main/workflow_spades/scaffolds",
                            "id": "#main/workflow_pilon/assembly"
                        },
                        {
                            "source": "#main/workflow_expressiontool_bam_index/hybrid_bamindex",
                            "id": "#main/workflow_pilon/bam_file"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_pilon/identifier"
                        }
                    ],
                    "out": [
                        "#main/workflow_pilon/pilon_polished_assembly",
                        "#main/workflow_pilon/pilon_vcf"
                    ],
                    "id": "#main/workflow_pilon"
                },
                {
                    "label": "Quality and filtering workflow",
                    "doc": "Quality assessment of illumina reads with rRNA filtering option",
                    "run": "#workflow_illumina_quality.cwl",
                    "in": [
                        {
                            "source": "#main/deduplicate",
                            "id": "#main/workflow_quality/deduplicate"
                        },
                        {
                            "source": "#main/filter_references",
                            "id": "#main/workflow_quality/filter_references"
                        },
                        {
                            "source": "#main/forward_reads",
                            "id": "#main/workflow_quality/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_quality/identifier"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/workflow_quality/memory"
                        },
                        {
                            "source": "#main/mapped_reads",
                            "id": "#main/workflow_quality/output_mapped"
                        },
                        {
                            "source": "#main/reverse_reads",
                            "id": "#main/workflow_quality/reverse_reads"
                        },
                        {
                            "default": 1,
                            "id": "#main/workflow_quality/step"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_quality/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_quality/QC_reverse_reads",
                        "#main/workflow_quality/QC_forward_reads",
                        "#main/workflow_quality/reports_folder"
                    ],
                    "id": "#main/workflow_quality"
                },
                {
                    "doc": "Genome assembly quality assessment using Quast",
                    "label": "Quast workflow",
                    "run": "#quast.cwl",
                    "in": [
                        {
                            "source": "#main/workflow_pilon/pilon_polished_assembly",
                            "id": "#main/workflow_quast/assembly"
                        }
                    ],
                    "out": [
                        "#main/workflow_quast/basicStats",
                        "#main/workflow_quast/icarusDir",
                        "#main/workflow_quast/icarusHtml",
                        "#main/workflow_quast/quastReport",
                        "#main/workflow_quast/quastLog",
                        "#main/workflow_quast/transposedReport"
                    ],
                    "id": "#main/workflow_quast"
                },
                {
                    "label": "samtools index",
                    "doc": "Index file generation for sorted bam file",
                    "run": "#samtools_index.cwl",
                    "in": [
                        {
                            "source": "#main/workflow_sam_to_sorted_bam/sortedbam",
                            "id": "#main/workflow_sam_index/bam_file"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_sam_index/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_sam_index/bam_index"
                    ],
                    "id": "#main/workflow_sam_index"
                },
                {
                    "label": "sam conversion to sorted bam",
                    "doc": "Sam file conversion to a sorted bam file",
                    "run": "#sam_to_sorted-bam.cwl",
                    "in": [
                        {
                            "source": "#main/identifier",
                            "id": "#main/workflow_sam_to_sorted_bam/identifier"
                        },
                        {
                            "source": "#main/workflow_bbmap/sam",
                            "id": "#main/workflow_sam_to_sorted_bam/sam"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_sam_to_sorted_bam/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_sam_to_sorted_bam/sortedbam"
                    ],
                    "id": "#main/workflow_sam_to_sorted_bam"
                },
                {
                    "doc": "Genome assembly using spades with illumina/pacbio reads",
                    "label": "SPAdes assembly",
                    "run": "#spades.cwl",
                    "in": [
                        {
                            "source": [
                                "#main/workflow_quality/QC_forward_reads"
                            ],
                            "linkMerge": "merge_nested",
                            "id": "#main/workflow_spades/forward_reads"
                        },
                        {
                            "default": true,
                            "id": "#main/workflow_spades/isolate"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/workflow_spades/memory"
                        },
                        {
                            "source": [
                                "#main/workflow_quality/QC_reverse_reads"
                            ],
                            "linkMerge": "merge_nested",
                            "id": "#main/workflow_spades/reverse_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/workflow_spades/threads"
                        }
                    ],
                    "out": [
                        "#main/workflow_spades/contigs",
                        "#main/workflow_spades/scaffolds",
                        "#main/workflow_spades/assembly_graph",
                        "#main/workflow_spades/contigs_assembly_paths",
                        "#main/workflow_spades/scaffolds_assembly_paths",
                        "#main/workflow_spades/contigs_before_rr",
                        "#main/workflow_spades/params",
                        "#main/workflow_spades/log",
                        "#main/workflow_spades/internal_config",
                        "#main/workflow_spades/internal_dataset"
                    ],
                    "id": "#main/workflow_spades"
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
            "https://schema.org/dateCreated": "2022-02-00",
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
