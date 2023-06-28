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
                    "outputSource": "#main/pilon/pilon_log",
                    "id": "#main/log"
                },
                {
                    "label": "Polished genome",
                    "type": "File",
                    "outputSource": "#main/pilon/pilon_polished_assembly",
                    "id": "#main/pilon_polished_assembly"
                },
                {
                    "label": "VCF file",
                    "doc": "Compressed VCF file containing the changes.",
                    "type": "File",
                    "outputSource": "#main/vcf_compress/outfile",
                    "id": "#main/vcf"
                }
            ],
            "inputs": [
                {
                    "type": "File",
                    "doc": "Assembly in fasta format",
                    "label": "Assembly",
                    "id": "#main/assembly"
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
                    "type": [
                        "null",
                        "string"
                    ],
                    "label": "Pilon fix list",
                    "doc": "A comma-separated list of categories of issues to try to fix",
                    "default": "all",
                    "id": "#main/fixlist"
                },
                {
                    "type": "string",
                    "doc": "Identifier for this dataset used in this workflow",
                    "label": "identifier used",
                    "id": "#main/identifier"
                },
                {
                    "type": "File",
                    "doc": "forward sequence file locally",
                    "label": "forward reads",
                    "id": "#main/illumina_forward_reads"
                },
                {
                    "type": "File",
                    "doc": "reverse sequence file locally",
                    "label": "reverse reads",
                    "id": "#main/illumina_reverse_reads"
                },
                {
                    "type": [
                        "null",
                        "int"
                    ],
                    "doc": "Maximum memory usage in megabytes",
                    "label": "Maximum memory in MB",
                    "default": 40000,
                    "id": "#main/memory"
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
                    "label": "samtools index",
                    "doc": "Index file generation for sorted bam file",
                    "run": "#samtools_index.cwl",
                    "in": [
                        {
                            "source": "#main/sam_to_sorted_bam/sortedbam",
                            "id": "#main/bam_index/bam_file"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/bam_index/threads"
                        }
                    ],
                    "out": [
                        "#main/bam_index/bam_index"
                    ],
                    "id": "#main/bam_index"
                },
                {
                    "label": "CWL hybrid bam/bai file",
                    "run": "#expression_bam_index.cwl",
                    "in": [
                        {
                            "source": "#main/sam_to_sorted_bam/sortedbam",
                            "id": "#main/expressiontool_bam_index/bam_file"
                        },
                        {
                            "source": "#main/bam_index/bam_index",
                            "id": "#main/expressiontool_bam_index/bam_index"
                        }
                    ],
                    "out": [
                        "#main/expressiontool_bam_index/hybrid_bamindex"
                    ],
                    "id": "#main/expressiontool_bam_index"
                },
                {
                    "label": "pilon",
                    "doc": "Pilon draft assembly polishing with the mapped reads",
                    "run": "#pilon.cwl",
                    "in": [
                        {
                            "source": "#main/assembly",
                            "id": "#main/pilon/assembly"
                        },
                        {
                            "source": "#main/expressiontool_bam_index/hybrid_bamindex",
                            "id": "#main/pilon/bam_file"
                        },
                        {
                            "source": "#main/fixlist",
                            "id": "#main/pilon/fixlist"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/pilon/identifier"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/pilon/memory"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/pilon/threads"
                        }
                    ],
                    "out": [
                        "#main/pilon/pilon_polished_assembly",
                        "#main/pilon/pilon_vcf",
                        "#main/pilon/pilon_log"
                    ],
                    "id": "#main/pilon"
                },
                {
                    "label": "BBMap read mapping",
                    "doc": "Illumina read mapping using BBmap on assembled contigs",
                    "run": "#bbmap.cwl",
                    "in": [
                        {
                            "source": "#main/illumina_forward_reads",
                            "id": "#main/readmapping_pilon/forward_reads"
                        },
                        {
                            "source": "#main/identifier",
                            "id": "#main/readmapping_pilon/identifier"
                        },
                        {
                            "source": "#main/memory",
                            "id": "#main/readmapping_pilon/memory"
                        },
                        {
                            "source": "#main/assembly",
                            "id": "#main/readmapping_pilon/reference"
                        },
                        {
                            "source": "#main/illumina_reverse_reads",
                            "id": "#main/readmapping_pilon/reverse_reads"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/readmapping_pilon/threads"
                        }
                    ],
                    "out": [
                        "#main/readmapping_pilon/sam",
                        "#main/readmapping_pilon/stats",
                        "#main/readmapping_pilon/covstats",
                        "#main/readmapping_pilon/log"
                    ],
                    "id": "#main/readmapping_pilon"
                },
                {
                    "label": "sam conversion to sorted bam",
                    "doc": "Sam file conversion to a sorted bam file",
                    "run": "#sam_to_sorted-bam.cwl",
                    "in": [
                        {
                            "source": "#main/identifier",
                            "id": "#main/sam_to_sorted_bam/identifier"
                        },
                        {
                            "source": "#main/readmapping_pilon/sam",
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
                    "run": "#pigz.cwl",
                    "in": [
                        {
                            "source": "#main/pilon/pilon_vcf",
                            "id": "#main/vcf_compress/inputfile"
                        },
                        {
                            "source": "#main/threads",
                            "id": "#main/vcf_compress/threads"
                        }
                    ],
                    "out": [
                        "#main/vcf_compress/outfile"
                    ],
                    "id": "#main/vcf_compress"
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
            "https://schema.org/dateCreated": "2022-04-00",
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential"
        }
    ],
    "cwlVersion": "v1.2",
    "$namespaces": {
        "s": "https://schema.org/"
    }
}
