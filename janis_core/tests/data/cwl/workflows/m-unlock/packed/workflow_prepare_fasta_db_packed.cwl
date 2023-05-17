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
                    "id": "#main/fasta_input"
                },
                {
                    "type": [
                        "null",
                        "boolean"
                    ],
                    "label": "Make headers unique",
                    "doc": "Make fasta headers unique avoiding same fasta headers, which can break some tools.",
                    "default": false,
                    "id": "#main/make_headers_unique"
                },
                {
                    "type": "string",
                    "doc": "Output name for this dataset used",
                    "label": "identifier used",
                    "id": "#main/output_name"
                }
            ],
            "outputs": [
                {
                    "type": "File",
                    "label": "Prepared fasta file",
                    "doc": "Prepared fasta file",
                    "outputSource": [
                        "#main/fasta_array_to_file/file",
                        "#main/merge_input/output",
                        "#main/prepare_fasta_db/fasta_db"
                    ],
                    "pickValue": "first_non_null",
                    "id": "#main/fasta_db"
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
                            "source": "#main/fasta_input",
                            "id": "#main/fasta_array_to_file/fasta_input"
                        },
                        {
                            "source": "#main/fasta_input",
                            "id": "#main/fasta_array_to_file/files"
                        },
                        {
                            "source": "#main/make_headers_unique",
                            "id": "#main/fasta_array_to_file/make_headers_unique"
                        }
                    ],
                    "out": [
                        "#main/fasta_array_to_file/file"
                    ],
                    "id": "#main/fasta_array_to_file"
                },
                {
                    "label": "Merge reference files",
                    "doc": "Only merge input when make unique is false.",
                    "when": "$(inputs.make_headers_unique === false && inputs.fasta_input.length > 1)",
                    "run": "#concatenate.cwl",
                    "in": [
                        {
                            "source": "#main/fasta_input",
                            "id": "#main/merge_input/fasta_input"
                        },
                        {
                            "source": "#main/fasta_input",
                            "id": "#main/merge_input/infiles"
                        },
                        {
                            "source": "#main/make_headers_unique",
                            "id": "#main/merge_input/make_headers_unique"
                        },
                        {
                            "valueFrom": "$(inputs.output_name)_filter-reference_merged.fa",
                            "id": "#main/merge_input/outname"
                        },
                        {
                            "source": "#main/output_name",
                            "id": "#main/merge_input/output_name"
                        }
                    ],
                    "out": [
                        "#main/merge_input/output"
                    ],
                    "id": "#main/merge_input"
                },
                {
                    "label": "Prepare references",
                    "doc": "Prepare references to a single fasta file and unique headers",
                    "when": "$(inputs.make_headers_unique)",
                    "run": "#prepare_fasta_db.cwl",
                    "in": [
                        {
                            "source": "#main/fasta_input",
                            "id": "#main/prepare_fasta_db/fasta_files"
                        },
                        {
                            "source": "#main/fasta_input",
                            "id": "#main/prepare_fasta_db/fasta_input"
                        },
                        {
                            "source": "#main/make_headers_unique",
                            "id": "#main/prepare_fasta_db/make_headers_unique"
                        },
                        {
                            "valueFrom": "$(inputs.output_name)_filter-reference_uniq.fa.gz",
                            "id": "#main/prepare_fasta_db/output_file_name"
                        },
                        {
                            "source": "#main/output_name",
                            "id": "#main/prepare_fasta_db/output_name"
                        }
                    ],
                    "out": [
                        "#main/prepare_fasta_db/fasta_db"
                    ],
                    "id": "#main/prepare_fasta_db"
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
