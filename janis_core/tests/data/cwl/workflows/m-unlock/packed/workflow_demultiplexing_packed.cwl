{
    "$graph": [
        {
            "class": "CommandLineTool",
            "label": "NGTax demultixplexing",
            "doc": "Runs NGTAX amplicon demultiplexing procedure\n",
            "requirements": [
                {
                    "class": "InlineJavascriptRequirement"
                }
            ],
            "hints": [
                {
                    "dockerPull": "docker-registry.wur.nl/m-unlock/docker/ngtax:2.2.6",
                    "class": "DockerRequirement"
                }
            ],
            "inputs": [
                {
                    "type": "string",
                    "doc": "Forward primer used",
                    "label": "The forward primer used",
                    "inputBinding": {
                        "prefix": "-for_p"
                    },
                    "id": "#ngtax_demultiplexing.cwl/forward_primer"
                },
                {
                    "type": "File",
                    "doc": "Forward library file",
                    "label": "The forward library used",
                    "id": "#ngtax_demultiplexing.cwl/forward_reads"
                },
                {
                    "type": "File",
                    "doc": "Mapping file containing barcode information",
                    "label": "The mapping file",
                    "inputBinding": {
                        "prefix": "-mapFile"
                    },
                    "id": "#ngtax_demultiplexing.cwl/mapping_file"
                },
                {
                    "type": "string",
                    "doc": "Reverse primer used",
                    "label": "The reverse primer used",
                    "inputBinding": {
                        "prefix": "-rev_p"
                    },
                    "id": "#ngtax_demultiplexing.cwl/reverse_primer"
                },
                {
                    "type": "File",
                    "doc": "Reverse library file",
                    "label": "The reverse library used",
                    "id": "#ngtax_demultiplexing.cwl/reverse_reads"
                }
            ],
            "baseCommand": [
                "-demultiplex",
                "-output",
                "demultiplexed"
            ],
            "arguments": [
                {
                    "prefix": "-fastQ",
                    "valueFrom": "$(inputs.forward_reads.path),$(inputs.reverse_reads.path)"
                }
            ],
            "id": "#ngtax_demultiplexing.cwl",
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
            "https://schema.org/license": "https://spdx.org/licenses/Apache-2.0",
            "https://schema.org/copyrightHolder": "UNLOCK - Unlocking Microbial Potential",
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputBinding": {
                        "glob": "demultiplexed/*"
                    },
                    "id": "#ngtax_demultiplexing.cwl/output"
                }
            ]
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
                    "type": "string",
                    "doc": "Forward primer used",
                    "label": "The forward primer used",
                    "id": "#main/forward_primer"
                },
                {
                    "type": "File",
                    "doc": "Forward library file",
                    "label": "The forward library used",
                    "id": "#main/forward_reads"
                },
                {
                    "type": "File",
                    "doc": "Mapping file containing barcode information",
                    "label": "The mapping file",
                    "id": "#main/mapping_file"
                },
                {
                    "type": "string",
                    "doc": "Reverse primer used",
                    "label": "The reverse primer used",
                    "id": "#main/reverse_primer"
                },
                {
                    "type": "File",
                    "doc": "Reverse library file",
                    "label": "The reverse library used",
                    "id": "#main/reverse_reads"
                }
            ],
            "steps": [
                {
                    "run": "#ngtax_demultiplexing.cwl",
                    "in": [
                        {
                            "source": "#main/forward_primer",
                            "id": "#main/ngtax/forward_primer"
                        },
                        {
                            "source": "#main/forward_reads",
                            "id": "#main/ngtax/forward_reads"
                        },
                        {
                            "source": "#main/mapping_file",
                            "id": "#main/ngtax/mapping_file"
                        },
                        {
                            "source": "#main/reverse_primer",
                            "id": "#main/ngtax/reverse_primer"
                        },
                        {
                            "source": "#main/reverse_reads",
                            "id": "#main/ngtax/reverse_reads"
                        }
                    ],
                    "out": [
                        "#main/ngtax/output"
                    ],
                    "id": "#main/ngtax"
                }
            ],
            "outputs": [
                {
                    "type": {
                        "type": "array",
                        "items": "File"
                    },
                    "outputSource": "#main/ngtax/output",
                    "id": "#main/demultiplex_output"
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
